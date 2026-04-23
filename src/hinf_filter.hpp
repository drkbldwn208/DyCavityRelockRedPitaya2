#pragma once

#include <ap_fixed.h>
#include <stdint.h>
#include "hinf_coeffs.h"

#define HINF_COEF_TOTAL_BITS (HINF_COEF_INT_BITS + HINF_COEF_FRAC_BITS)

// 1. Widen types and enforce Saturation (AP_SAT)
typedef ap_fixed<16, 1, AP_TRN, AP_SAT> sig_t;     // Q1.15 signal
typedef ap_fixed<HINF_COEF_TOTAL_BITS, HINF_COEF_INT_BITS> coeff_t; 

// The inter-stage pipe. ±2048 headroom.
typedef ap_fixed<32, 12, AP_TRN, AP_SAT> pipe_t;   

// The master accumulator. ±16 million headroom. 
// This Maps perfectly to the Xilinx DSP48 slice.
typedef ap_fixed<64, 24, AP_TRN, AP_SAT> acc_t;    

class HinfFilter {
  private:
    // DF-I requires us to remember the last two inputs (x) and outputs (y) per section
    pipe_t x_hist[HINF_N_SECTIONS][2];
    pipe_t y_hist[HINF_N_SECTIONS][2];

    static coeff_t int_to_coeff(int32_t raw) {
      #pragma HLS INLINE
      coeff_t c;
      c.range() = raw;
      return c;
    }

    // Direct Form I Biquad
    void biquad_df1(pipe_t x_in, pipe_t &y_out, int sec,
                    coeff_t b0, coeff_t b1, coeff_t b2,
                    coeff_t a1, coeff_t a2)
    {
      #pragma HLS INLINE
      
      // Compute the entire sum in one massive 64-bit register.
      // Because 'acc_t' is so wide, internal overflow during the addition is impossible.
      acc_t acc = (acc_t)(b0 * x_in) 
                + (acc_t)(b1 * x_hist[sec][0]) 
                + (acc_t)(b2 * x_hist[sec][1]) 
                - (acc_t)(a1 * y_hist[sec][0]) 
                - (acc_t)(a2 * y_hist[sec][1]);

      // Saturate and truncate down to the inter-stage pipe width
      y_out = pipe_t(acc);

      // Shift the history registers for the next clock cycle
      x_hist[sec][1] = x_hist[sec][0];
      x_hist[sec][0] = x_in;
      
      y_hist[sec][1] = y_hist[sec][0];
      y_hist[sec][0] = y_out;
    }

  public:
    HinfFilter() {
      // Completely partition the history arrays so they become individual flip-flops
      #pragma HLS ARRAY_PARTITION variable=x_hist complete dim=0
      #pragma HLS ARRAY_PARTITION variable=y_hist complete dim=0

      for (int i = 0; i < HINF_N_SECTIONS; i++) {
        #pragma HLS UNROLL
        x_hist[i][0] = 0; x_hist[i][1] = 0;
        y_hist[i][0] = 0; y_hist[i][1] = 0;
      }
    }

    sig_t process(sig_t u_in) {
      #pragma HLS INLINE

      pipe_t pipe[HINF_N_SECTIONS + 1];
      #pragma HLS ARRAY_PARTITION variable=pipe complete
      pipe[0] = pipe_t(u_in);  

      CASCADE:
      for (int i = 0; i < HINF_N_SECTIONS; i++) {
        #pragma HLS UNROLL
        coeff_t cb0 = int_to_coeff(HINF_SOS[i].b0);
        coeff_t cb1 = int_to_coeff(HINF_SOS[i].b1);
        coeff_t cb2 = int_to_coeff(HINF_SOS[i].b2);
        coeff_t ca1 = int_to_coeff(HINF_SOS[i].a1);
        coeff_t ca2 = int_to_coeff(HINF_SOS[i].a2);

        // Call the new DF-I structure
        biquad_df1(pipe[i], pipe[i + 1], i, cb0, cb1, cb2, ca1, ca2);
      }
      
      // Final output saturation down to Q1.15
      return sig_t(pipe[HINF_N_SECTIONS]);  
    }
};