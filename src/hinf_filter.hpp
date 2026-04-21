#pragma once

#include <ap_fixed.h>
#include <stdint.h>
#include "hinf_coeffs.h"

#define HINF_COEF_TOTAL_BITS (HINF_COEF_INT_BITS + HINF_COEF_FRAC_BITS)

typedef ap_fixed<16, 1> sig_t; // Q1.15 signal
typedef ap_fixed<HINF_COEF_TOTAL_BITS, HINF_COEF_INT_BITS> coeff_t; 
typedef ap_fixed<40, 8> state_t; // Q8.32 state
typedef ap_fixed<64, 12> acc_t; // Q12.52 accumulator

class HinfFilter {
  private:
    state_t w1[HINF_N_SECTIONS];
    state_t w2[HINF_N_SECTIONS];

    static coeff_t int_to_coeff(int32_t raw) {
      #pragma HLS INLINE
      coeff_t c;
      c.range() = raw;
      return c;
    }

    static void biquad(sig_t x, sig_t &y,
                    state_t &w1, state_t &w2,
                    coeff_t b0, coeff_t b1, coeff_t b2,
                    coeff_t a1, coeff_t a2)
    {
    #pragma HLS INLINE
    acc_t y_acc = acc_t(b0 * x) + (acc_t)w1;
    sig_t y_out = sig_t(y_acc);

    acc_t w1_new = acc_t(b1 * x) - acc_t(a1 * y_out) + (acc_t)w2;
    acc_t w2_new = acc_t(b2 * x) - acc_t(a2 * y_out);

    w1 = state_t(w1_new);
    w2 = state_t(w2_new);
    y = y_out;
    }

  public:
    HinfFilter() {
      #pragma HLS ARRAY_PARTITION variable = w1 complete
      #pragma HLS ARRAY_PARTITION variable = w2 complete

      #pragma HLS DEPENDENCE variable = w1 type=inter direction=RAW distance=4
      #pragma HLS DEPENDENCE variable = w2 type=inter direction=RAW distance=4

      for (int i = 0; i < HINF_N_SECTIONS; i++) {
        #pragma HLS UNROLL
        w1[i] = 0;
        w2[i] = 0;
      }
    }

    sig_t process(sig_t u_in) {
      #pragma HLS INLINE

      sig_t pipe[HINF_N_SECTIONS + 1];
      #pragma HLS ARRAY_PARTITION variable = pipe complete
      pipe[0] = u_in;

      CASCADE: 
        for (int i = 0; i < HINF_N_SECTIONS; i++) {
          #pragma HLS UNROLL
          coeff_t cb0 = int_to_coeff(HINF_SOS[i].b0);
          coeff_t cb1 = int_to_coeff(HINF_SOS[i].b1);
          coeff_t cb2 = int_to_coeff(HINF_SOS[i].b2);
          coeff_t ca1 = int_to_coeff(HINF_SOS[i].a1);
          coeff_t ca2 = int_to_coeff(HINF_SOS[i].a2);

          biquad(pipe[i], pipe[i + 1], w1[i], w2[i], cb0, cb1, cb2, ca1, ca2);
        }
      return pipe[HINF_N_SECTIONS];
    }
};