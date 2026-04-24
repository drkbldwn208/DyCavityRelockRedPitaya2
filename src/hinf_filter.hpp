#pragma once

#include <ap_fixed.h>
#include <stdint.h>
#include "hinf_coeffs.h"

// =============================================================================
// Direct Form I biquad cascade — fixed-point, saturating.
//
// Fixed-point formats:
//   sig_t   Q1.15  (ADC/DAC sample)
//   pipe_t  Q12.20 (inter-stage signal; real headroom ±2048)
//   acc_t   Q24.40 (per-biquad accumulator; real headroom ±16 million)
//   coeff_t Q3.29  (loaded from hinf_coeffs.h)
//
// On the DF-I overflow-cancellation property:
//   DF-I has a known property that transient overflows of the internal
//   accumulator cancel algebraically, PROVIDED the accumulator uses modulo
//   (wrap) arithmetic AND the biquad's final output is actually representable
//   in its output format. That property saves you when you want a tight
//   accumulator — it does NOT let you shrink the INTER-STAGE pipe.
//
//   Here we don't rely on it: acc_t has ~24 bits of real-value headroom
//   (millions of ×), so arithmetic overflow of the accumulator is impossible
//   for any reasonable signal. AP_SAT on acc_t is defensive only; it never
//   fires.
//
//   AP_SAT on pipe_t is the MEANINGFUL protection. If the true filter output
//   between biquads exceeded ±2048, AP_WRAP there would fold the signal and
//   destroy the downstream biquads' state irreversibly; AP_SAT instead clips
//   gracefully. For the current controller the peak pipe value is ~1.5e-3
//   out of ±2048 real — >1e6× headroom, so pipe_t never saturates either.
// =============================================================================

#define HINF_COEF_TOTAL_BITS (HINF_COEF_INT_BITS + HINF_COEF_FRAC_BITS)

typedef ap_fixed<16, 1,  AP_TRN, AP_SAT> sig_t;
typedef ap_fixed<32, 12, AP_TRN, AP_SAT> pipe_t;
typedef ap_fixed<64, 24, AP_TRN, AP_SAT> acc_t;
typedef ap_fixed<HINF_COEF_TOTAL_BITS, HINF_COEF_INT_BITS> coeff_t;


class HinfFilter {
  private:
    pipe_t x_hist[HINF_N_SECTIONS][2];
    pipe_t y_hist[HINF_N_SECTIONS][2];

    static coeff_t int_to_coeff(int32_t raw) {
      #pragma HLS INLINE
      coeff_t c;
      c.range() = raw;
      return c;
    }

    void biquad_df1(pipe_t x_in, pipe_t &y_out, int sec,
                    coeff_t b0, coeff_t b1, coeff_t b2,
                    coeff_t a1, coeff_t a2)
    {
      #pragma HLS INLINE
      acc_t acc = (acc_t)(b0 * x_in)
                + (acc_t)(b1 * x_hist[sec][0])
                + (acc_t)(b2 * x_hist[sec][1])
                - (acc_t)(a1 * y_hist[sec][0])
                - (acc_t)(a2 * y_hist[sec][1]);

      y_out = pipe_t(acc);                // Q24.40 -> Q12.20, saturating

      x_hist[sec][1] = x_hist[sec][0];
      x_hist[sec][0] = x_in;
      y_hist[sec][1] = y_hist[sec][0];
      y_hist[sec][0] = y_out;
    }

  public:
    HinfFilter() {
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

      pipe[0] = pipe_t(u_in);             // Q1.15 -> Q12.20

      CASCADE:
      for (int i = 0; i < HINF_N_SECTIONS; i++) {
        #pragma HLS UNROLL
        biquad_df1(pipe[i], pipe[i + 1], i,
                   int_to_coeff(HINF_SOS[i].b0),
                   int_to_coeff(HINF_SOS[i].b1),
                   int_to_coeff(HINF_SOS[i].b2),
                   int_to_coeff(HINF_SOS[i].a1),
                   int_to_coeff(HINF_SOS[i].a2));
      }
      return sig_t(pipe[HINF_N_SECTIONS]); // Q12.20 -> Q1.15, saturating
    }
};
