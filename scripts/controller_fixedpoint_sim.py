"""
Bit-accurate Python model of the HLS Direct-Form-I biquad cascade.

Mirrors hinf_filter.hpp exactly:
  sig_t   Q1.15   (ADC/DAC words)
  pipe_t  Q12.20  (inter-stage; ±2048 real headroom)
  acc_t   Q24.40  (accumulator; ±16M headroom)

The controller operates on pipe-scale integers internally; sig_t values are
shifted up by PIPE_SHIFT = 20-15 = 5 on the way in, and shifted down +
saturated to ±32767 on the way out.

Usage:
  python3 scripts/controller_fixedpoint_sim.py       # step response
"""

import os
import re
import numpy as np
import scipy.signal as sig

PIPE_SHIFT = 5                     # Q12.20 - Q1.15
SIG_MAX,  SIG_MIN  = 32767, -32768


def parse_hinf_coeffs(path):
    with open(path) as f:
        text = f.read()
    def grab(name):
        m = re.search(rf"#define\s+{name}\s+(-?\d+)L?", text)
        if not m:
            raise RuntimeError(f"{name} not found in {path}")
        return int(m.group(1))
    n_sec    = grab("HINF_N_SECTIONS")
    int_bits = grab("HINF_COEF_INT_BITS")
    frac_bits= grab("HINF_COEF_FRAC_BITS")
    scale    = grab("HINF_COEF_SCALE")
    rows = re.findall(
        r"\{\s*(-?\d+)\s*,\s*(-?\d+)\s*,\s*(-?\d+)\s*,\s*(-?\d+)\s*,\s*(-?\d+)\s*\}",
        text.split("HINF_SOS")[-1])
    sos = [tuple(int(v) for v in r) for r in rows[:n_sec]]
    return sos, scale, int_bits, frac_bits


def run_fixedpoint(sos_int, scale, x_sig, saturate_output=True):
    """Input and output are sig_t (Q1.15) raw integers, matching HinfFilter::process()."""
    n_sec  = len(sos_int)
    x_hist = np.zeros((n_sec, 2), dtype=np.int64)
    y_hist = np.zeros((n_sec, 2), dtype=np.int64)
    y_sig  = np.zeros_like(x_sig, dtype=np.int64)

    for n in range(len(x_sig)):
        pipe = int(x_sig[n]) << PIPE_SHIFT
        for i, (b0, b1, b2, a1, a2) in enumerate(sos_int):
            acc = (b0 * pipe
                   + b1 * int(x_hist[i, 0])
                   + b2 * int(x_hist[i, 1])
                   - a1 * int(y_hist[i, 0])
                   - a2 * int(y_hist[i, 1]))
            pipe_out = acc // scale           # AP_TRN: floor toward -inf
            x_hist[i, 1] = x_hist[i, 0]
            x_hist[i, 0] = pipe
            y_hist[i, 1] = y_hist[i, 0]
            y_hist[i, 0] = pipe_out
            pipe = int(pipe_out)
        out = pipe >> PIPE_SHIFT
        if saturate_output:
            out = max(SIG_MIN, min(SIG_MAX, out))
        y_sig[n] = out
    return y_sig


def run_float_reference(sos_int, scale, x_sig):
    """Same coefficients, float64 SOS, sig_t-scale inputs/outputs."""
    sos = np.array([[b0/scale, b1/scale, b2/scale, 1.0, a1/scale, a2/scale]
                    for (b0, b1, b2, a1, a2) in sos_int])
    return sig.sosfilt(sos, x_sig.astype(float))


def main():
    path = os.path.join(os.path.dirname(__file__), "..", "src", "hinf_coeffs.h")
    sos_int, scale, ib, fb = parse_hinf_coeffs(os.path.normpath(path))
    print(f"Loaded {len(sos_int)} SOS sections from hinf_coeffs.h  (Q{ib}.{fb})")

    x = np.full(2000, 50, dtype=np.int64)     # step in sig_t counts
    y_fx = run_fixedpoint(sos_int, scale, x)
    y_fp = run_float_reference(sos_int, scale, x)

    print("\nStep response (sig_t counts):")
    print(f"  fixed-point: last value = {y_fx[-1]}")
    print(f"  float64 ref: last value = {y_fp[-1]:.3f}")
    print(f"  max |fx - fp| = {np.max(np.abs(y_fx.astype(float) - y_fp)):.3f}")


if __name__ == "__main__":
    main()
