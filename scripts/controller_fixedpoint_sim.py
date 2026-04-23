"""
Bit-accurate Python model of the HLS Direct-Form-I biquad cascade.

Parses src/hinf_coeffs.h directly so it is always in sync with the hardware.
Compares against a Float-64 SOS cascade to report the truncation error.
"""

import os
import re
import sys
import numpy as np


def parse_hinf_coeffs(path):
    with open(path) as f:
        text = f.read()
    def grab_int(name):
        m = re.search(rf"#define\s+{name}\s+(-?\d+)", text)
        if not m:
            raise RuntimeError(f"{name} not found in {path}")
        return int(m.group(1))
    def grab_scale(name):
        m = re.search(rf"#define\s+{name}\s+(-?\d+)L?", text)
        return int(m.group(1))

    n_sec      = grab_int("HINF_N_SECTIONS")
    int_bits   = grab_int("HINF_COEF_INT_BITS")
    frac_bits  = grab_int("HINF_COEF_FRAC_BITS")
    scale      = grab_scale("HINF_COEF_SCALE")

    sos_ints = re.findall(
        r"\{\s*(-?\d+)\s*,\s*(-?\d+)\s*,\s*(-?\d+)\s*,\s*(-?\d+)\s*,\s*(-?\d+)\s*\}",
        text.split("HINF_SOS")[-1]
    )
    sos = [tuple(int(v) for v in row) for row in sos_ints[:n_sec]]
    return sos, scale, int_bits, frac_bits


def run_fixedpoint(sos_int, scale, x_in):
    """Bit-accurate DF-I cascade. Each biquad's accumulator is int64."""
    n_sec = len(sos_int)
    x_hist = np.zeros((n_sec, 2), dtype=np.int64)
    y_hist = np.zeros((n_sec, 2), dtype=np.int64)
    y_out = np.zeros_like(x_in, dtype=np.int64)

    for n in range(len(x_in)):
        pipe = int(x_in[n])
        for i, (b0, b1, b2, a1, a2) in enumerate(sos_int):
            acc = (b0 * pipe
                   + b1 * int(x_hist[i, 0])
                   + b2 * int(x_hist[i, 1])
                   - a1 * int(y_hist[i, 0])
                   - a2 * int(y_hist[i, 1]))
            pipe_out = acc // scale
            x_hist[i, 1] = x_hist[i, 0]
            x_hist[i, 0] = pipe
            y_hist[i, 1] = y_hist[i, 0]
            y_hist[i, 0] = pipe_out
            pipe = int(pipe_out)
        y_out[n] = pipe
    return y_out


def run_float(sos_int, scale, x_in):
    """Float-64 SOS reference using the same coefficients (no quantization effect)."""
    import scipy.signal as sig
    sos = np.array([[b0/scale, b1/scale, b2/scale, 1.0, a1/scale, a2/scale]
                    for (b0, b1, b2, a1, a2) in sos_int])
    return sig.sosfilt(sos, x_in)


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    coeffs_path = os.path.normpath(os.path.join(here, "..", "src", "hinf_coeffs.h"))
    sos_int, scale, ib, fb = parse_hinf_coeffs(coeffs_path)
    print(f"Loaded {len(sos_int)} SOS sections from {coeffs_path}  (Q{ib}.{fb})")

    x = np.full(200, 50.0)
    y_int   = run_fixedpoint(sos_int, scale, x)
    y_float = run_float(sos_int, scale, x)

    print("\nStep response (last 10 samples):")
    print("  n   fixed-point   float")
    for i in range(190, 200):
        print(f"  {i:3d}   {int(y_int[i]):10d}   {y_float[i]:10.3f}")

    diff = y_int.astype(float) - y_float
    print(f"\nMax |fixed-point - float| = {np.max(np.abs(diff)):.3f}")


if __name__ == "__main__":
    main()
