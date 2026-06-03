"""
Quantize the discrete controller (K_zpk.npz) to Q3.29 SOS and write
src/hinf_coeffs.h if the quantization error is small enough.

  python3 scripts/coeffs_analyze.py             # defaults
  python3 scripts/coeffs_analyze.py --no-show

Reports four error budgets:
  1. SOS decomposition vs. ZPK truth  (≈ float64 machine epsilon)
  2. Quantized SOS vs. ZPK truth       (total error the FPGA will see)
  3. Quantized SOS vs. float SOS       (isolated quantization effect)
  4. Fixed-point arithmetic step test  (HLS realization vs. quantized float SOS)

The fourth check matters for controllers with low-frequency pole/zero
cancellations.  Those cancellations can look harmless in the quantized
frequency response while still producing DC drift or limit-cycle artifacts in
the actual fixed-point direct-form realization.
"""

import argparse
import numpy as np
import scipy.signal as sig
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--zpk",       default="K_zpk.npz")
    p.add_argument("--header",    default="src/hinf_coeffs.h")
    p.add_argument("--total-bits", type=int, default=32)
    p.add_argument("--int-bits",  type=int, default=3)
    p.add_argument("--max-qerr-db", type=float, default=-10.0,
                   help="don't write header if pure-quant error exceeds this (dB)")
    p.add_argument("--fixedpoint-step-counts", type=int, default=50,
                   help="signed ADC-count step used for fixed-point arithmetic verification")
    p.add_argument("--fixedpoint-samples", type=int, default=20000,
                   help="number of controller samples in fixed-point arithmetic verification")
    p.add_argument("--max-fixedpoint-error-counts", type=float, default=500.0,
                   help="don't write header if fixed-point step differs from quantized float by more than this many output counts")
    p.add_argument("--no-fixedpoint-check", action="store_true",
                   help="skip fixed-point arithmetic verification")
    p.add_argument("--no-show", action="store_true")
    return p.parse_args()


SIG_FRAC = 15
PIPE_FRAC = 20
ACC_FRAC = 40
SIG_BITS = 16
PIPE_BITS = 32


def sat_signed(x, bits):
    return max(-(1 << (bits - 1)), min((1 << (bits - 1)) - 1, int(x)))


def qshift(value, shift, mode):
    """Quantize integer raw fixed-point value by 2**shift."""
    value = int(value)
    if shift <= 0:
        return value << (-shift)

    scale = 1 << shift
    if mode == "trn_zero":
        return (abs(value) // scale) * (1 if value >= 0 else -1)
    if mode == "trn":
        return value // scale
    if mode == "rnd_conv":
        q = value // scale
        r = value - q * scale
        twice = 2 * r
        if twice < scale:
            return q
        if twice > scale:
            return q + 1
        return q if (q & 1) == 0 else q + 1
    raise ValueError(f"unknown quantization mode {mode}")


def run_fixedpoint_df1(sos_int, coef_frac_bits, x_sig, qmode="trn_zero"):
    """Approximate the HLS DF-I cascade in raw sig_t counts.

    This mirrors hinf_filter.hpp after the AP_TRN_ZERO change:
      sig_t  Q1.15, saturating
      pipe_t Q12.20, saturating
      acc_t  Q24.40, wrapping in HLS; no practical wrap expected here

    The model quantizes each product into acc_t before summing, then quantizes
    acc_t into pipe_t, matching the two places where the HLS realization loses
    fractional bits.
    """
    n_sec = len(sos_int)
    x_hist = np.zeros((n_sec, 2), dtype=object)
    y_hist = np.zeros((n_sec, 2), dtype=object)
    y_sig = np.zeros_like(x_sig, dtype=np.int64)

    product_shift = coef_frac_bits + PIPE_FRAC - ACC_FRAC
    pipe_shift = ACC_FRAC - PIPE_FRAC
    sig_shift = PIPE_FRAC - SIG_FRAC

    for n, x in enumerate(x_sig):
        pipe = int(x) << sig_shift
        for i, (b0, b1, b2, a1, a2) in enumerate(sos_int):
            terms = (
                b0 * pipe,
                b1 * int(x_hist[i, 0]),
                b2 * int(x_hist[i, 1]),
                -a1 * int(y_hist[i, 0]),
                -a2 * int(y_hist[i, 1]),
            )
            acc = sum(qshift(term, product_shift, qmode) for term in terms)
            pipe_out = sat_signed(qshift(acc, pipe_shift, qmode), PIPE_BITS)

            x_hist[i, 1] = x_hist[i, 0]
            x_hist[i, 0] = pipe
            y_hist[i, 1] = y_hist[i, 0]
            y_hist[i, 0] = pipe_out
            pipe = pipe_out

        y_sig[n] = sat_signed(qshift(pipe, sig_shift, qmode), SIG_BITS)
    return y_sig


def run_float_sos_counts(sos_int, scale, x_sig):
    sos = np.array([[b0 / scale, b1 / scale, b2 / scale, 1.0, a1 / scale, a2 / scale]
                    for (b0, b1, b2, a1, a2) in sos_int])
    return sig.sosfilt(sos, x_sig.astype(float))


def main():
    args = parse_args()
    FRAC_BITS = args.total_bits - args.int_bits
    SCALE     = 2 ** FRAC_BITS
    COEF_MAX  =  2 ** (args.total_bits - 1) - 1
    COEF_MIN  = -2 ** (args.total_bits - 1)

    d = np.load(args.zpk)
    z_d, p_d, k_d = d["z"], d["p"], d["k"]
    fs = float(d["fs"])
    
    # --- ADDED: Explicitly print exact continuous-equivalent locations in Hz ---
    print("\n--- Exact Pole/Zero Locations (Hz) ---")
    print("Poles:")
    for i, p in enumerate(p_d):
        f_hz = np.log(complex(p)) * fs / (2 * np.pi)
        print(f"  p_{i}: {f_hz.real:+.3f} {f_hz.imag:+.3f}j Hz (mag: {abs(f_hz):.3f} Hz)")
        
    print("Zeros:")
    for i, z in enumerate(z_d):
        f_hz = np.log(complex(z)) * fs / (2 * np.pi)
        print(f"  z_{i}: {f_hz.real:+.3f} {f_hz.imag:+.3f}j Hz (mag: {abs(f_hz):.3f} Hz)")
    print("--------------------------------------\n")

    sos = sig.zpk2sos(z_d, p_d, k_d, pairing='nearest')
    n_sec = sos.shape[0]
    print(f"Controller: {len(p_d)} poles  fs={fs/1e6:.4f} MHz  → {n_sec} SOS sections")

    for i, s in enumerate(sos):
        fp = sorted(abs(np.angle(q)) * fs / (2*np.pi)
                    for q in np.roots([1.0, s[4], s[5]]))
        print(f"  Sec {i}: pole freqs ≈ {fp[0]:.0f} / {fp[1]:.0f} Hz")

    max_coef = max(abs(c) for s in sos for c in (s[0], s[1], s[2], s[4], s[5])
                   if abs(c) > 1e-15)
    headroom = 2**(args.int_bits - 1) / max_coef
    print(f"Q{args.int_bits}.{FRAC_BITS}  max |coef|={max_coef:.3f}  headroom={headroom:.2f}×")
    if headroom < 1.0:
        print("  WARNING: a coefficient exceeds the integer range — raise --int-bits")

    def Q(x):
        return max(COEF_MIN, min(COEF_MAX, int(round(x * SCALE))))
    q_sos = [(Q(s[0]), Q(s[1]), Q(s[2]), Q(s[4]), Q(s[5])) for s in sos]

    # Three frequency responses for error budgeting
    w_hz = np.logspace(1, np.log10(fs * 0.49), 3000)
    z_e  = np.exp(1j * w_hz / fs * 2*np.pi)

    H_truth = np.full(len(z_e), k_d, dtype=complex)
    for z in z_d: H_truth *= (z_e - z)
    for p in p_d: H_truth /= (z_e - p)

    H_sos = np.ones(len(z_e), dtype=complex)
    for s in sos:
        H_sos *= np.polyval([s[0], s[1], s[2]], z_e) / np.polyval([1.0, s[4], s[5]], z_e)

    H_quant = np.ones(len(z_e), dtype=complex)
    H_secs  = []
    for b0, b1, b2, a1, a2 in q_sos:
        H_s = (np.polyval([b0, b1, b2], z_e) / SCALE) / \
              np.polyval([1.0, a1/SCALE, a2/SCALE], z_e)
        H_quant *= H_s
        H_secs.append(H_s)

    mask = w_hz < fs * 0.4
    def db_err(a, b):
        return 20 * np.log10(np.max(np.abs(a - b) / (np.abs(b) + 1e-30))[mask.nonzero()] + 1e-30)
    err_sos  = np.abs(H_sos   - H_truth) / (np.abs(H_truth) + 1e-30)
    err_quan = np.abs(H_quant - H_truth) / (np.abs(H_truth) + 1e-30)
    err_q    = np.abs(H_quant - H_sos)   / (np.abs(H_sos)   + 1e-30)
    max_qerr_db = 20*np.log10(np.max(err_q[mask]) + 1e-30)

    print(f"Errors (f < 0.4·fs):")
    print(f"  SOS vs truth:       {20*np.log10(np.max(err_sos[mask])+1e-30):+.1f} dB")
    print(f"  Quantized vs truth: {20*np.log10(np.max(err_quan[mask])+1e-30):+.1f} dB")
    print(f"  Pure quantization:  {max_qerr_db:+.1f} dB")

    print("\nLow-frequency section sums:")
    fragile_sections = []
    for i, (b0, b1, b2, a1, a2) in enumerate(q_sos):
        b_sum = b0 + b1 + b2
        den_sum = SCALE + a1 + a2
        if abs(den_sum) <= 4096 or b_sum == 0:
            fragile_sections.append(i)
        print(f"  Sec {i}: b0+b1+b2={b_sum:+d},  1+a1+a2={den_sum:+d} raw")
    if fragile_sections:
        print("  NOTE: sections", fragile_sections,
              "have exact/near DC cancellation or a very small DC denominator.")
        print("        These are physically plausible for an integrator/leaky-integrator model,")
        print("        but they are fragile in fixed-point direct-form arithmetic.")

    fixedpoint_failed = False
    max_fixedpoint_error = 0.0
    fixedpoint_summary = "skipped"
    if not args.no_fixedpoint_check:
        print("\nFixed-point arithmetic check (matches hinf_filter.hpp AP_TRN_ZERO):")
        step = int(args.fixedpoint_step_counts)
        fixedpoint_summary_rows = []
        for signed_step in (step, -step):
            x_step = np.full(args.fixedpoint_samples, signed_step, dtype=np.int64)
            y_fx = run_fixedpoint_df1(q_sos, FRAC_BITS, x_step, qmode="trn_zero")
            y_fp = run_float_sos_counts(q_sos, SCALE, x_step)
            err = np.abs(y_fx.astype(float) - y_fp)
            max_err = float(np.max(err))
            max_fixedpoint_error = max(max_fixedpoint_error, max_err)
            row = (signed_step, int(y_fx[-1]), float(y_fp[-1]), max_err)
            fixedpoint_summary_rows.append(row)
            print(f"  step {signed_step:+d} counts: fixed final={row[1]:+d}, "
                  f"float final={row[2]:+.3f}, max |err|={max_err:.3f} counts")

        fixedpoint_summary = "; ".join(
            f"step {s:+d}: fx {fx:+d}, float {fp:+.2f}, maxerr {err:.1f}"
            for s, fx, fp, err in fixedpoint_summary_rows)
        if max_fixedpoint_error > args.max_fixedpoint_error_counts:
            fixedpoint_failed = True
            print(f"  FAIL: fixed-point arithmetic error {max_fixedpoint_error:.1f} counts "
                  f"> {args.max_fixedpoint_error_counts:.1f} count limit")
        else:
            print(f"  PASS: max fixed-point arithmetic error {max_fixedpoint_error:.1f} counts "
                  f"≤ {args.max_fixedpoint_error_counts:.1f} count limit")

    # Plots
    fig, axs = plt.subplots(2, 2, figsize=(13, 9))
    fig.suptitle(f"SOS verification  |  {n_sec} sections  Q{args.int_bits}.{FRAC_BITS}")
    axs[0, 0].semilogx(w_hz, 20*np.log10(np.abs(H_truth)+1e-30), "C0", lw=2.5, label="Truth")
    axs[0, 0].semilogx(w_hz, 20*np.log10(np.abs(H_sos)+1e-30),   "C2", lw=1.2, label="SOS float")
    axs[0, 0].semilogx(w_hz, 20*np.log10(np.abs(H_quant)+1e-30), "C3", lw=1.0, ls="--", label="SOS quantized")
    axs[0, 0].set(title="Magnitude", xlabel="Hz", ylabel="dB"); axs[0, 0].grid(True, which="both", alpha=0.3); axs[0, 0].legend(fontsize=8)

    axs[0, 1].semilogx(w_hz, np.degrees(np.unwrap(np.angle(H_truth))), "C0", lw=2.5, label="Truth")
    axs[0, 1].semilogx(w_hz, np.degrees(np.unwrap(np.angle(H_sos))),   "C2", lw=1.2, label="SOS float")
    axs[0, 1].semilogx(w_hz, np.degrees(np.unwrap(np.angle(H_quant))), "C3", lw=1.0, ls="--", label="SOS quantized")
    axs[0, 1].set(title="Phase", xlabel="Hz", ylabel="degrees"); axs[0, 1].grid(True, which="both", alpha=0.3); axs[0, 1].legend(fontsize=8)

    axs[1, 0].semilogx(w_hz, 20*np.log10(err_sos+1e-30),  "C2", label="SOS vs truth")
    axs[1, 0].semilogx(w_hz, 20*np.log10(err_quan+1e-30), "C3", label="Quantized vs truth")
    axs[1, 0].semilogx(w_hz, 20*np.log10(err_q+1e-30),    "C1", ls="--", label="Pure quantization")
    axs[1, 0].set(title="Error budget", xlabel="Hz", ylabel="dB"); axs[1, 0].set_ylim([-160, 10])
    axs[1, 0].grid(True, which="both", alpha=0.3); axs[1, 0].legend(fontsize=8)

    for i, H_s in enumerate(H_secs):
        fp = sorted(abs(np.angle(q)) * fs / (2*np.pi) for q in np.roots([1.0, sos[i, 4], sos[i, 5]]))
        axs[1, 1].semilogx(w_hz, 20*np.log10(np.abs(H_s)+1e-30), lw=1.5, label=f"Sec {i} ({fp[0]:.0f}/{fp[1]:.0f} Hz)")
    axs[1, 1].set(title="Per-section", xlabel="Hz", ylabel="dB"); axs[1, 1].grid(True, which="both", alpha=0.3); axs[1, 1].legend(fontsize=8)

    plt.tight_layout()
    plt.savefig("sos_verification.png", dpi=120)
    print("Saved sos_verification.png")

    # Write header
    if max_qerr_db > args.max_qerr_db:
        print(f"Pure quant error {max_qerr_db:.1f} dB > {args.max_qerr_db} dB — header NOT written")
    elif fixedpoint_failed:
        print(f"Fixed-point arithmetic check failed — header NOT written")
    else:
        lines = [
            "// hinf_coeffs.h — auto-generated by coeffs_analyze.py",
            f"// Q{args.int_bits}.{FRAC_BITS}  {n_sec} sections  fs={fs/1e6:.4f} MHz  "
            f"quant err {max_qerr_db:.1f} dB",
            f"// HLS arithmetic check: AP_TRN_ZERO, {fixedpoint_summary}",
            "#ifndef HINF_COEFFS_H",
            "#define HINF_COEFFS_H",
            "#include <stdint.h>",
            "",
            f"#define HINF_N_SECTIONS     {n_sec}",
            f"#define HINF_COEF_INT_BITS  {args.int_bits}",
            f"#define HINF_COEF_FRAC_BITS {FRAC_BITS}",
            f"#define HINF_COEF_SCALE     {SCALE}L",
            "",
            "typedef struct { int32_t b0,b1,b2,a1,a2; } hinf_sos_t;",
            "",
            f"static const hinf_sos_t HINF_SOS[{n_sec}] = {{",
        ]
        for i, (b0, b1, b2, a1, a2) in enumerate(q_sos):
            comma = "," if i < n_sec - 1 else ""
            lines.append(f"    {{{b0:11d},{b1:11d},{b2:11d},{a1:11d},{a2:11d}}}{comma}")
        lines += ["};", "", "#endif", ""]
        with open(args.header, "w") as f:
            f.write("\n".join(lines))
        print(f"Wrote {args.header}")

    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
