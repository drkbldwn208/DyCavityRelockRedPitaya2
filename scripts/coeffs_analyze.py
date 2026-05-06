"""
Quantize the discrete controller (K_zpk.npz) to Q3.29 SOS and write
src/hinf_coeffs.h if the quantization error is small enough.

  python3 scripts/coeffs_analyze.py             # defaults
  python3 scripts/coeffs_analyze.py --no-show

Reports three error budgets:
  1. SOS decomposition vs. ZPK truth  (≈ float64 machine epsilon)
  2. Quantized SOS vs. ZPK truth       (total error the FPGA will see)
  3. Quantized SOS vs. float SOS       (isolated quantization effect)
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
    p.add_argument("--no-show", action="store_true")
    return p.parse_args()


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
    else:
        lines = [
            "// hinf_coeffs.h — auto-generated by coeffs_analyze.py",
            f"// Q{args.int_bits}.{FRAC_BITS}  {n_sec} sections  fs={fs/1e6:.4f} MHz  "
            f"quant err {max_qerr_db:.1f} dB",
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