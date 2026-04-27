"""
Fit a measured open-loop Bode plot with a rational IIR model.

Usage:
    python3 scripts/bode_fit.py <data.csv> [--poles N] [--zeros M] [--out plant.npz]

    data.csv must have columns:  freq_hz, mag_db, phase_deg

Preferred backend is IIRrational (pip install iirrational) if installed;
otherwise falls back to an inline Levi least-squares fit with one
Sanathanan-Koerner re-weighting pass.

Outputs:
    bode_fit.png  — measurement vs fitted model
    <out>.npz     — ZPK: z, p, k, fs (fs=0 marks continuous-time)
"""

import argparse
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig


def load_bode_csv(path):
    data = np.genfromtxt(path, delimiter=',', names=True, skip_header=0)
    try:
        f_hz      = np.asarray(data["freq_hz"],   dtype=float)
        mag_db    = np.asarray(data["mag_db"],    dtype=float)
        phase_deg = np.asarray(data["phase_deg"], dtype=float)
    except ValueError:
        # fall back: positional columns freq_hz, mag_db, phase_deg
        arr = np.loadtxt(path, delimiter=',', skiprows=1)
        f_hz, mag_db, phase_deg = arr[:, 0], arr[:, 1], arr[:, 2]
    mag   = 10 ** (mag_db / 20.0)
    phase = np.radians(phase_deg)
    return f_hz, mag * np.exp(1j * phase)

def fit_levi_mag(f_hz, H, n_poles, n_zeros, n_iter=8):
    """
    Magnitude-only rational fit via iterated Levi.

    Measured |H| is fixed throughout; the phase target is bootstrapped from
    the current model each iteration so the measured phase is never used.
    Typically converges in 3-5 iterations.
    """
    mag = np.abs(H)
    H_target = mag.astype(complex)          # seed: zero phase
    z, p, k = fit_levi(f_hz, H_target, n_poles, n_zeros)
    for _ in range(n_iter - 1):
        phase_model = np.angle(eval_zpk(z, p, k, f_hz))
        H_target = mag * np.exp(1j * phase_model)
        z, p, k = fit_levi(f_hz, H_target, n_poles, n_zeros)
    return z, p, k

def fit_iirrational(f_hz, H, order):
    import IIRrational.v2 as iirr
    fit = iirr.data2filter(data=H, F_Hz=f_hz)
    fit.choose(order=order)
    z = np.asarray(fit.zeros,  dtype=complex)
    p = np.asarray(fit.poles,  dtype=complex)
    k = float(fit.gain)
    return z, p, k


def fit_levi(f_hz, H, n_poles, n_zeros):
    """Levi least-squares rational fit in the s-domain.

    Fits H(jω) ≈ (b0 + b1·s + … + bm·s^m) / (1 + a1·s + … + an·s^n)
    by solving a single linear least-squares problem. One Sanathanan-Koerner
    iteration is applied to de-bias the fit from the implicit |D|² weighting.
    """
    w = 2 * np.pi * f_hz
    s = 1j * w

    def solve(weight):
        # Columns: [s^0 .. s^m, -H·s^1 .. -H·s^n]; target: H
        rows_re, rows_im, b_re, b_im = [], [], [], []
        for k in range(len(w)):
            row = ([s[k]**i for i in range(n_zeros + 1)] +
                   [-H[k] * s[k]**i for i in range(1, n_poles + 1)])
            rows_re.append([np.real(c) * weight[k] for c in row])
            rows_im.append([np.imag(c) * weight[k] for c in row])
            b_re.append(np.real(H[k]) * weight[k])
            b_im.append(np.imag(H[k]) * weight[k])
        A = np.vstack([np.array(rows_re), np.array(rows_im)])
        b = np.concatenate([b_re, b_im])
        x, *_ = np.linalg.lstsq(A, b, rcond=None)
        num = x[:n_zeros + 1][::-1]                     # highest power first
        den = np.concatenate(([1.0], x[n_zeros + 1:]))[::-1]
        return np.real(num), np.real(den)

    weight = np.ones_like(w)
    num, den = solve(weight)
    # Sanathanan-Koerner re-weighting: 1/|D(jω)|
    D = np.polyval(den, s)
    weight = 1.0 / np.maximum(np.abs(D), 1e-12)
    num, den = solve(weight)

    z, p, k = sig.tf2zpk(num, den)
    return z, p, k


def eval_zpk(z, p, k, f_hz):
    s = 1j * 2 * np.pi * f_hz
    num = np.ones_like(s, dtype=complex) * k
    for zi in z: num *= (s - zi)
    den = np.ones_like(s, dtype=complex)
    for pi in p: den *= (s - pi)
    return num / den


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("csv", help="CSV with columns: freq_hz, mag_db, phase_deg")
    ap.add_argument("--poles",   type=int, default=4)
    ap.add_argument("--zeros",   type=int, default=None,
                    help="scipy backend only (default: poles-1)")
    ap.add_argument("--backend", choices=("auto", "iirrational", "scipy"),
                    default="auto")
    ap.add_argument("--out",     default="plant_fit.npz",
                    help="NPZ output for ZPK (fs=0 denotes continuous-time)")
    ap.add_argument("--mag-only", action="store_true",
                    help="fit magnitude only; ignore measured phase (scipy backend)")
    args = ap.parse_args()

    if not os.path.exists(args.csv):
        sys.exit(f"ERROR: {args.csv} not found")

    f_hz, H_meas = load_bode_csv(args.csv)
    print(f"Loaded {len(f_hz)} points from {args.csv}  "
          f"({f_hz[0]:.1f} Hz – {f_hz[-1]/1e3:.1f} kHz)")

    zeros_order = args.zeros if args.zeros is not None else max(1, args.poles - 1)

    backend = args.backend
    if backend == "auto":
        try:
            import IIRrational  # noqa: F401
            backend = "iirrational"
        except ImportError:
            backend = "scipy"
    print(f"Backend: {backend}")

    if backend == "iirrational":
        z, p, k = fit_iirrational(f_hz, H_meas, order=args.poles)
    else:
        if args.mag_only:
            z, p, k = fit_levi_mag(f_hz, H_meas, args.poles, zeros_order)
        else:
            z, p, k = fit_levi(f_hz, H_meas, args.poles, zeros_order)
    print(f"Fit: {len(p)} poles, {len(z)} zeros, k={k:.4g}")
    print("  Poles (rad/s):", p)
    print("  Zeros (rad/s):", z)

    H_fit = eval_zpk(z, p, k, f_hz)
    err   = np.abs(H_fit - H_meas) / (np.abs(H_meas) + 1e-30)
    print(f"Max relative error: {20*np.log10(np.max(err)+1e-30):.1f} dB")

    np.savez(args.out, z=z, p=p, k=k, fs=0.0)
    print(f"Saved ZPK to {args.out}")

    # --- plot -------------------------------------------------------------
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    fig.suptitle(f"Bode fit: {os.path.basename(args.csv)}  "
                 f"({len(p)}p / {len(z)}z, {backend})")

    a1.semilogx(f_hz, 20*np.log10(np.abs(H_meas)+1e-30),
                "o", ms=3, alpha=0.5, label="measured")
    a1.semilogx(f_hz, 20*np.log10(np.abs(H_fit) +1e-30),
                lw=2, label="fit")
    a1.set_ylabel("Magnitude (dB)"); a1.grid(True, which="both", alpha=0.3); a1.legend()

    a2.semilogx(f_hz, np.degrees(np.unwrap(np.angle(H_meas))),
                "o", ms=3, alpha=0.5, label="measured")
    a2.semilogx(f_hz, np.degrees(np.unwrap(np.angle(H_fit))),
                lw=2, label="fit")
    a2.set_xlabel("Frequency (Hz)"); a2.set_ylabel("Phase (deg)")
    a2.grid(True, which="both", alpha=0.3); a2.legend()

    plt.tight_layout()
    plt.savefig("bode_fit.png", dpi=120)
    print("Saved bode_fit.png")

    plt.show()


if __name__ == "__main__":
    main()
