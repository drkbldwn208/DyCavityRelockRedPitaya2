"""
Fit a measured open-loop Bode plot with a rational IIR model.

Usage:
    python3 scripts/bode_fit.py <data.csv> [--poles N] [--zeros M] [--out plant.npz]

    data.csv must have columns:  freq_hz, mag_db, phase_deg

Preferred backend is IIRrational (pip install iirrational) if installed.
The inline SciPy/Levi backend is a fallback sanity check, not a substitute for
IIRrational when sharp resonances or near-canceling pole-zero pairs matter.

Outputs:
    bode_fit.png  — measurement vs fitted model
    <out>.npz     — ZPK: z, p, k, fs (fs=0 marks continuous-time)
"""

import argparse
import collections
import collections.abc
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig

if not hasattr(collections, "Mapping"):
    collections.Mapping = collections.abc.Mapping
if not hasattr(collections, "MutableMapping"):
    collections.MutableMapping = collections.abc.MutableMapping
if not hasattr(collections, "Sequence"):
    collections.Sequence = collections.abc.Sequence


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
    keep = np.isfinite(f_hz) & np.isfinite(mag_db) & np.isfinite(phase_deg) & (f_hz > 0)
    f_hz, mag_db, phase_deg = f_hz[keep], mag_db[keep], phase_deg[keep]
    order = np.argsort(f_hz)
    f_hz, mag_db, phase_deg = f_hz[order], mag_db[order], phase_deg[order]
    mag   = 10 ** (mag_db / 20.0)
    phase = np.radians(phase_deg)
    return f_hz, mag * np.exp(1j * phase), phase_deg

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

def import_iirrational_v2():
    try:
        import IIRrational.v2 as iirr
        return iirr, "IIRrational.v2"
    except ImportError as legacy_err:
        try:
            import wield.iirrational.v2 as iirr
            return iirr, "wield.iirrational.v2"
        except ImportError as wield_err:
            raise ImportError(
                "Could not import either IIRrational.v2 or wield.iirrational.v2. "
                f"Legacy import failed with: {legacy_err}. "
                f"Wield import failed with: {wield_err}."
            ) from wield_err


def resonance_pair_rad_s(freq_hz, q):
    """Return a stable continuous-time complex pair from resonance frequency/Q."""
    if freq_hz <= 0:
        raise ValueError("resonance frequency must be positive")
    if q <= 0:
        raise ValueError("resonance Q must be positive")
    w0 = 2 * np.pi * freq_hz
    sigma = -w0 / (2 * q)
    wd_sq = max(w0 * w0 - sigma * sigma, 0.0)
    wd = np.sqrt(wd_sq)
    return np.array([sigma + 1j * wd, sigma - 1j * wd], dtype=complex)


def pairs_from_args(pair_args):
    roots = []
    for freq_hz, q in pair_args or []:
        roots.extend(resonance_pair_rad_s(freq_hz, q))
    return np.asarray(roots, dtype=complex)


def fit_iirrational(
    f_hz,
    H,
    order,
    mode="full",
    choose_direction="<=",
    order_min=None,
    order_max=None,
    snr=None,
    snr_estimate_width=None,
    snr_min=None,
    emphasis=None,
    total_degree_min=None,
    force_poles=None,
    force_zeros=None,
    seed_poles=None,
    seed_zeros=None,
    suggest=False,
    never_unstable_poles=False,
    never_unstable_zeros=False,
):
    iirr, _ = import_iirrational_v2()
    fit_kwargs = {"data": H, "F_Hz": f_hz, "mode": mode}
    if order_min is not None:
        fit_kwargs["order_min"] = order_min
    if order_max is not None:
        fit_kwargs["order_max"] = order_max
    if total_degree_min is not None:
        fit_kwargs["total_degree_min"] = total_degree_min
    if snr is not None:
        fit_kwargs["SNR"] = snr
    if snr_estimate_width is not None:
        fit_kwargs["SNR_estimate_width"] = snr_estimate_width
    if snr_min is not None:
        fit_kwargs["SNR_min"] = snr_min
    if emphasis is not None:
        fit_kwargs["emphasis"] = emphasis
    if force_poles is not None and len(force_poles):
        fit_kwargs["poles_overlay"] = force_poles
    if force_zeros is not None and len(force_zeros):
        fit_kwargs["zeros_overlay"] = force_zeros
    if seed_poles is not None and len(seed_poles):
        fit_kwargs["poles"] = seed_poles
        suggest = True
    if seed_zeros is not None and len(seed_zeros):
        fit_kwargs["zeros"] = seed_zeros
        suggest = True
    if suggest:
        fit_kwargs["suggest"] = True
    if never_unstable_poles:
        fit_kwargs["never_unstable_poles"] = True
    if never_unstable_zeros:
        fit_kwargs["never_unstable_zeros"] = True

    fit = iirr.data2filter(**fit_kwargs)
    try:
        orders, residuals = fit.residuals_by_order()
        print("IIRrational explored orders:")
        for oi, ri in zip(orders, residuals):
            print(f"  order {int(oi):2d}: residual {ri:.4g}")
    except Exception:
        pass

    chosen = fit.choose(order=order, direction=choose_direction)
    if hasattr(fit, "as_scipy_signal_ZPKsw"):
        z, p, k = fit.as_scipy_signal_ZPKsw()
    elif all(hasattr(fit, name) for name in ("zeros", "poles", "gain")):
        z, p, k = fit.zeros, fit.poles, fit.gain
    elif chosen is not None and all(hasattr(chosen, name) for name in ("zeros", "poles", "gain")):
        z, p, k = chosen.zeros, chosen.poles, chosen.gain
    else:
        raise AttributeError("Could not extract ZPK from IIRrational result object")
    z = np.asarray(z, dtype=complex)
    p = np.asarray(p, dtype=complex)
    k = float(k)
    return z, p, k


def build_emphasis(f_hz, args):
    if args.emphasis_min is None and args.emphasis_max is None:
        return None
    emphasis = np.ones_like(f_hz, dtype=float)
    if args.emphasis_min is not None:
        emphasis[f_hz < args.emphasis_min] = 1.0
    if args.emphasis_max is not None:
        emphasis[f_hz > args.emphasis_max] = 1.0
    if args.emphasis_min is not None or args.emphasis_max is not None:
        mask = np.ones_like(f_hz, dtype=bool)
        if args.emphasis_min is not None:
            mask &= f_hz >= args.emphasis_min
        if args.emphasis_max is not None:
            mask &= f_hz <= args.emphasis_max
        emphasis[mask] = args.emphasis_gain
    return emphasis


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
    ap.add_argument("--iirr-mode",
                    choices=("full", "rational", "rational2x", "full2x", "AAA", "fit"),
                    default="full",
                    help="IIRrational fitting mode")
    ap.add_argument("--iirr-choose-direction",
                    choices=("<=", "le", ">=", "ge", "<", "l", ">", "g"),
                    default="<=",
                    help="IIRrational order choice direction; use ge/>= to request at least --poles")
    ap.add_argument("--order-min", type=int, default=None,
                    help="IIRrational initial fit minimum order hint")
    ap.add_argument("--order-max", type=int, default=None,
                    help="IIRrational initial fit maximum order hint")
    ap.add_argument("--total-degree-min", type=int, default=None,
                    help="IIRrational minimum degree during order reduction")
    ap.add_argument("--flat-snr", action="store_true",
                    help="IIRrational: use uniform SNR weights instead of estimating/dropping low-SNR points")
    ap.add_argument("--snr-estimate-width", type=int, default=None,
                    help="IIRrational SNR_estimate_width; use 0 to disable SNR estimation")
    ap.add_argument("--snr-min", type=float, default=None,
                    help="IIRrational SNR_min threshold; use 0 to avoid dropping low-SNR points")
    ap.add_argument("--emphasis-min", type=float, default=None,
                    help="minimum Hz of an emphasized fit band")
    ap.add_argument("--emphasis-max", type=float, default=None,
                    help="maximum Hz of an emphasized fit band")
    ap.add_argument("--emphasis-gain", type=float, default=5.0,
                    help="IIRrational weight multiplier inside the emphasized band")
    ap.add_argument("--force-pole-pair", type=float, nargs=2, action="append",
                    metavar=("F_HZ", "Q"),
                    help="IIRrational: fixed/overlay stable pole pair at resonance frequency F_HZ with Q")
    ap.add_argument("--force-zero-pair", type=float, nargs=2, action="append",
                    metavar=("F_HZ", "Q"),
                    help="IIRrational: fixed/overlay stable zero pair at resonance frequency F_HZ with Q")
    ap.add_argument("--seed-pole-pair", type=float, nargs=2, action="append",
                    metavar=("F_HZ", "Q"),
                    help="IIRrational: initial stable pole-pair suggestion near F_HZ with Q")
    ap.add_argument("--seed-zero-pair", type=float, nargs=2, action="append",
                    metavar=("F_HZ", "Q"),
                    help="IIRrational: initial stable zero-pair suggestion near F_HZ with Q")
    ap.add_argument("--iirr-suggest", action="store_true",
                    help="IIRrational: treat provided initial ZPK roots as suggestions")
    ap.add_argument("--never-unstable-poles", action="store_true",
                    help="IIRrational: prevent phase patching from adding unstable poles")
    ap.add_argument("--never-unstable-zeros", action="store_true",
                    help="IIRrational: prevent phase patching from adding unstable zeros")
    ap.add_argument("--freq-min", type=float, default=None,
                    help="minimum fitted frequency, Hz")
    ap.add_argument("--freq-max", type=float, default=None,
                    help="maximum fitted frequency, Hz")
    ap.add_argument("--out",     default="plant_fit.npz",
                    help="NPZ output for ZPK (fs=0 denotes continuous-time)")
    ap.add_argument("--plot-out", default="bode_fit.png",
                    help="diagnostic plot path")
    ap.add_argument("--mag-only", action="store_true",
                    help="fit magnitude only; ignore measured phase (scipy backend)")
    ap.add_argument("--no-show", action="store_true",
                    help="save plot without opening a GUI window")
    args = ap.parse_args()

    if not os.path.exists(args.csv):
        sys.exit(f"ERROR: {args.csv} not found")

    f_hz, H_meas, phase_meas_deg = load_bode_csv(args.csv)
    fit_mask = np.ones_like(f_hz, dtype=bool)
    if args.freq_min is not None:
        fit_mask &= f_hz >= args.freq_min
    if args.freq_max is not None:
        fit_mask &= f_hz <= args.freq_max
    f_hz, H_meas, phase_meas_deg = f_hz[fit_mask], H_meas[fit_mask], phase_meas_deg[fit_mask]
    if len(f_hz) < max(args.poles + 1, 4):
        sys.exit("ERROR: not enough points remain after frequency filtering")
    print(f"Loaded {len(f_hz)} points from {args.csv}  "
          f"({f_hz[0]:.1f} Hz – {f_hz[-1]/1e3:.1f} kHz)")

    zeros_order = args.zeros if args.zeros is not None else max(1, args.poles - 1)

    backend = args.backend
    if backend == "auto":
        try:
            _, iirr_module = import_iirrational_v2()
            backend = "iirrational"
            print(f"IIRrational module: {iirr_module}")
        except ImportError:
            backend = "scipy"
            print("WARNING: IIRrational is not installed; falling back to simple SciPy/Levi fit.")
            print("         Expect poor fits for sharp resonances or near pole-zero cancellations.")
    print(f"Backend: {backend}")

    if backend == "iirrational":
        try:
            snr = np.ones_like(f_hz) if args.flat_snr else None
            snr_estimate_width = 0 if args.flat_snr else args.snr_estimate_width
            snr_min = 0 if args.flat_snr else args.snr_min
            emphasis = build_emphasis(f_hz, args)
            force_poles = pairs_from_args(args.force_pole_pair)
            force_zeros = pairs_from_args(args.force_zero_pair)
            seed_poles = pairs_from_args(args.seed_pole_pair)
            seed_zeros = pairs_from_args(args.seed_zero_pair)
            z, p, k = fit_iirrational(
                f_hz,
                H_meas,
                order=args.poles,
                mode=args.iirr_mode,
                choose_direction=args.iirr_choose_direction,
                order_min=args.order_min,
                order_max=args.order_max,
                snr=snr,
                snr_estimate_width=snr_estimate_width,
                snr_min=snr_min,
                emphasis=emphasis,
                total_degree_min=args.total_degree_min,
                force_poles=force_poles,
                force_zeros=force_zeros,
                seed_poles=seed_poles,
                seed_zeros=seed_zeros,
                suggest=args.iirr_suggest,
                never_unstable_poles=args.never_unstable_poles,
                never_unstable_zeros=args.never_unstable_zeros,
            )
        except ImportError as exc:
            sys.exit(f"ERROR: IIRrational backend requested but not installed: {exc}")
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

    phase_fit_deg = np.degrees(np.unwrap(np.angle(H_fit)))
    phase_fit_deg += 360.0 * np.round((phase_meas_deg[0] - phase_fit_deg[0]) / 360.0)
    a2.semilogx(f_hz, phase_meas_deg,
                "o", ms=3, alpha=0.5, label="measured")
    a2.semilogx(f_hz, phase_fit_deg,
                lw=2, label="fit")
    a2.set_xlabel("Frequency (Hz)"); a2.set_ylabel("Phase (deg)")
    a2.grid(True, which="both", alpha=0.3); a2.legend()

    plt.tight_layout()
    plt.savefig(args.plot_out, dpi=120)
    print(f"Saved {args.plot_out}")

    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
