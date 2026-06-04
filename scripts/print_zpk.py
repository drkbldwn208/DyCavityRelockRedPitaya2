"""
Print poles, zeros, and gain from a continuous-time or discrete-time ZPK NPZ.

Expected NPZ keys:
    z, p, k

Usage:
    python3 scripts/print_zpk.py plant_fit.npz
    python3 scripts/print_zpk.py K_zpk.npz --discrete
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


TWO_PI = 2.0 * np.pi


def cfmt(x: complex, precision: int = 8) -> str:
    sign = "+" if x.imag >= 0 else "-"
    return f"{x.real:.{precision}g} {sign} {abs(x.imag):.{precision}g}j"


def root_rows(roots: np.ndarray, discrete: bool) -> list[tuple[float, complex, str]]:
    rows = []
    for root in np.asarray(roots, dtype=complex).flatten():
        if discrete:
            key = abs(root)
            detail = f"|z|={abs(root):.8g}, angle={np.angle(root):+.8g} rad"
        else:
            key = abs(root)
            real_hz = root.real / TWO_PI
            imag_hz = root.imag / TWO_PI
            fn_hz = abs(root) / TWO_PI
            if abs(root.imag) > 1e-12 and root.real < 0:
                q = abs(root) / (-2.0 * root.real)
                detail = f"Re={real_hz:+.8g} Hz, Im={imag_hz:+.8g} Hz, fn={fn_hz:.8g} Hz, Q={q:.8g}"
            else:
                detail = f"Re={real_hz:+.8g} Hz, Im={imag_hz:+.8g} Hz, |f|={fn_hz:.8g} Hz"
        rows.append((key, root, detail))
    return sorted(rows, key=lambda row: row[0])


def dc_gain(z: np.ndarray, p: np.ndarray, k: complex, discrete: bool) -> complex | None:
    try:
        if discrete:
            den = np.prod(1.0 - p) if len(p) else 1.0
            if abs(den) == 0:
                return None
            return k * (np.prod(1.0 - z) if len(z) else 1.0) / den
        den = np.prod(-p) if len(p) else 1.0
        if abs(den) == 0:
            return None
        return k * (np.prod(-z) if len(z) else 1.0) / den
    except FloatingPointError:
        return None


def print_roots(name: str, roots: np.ndarray, discrete: bool) -> None:
    print(f"\n{name}: {len(roots)}")
    for i, (_, root, detail) in enumerate(root_rows(roots, discrete)):
        print(f"  {i:2d}: {cfmt(root):>28s}    {detail}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("npz", type=Path)
    parser.add_argument("--discrete", action="store_true", help="interpret roots as z-plane roots")
    args = parser.parse_args()

    data = np.load(args.npz, allow_pickle=False)
    z = np.asarray(data["z"], dtype=complex).flatten()
    p = np.asarray(data["p"], dtype=complex).flatten()
    k = complex(np.asarray(data["k"]).item())

    if "fs" in data.files and not args.discrete:
        roots = np.concatenate([z, p]) if len(z) or len(p) else np.asarray([], dtype=complex)
        if len(roots) and np.nanmax(np.abs(roots)) < 2.0:
            print("NOTE: roots look like discrete z-plane values; consider --discrete.")

    print(f"File: {args.npz}")
    print(f"Keys: {', '.join(data.files)}")
    print(f"Interpretation: {'discrete z-plane' if args.discrete else 'continuous rad/s'}")
    print(f"Gain k: {cfmt(k)}")

    dc = dc_gain(z, p, k, args.discrete)
    if dc is None:
        print("DC gain: undefined/infinite")
    else:
        mag_db = 20.0 * np.log10(abs(dc) + 1e-300)
        print(f"DC gain: {cfmt(dc)}    |DC|={abs(dc):.8g}, {mag_db:+.3f} dB")

    print_roots("Zeros", z, args.discrete)
    print_roots("Poles", p, args.discrete)


if __name__ == "__main__":
    main()
