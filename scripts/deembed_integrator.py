#!/usr/bin/env python3
"""
Remove a known pure-integrator controller from measured open-loop Bode data.

This is intended for the workflow:

    measured loop L(jw) ~= C_int(jw) * P(jw)
    C_int(s) = sign * K / s

If K is not known, the script estimates the composite low-frequency constant
K*P0 from the median of |L(jw)|*w over an integrator-dominated band. The output
is a normalized plant CSV for scripts/bode_fit.py:

    freq_hz, mag_db, phase_deg

Because K and the plant DC gain P0 are not separately identifiable from open
loop data alone, the default behavior normalizes the de-embedded plant to 0 dB
over the estimate band. That is usually the right input to H-infinity synthesis;
the controller produced by ctrl.py then carries the compensator gain.

Examples:
    python3 scripts/deembed_integrator.py openloop.csv plant_for_fit.csv \\
        --estimate-min 2000 --estimate-max 5000 --no-show

    python3 scripts/bode_fit.py plant_for_fit.csv --poles 4 --zeros 3 \\
        --backend iirrational --out plant_fit.npz
"""

from __future__ import annotations

import argparse
import csv
import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib.pyplot as plt
import numpy as np


def normalized_name(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", name.lower())


def pick_column(names: tuple[str, ...], aliases: set[str], contains: tuple[str, ...]) -> str:
    normalized = {normalized_name(name): name for name in names}
    for alias in aliases:
        if alias in normalized:
            return normalized[alias]
    for key, original in normalized.items():
        if all(piece in key for piece in contains):
            return original
    raise KeyError(f"could not find a column matching {sorted(aliases)}")


def read_bode_csv(path: Path, freq_unit: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    data = np.genfromtxt(path, delimiter=",", names=True, dtype=float, encoding=None)
    if data.dtype.names is None:
        raise ValueError(f"{path}: missing CSV header")

    names = data.dtype.names
    freq_col = pick_column(
        names,
        {"freqhz", "frequencyhz", "frequencyinhz", "frequency"},
        contains=("freq",),
    )
    mag_col = pick_column(
        names,
        {"magdb", "gaindb", "magnitudedb", "gain"},
        contains=("db",),
    )
    phase_col = pick_column(
        names,
        {"phasedeg", "phase", "phasedegrees"},
        contains=("phase",),
    )

    freq_in = np.asarray(data[freq_col], dtype=float)
    mag_db = np.asarray(data[mag_col], dtype=float)
    phase_deg = np.asarray(data[phase_col], dtype=float)

    freq_key = normalized_name(freq_col)
    if freq_unit == "log10hz" or (freq_unit == "auto" and "log" in freq_key):
        freq_hz = np.power(10.0, freq_in)
    elif freq_unit == "auto" and np.nanmax(freq_in) < 20.0 and np.nanmin(freq_in) > 0.0:
        raise ValueError(
            f"{path}: frequency values look logarithmic, but column {freq_col!r} "
            "does not say Log. Re-run with --freq-unit log10hz if that is right."
        )
    else:
        freq_hz = freq_in

    keep = np.isfinite(freq_hz) & np.isfinite(mag_db) & np.isfinite(phase_deg) & (freq_hz > 0)
    freq_hz, mag_db, phase_deg = freq_hz[keep], mag_db[keep], phase_deg[keep]
    order = np.argsort(freq_hz)
    return freq_hz[order], mag_db[order], phase_deg[order]


def complex_from_bode(mag_db: np.ndarray, phase_deg: np.ndarray) -> np.ndarray:
    return np.power(10.0, mag_db / 20.0) * np.exp(1j * np.deg2rad(phase_deg))


def principal_phase_deg(phase_deg: np.ndarray) -> np.ndarray:
    return (phase_deg + 180.0) % 360.0 - 180.0


def phase_for_output(H: np.ndarray, mode: str) -> np.ndarray:
    phase = np.rad2deg(np.angle(H))
    if mode == "principal":
        return principal_phase_deg(phase)
    if mode == "unwrap":
        return np.rad2deg(np.unwrap(np.angle(H)))
    raise ValueError(f"unknown phase output mode: {mode}")


def write_bode_csv(path: Path, freq_hz: np.ndarray, H: np.ndarray, phase_mode: str) -> None:
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["freq_hz", "mag_db", "phase_deg"])
        mag_db = 20.0 * np.log10(np.abs(H) + 1e-300)
        phase_deg = phase_for_output(H, phase_mode)
        for row in zip(freq_hz, mag_db, phase_deg):
            writer.writerow([f"{row[0]:.12g}", f"{row[1]:.12g}", f"{row[2]:.12g}"])


def masked_range(freq_hz: np.ndarray, f_min: float | None, f_max: float | None) -> np.ndarray:
    mask = np.ones_like(freq_hz, dtype=bool)
    if f_min is not None:
        mask &= freq_hz >= f_min
    if f_max is not None:
        mask &= freq_hz <= f_max
    return mask


def estimate_integrator_gain_db(
    freq_hz: np.ndarray,
    mag_db: np.ndarray,
    mask: np.ndarray,
    target_plant_dc_db: float,
) -> tuple[float, float]:
    if np.count_nonzero(mask) < 3:
        raise ValueError("estimate band has fewer than three points")

    w = 2.0 * np.pi * freq_hz[mask]
    flat_db = mag_db[mask] + 20.0 * np.log10(w)
    gain_db = float(np.median(flat_db) - target_plant_dc_db)

    x = np.log10(freq_hz[mask])
    slope_db_per_dec, _ = np.polyfit(x, mag_db[mask], 1)
    return gain_db, float(slope_db_per_dec)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="De-embed C(s)=sign*K/s from open-loop Bode data.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("input_csv", type=Path, help="measured open-loop CSV")
    parser.add_argument("output_csv", type=Path, help="plant CSV for scripts/bode_fit.py")
    parser.add_argument(
        "--freq-unit",
        choices=("auto", "hz", "log10hz"),
        default="auto",
        help="input frequency unit",
    )
    parser.add_argument("--freq-min", type=float, default=None, help="minimum output Hz")
    parser.add_argument("--freq-max", type=float, default=None, help="maximum output Hz")
    parser.add_argument("--estimate-min", type=float, default=None, help="minimum Hz for estimating K")
    parser.add_argument("--estimate-max", type=float, default=None, help="maximum Hz for estimating K")
    parser.add_argument(
        "--integrator-gain",
        type=float,
        default=None,
        help="known K in C(s)=sign*K/s. If omitted, estimate composite K*P0.",
    )
    parser.add_argument(
        "--controller-sign",
        choices=("+1", "-1"),
        default="+1",
        help="sign of the known integrator controller",
    )
    parser.add_argument(
        "--target-plant-dc-db",
        type=float,
        default=0.0,
        help="plant level imposed in estimate band when K is estimated",
    )
    parser.add_argument(
        "--phase-output",
        choices=("principal", "unwrap"),
        default="unwrap",
        help="phase branch written to output CSV",
    )
    parser.add_argument("--plot", default="deembedded_integrator.png", help="diagnostic plot path")
    parser.add_argument("--no-show", action="store_true", help="save plot without opening a GUI window")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    freq_hz, loop_mag_db, loop_phase_deg = read_bode_csv(args.input_csv, args.freq_unit)

    output_mask = masked_range(freq_hz, args.freq_min, args.freq_max)
    estimate_mask = masked_range(freq_hz, args.estimate_min, args.estimate_max) & output_mask
    if not np.any(output_mask):
        raise ValueError("no points remain after output frequency filtering")

    if args.integrator_gain is None:
        gain_db, slope_db_per_dec = estimate_integrator_gain_db(
            freq_hz,
            loop_mag_db,
            estimate_mask,
            args.target_plant_dc_db,
        )
        integrator_gain = 10.0 ** (gain_db / 20.0)
        estimate_note = "estimated"
    else:
        if args.integrator_gain <= 0:
            raise ValueError("--integrator-gain must be positive")
        integrator_gain = args.integrator_gain
        gain_db = 20.0 * np.log10(integrator_gain)
        slope_db_per_dec = float("nan")
        estimate_note = "provided"

    sign = 1.0 if args.controller_sign == "+1" else -1.0
    H_loop = complex_from_bode(loop_mag_db, loop_phase_deg)
    s = 1j * 2.0 * np.pi * freq_hz
    H_plant = H_loop * s / (sign * integrator_gain)

    freq_out = freq_hz[output_mask]
    loop_out = H_loop[output_mask]
    plant_out = H_plant[output_mask]
    write_bode_csv(args.output_csv, freq_out, plant_out, args.phase_output)

    print(f"Loaded {len(freq_hz)} points from {args.input_csv}")
    print(f"Wrote {len(freq_out)} plant points to {args.output_csv}")
    print(f"Integrator gain ({estimate_note}): K = {integrator_gain:.6g} rad/s  ({gain_db:.2f} dB re rad/s)")
    if np.isfinite(slope_db_per_dec):
        print(f"Estimate-band loop slope: {slope_db_per_dec:.2f} dB/dec  (ideal integrator: -20 dB/dec)")
        if abs(slope_db_per_dec + 20.0) > 5.0:
            print("WARNING: estimate band is not very integrator-like; choose a cleaner --estimate-min/--estimate-max.")

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True)
    fig.suptitle(f"Integrator de-embedding: {args.input_csv.name}")

    axes[0, 0].semilogx(freq_out, 20.0 * np.log10(np.abs(loop_out) + 1e-300), ".", ms=2.5)
    axes[0, 0].set_title("Measured loop L")
    axes[0, 0].set_ylabel("Magnitude (dB)")
    axes[0, 0].grid(True, which="both", alpha=0.35)

    axes[1, 0].semilogx(freq_out, phase_for_output(loop_out, "principal"), ".", ms=2.5)
    axes[1, 0].set_ylabel("Phase (deg)")
    axes[1, 0].set_xlabel("Frequency (Hz)")
    axes[1, 0].set_ylim(-190, 190)
    axes[1, 0].grid(True, which="both", alpha=0.35)

    axes[0, 1].semilogx(freq_out, 20.0 * np.log10(np.abs(plant_out) + 1e-300), ".", ms=2.5, color="C1")
    axes[0, 1].axhline(args.target_plant_dc_db, color="0.35", lw=0.8, ls=":")
    axes[0, 1].set_title("After dividing by sign*K/s")
    axes[0, 1].set_ylabel("Magnitude (dB)")
    axes[0, 1].grid(True, which="both", alpha=0.35)

    axes[1, 1].semilogx(freq_out, phase_for_output(plant_out, "principal"), ".", ms=2.5, color="C1")
    axes[1, 1].set_ylabel("Phase (deg)")
    axes[1, 1].set_xlabel("Frequency (Hz)")
    axes[1, 1].set_ylim(-190, 190)
    axes[1, 1].grid(True, which="both", alpha=0.35)

    if np.any(estimate_mask & output_mask):
        f_est = freq_hz[estimate_mask & output_mask]
        for ax in axes.flat:
            ax.axvspan(f_est.min(), f_est.max(), color="C2", alpha=0.08)

    fig.tight_layout()
    fig.savefig(args.plot, dpi=150)
    print(f"Saved {args.plot}")
    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
