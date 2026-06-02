#!/usr/bin/env python3
"""
Plot measured Bode data and optionally export it for scripts/bode_fit.py.

The instrument CSVs seen so far use:
    Frequency Log(Hz), Gain (dB), Phase (deg)

The fitter uses:
    freq_hz, mag_db, phase_deg

This script accepts either format. By default it plots unsmoothed data with a
raw connecting line and principal phase in the [-180, 180) degree interval,
which keeps the display scale reasonable when the recorded phase branch count
has wandered. Median trend overlays are opt-in and never affect exported data. Use
--phase-mode unwrap to show a continuous trace; this converts the input degrees
onto the unit circle and uses numpy.unwrap, repairing +/-180 deg branch cuts
while preserving an already-continuous phase trace.

Examples:
    python3 scripts/plot_bode_data.py ~/Downloads/CavityOpenLoopXfer2k15k999Steps.csv

    python3 scripts/plot_bode_data.py ~/Downloads/CavityOpenLoopXfer*.csv \\
        --out cavity_open_loop_bode.png

    python3 scripts/plot_bode_data.py ~/Downloads/CavityOpenLoopXfer2k15k999Steps.csv \\
        --fit-csv openloop_for_fit.csv --freq-min 2000 --freq-max 15000
"""

from __future__ import annotations

import argparse
import csv
import os
import re
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib.pyplot as plt
import numpy as np


@dataclass
class BodeTrace:
    path: Path
    freq_hz: np.ndarray
    mag_db: np.ndarray
    phase_raw_deg: np.ndarray
    phase_plot_deg: np.ndarray
    phase_fit_deg: np.ndarray


def normalized_name(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", name.lower())


def pick_column(names: tuple[str, ...], aliases: set[str], contains: tuple[str, ...] = ()) -> str:
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


def unwrap_phase_deg(phase_deg: np.ndarray, jump_deg: float) -> np.ndarray:
    return np.rad2deg(np.unwrap(np.deg2rad(phase_deg), discont=np.deg2rad(jump_deg)))


def principal_phase_deg(phase_deg: np.ndarray) -> np.ndarray:
    return (phase_deg + 180.0) % 360.0 - 180.0


def continuous_principal_phase_deg(phase_deg: np.ndarray, jump_deg: float) -> np.ndarray:
    return unwrap_phase_deg(principal_phase_deg(phase_deg), jump_deg)


def phase_reference_index(freq_hz: np.ndarray, spec: str) -> int:
    if spec == "first":
        return 0
    if spec == "last":
        return len(freq_hz) - 1
    target = float(spec)
    return int(np.argmin(np.abs(freq_hz - target)))


def shift_phase_branch(
    phase_deg: np.ndarray,
    freq_hz: np.ndarray,
    target_deg: float | None,
    reference: str,
) -> np.ndarray:
    if target_deg is None:
        return phase_deg
    idx = phase_reference_index(freq_hz, reference)
    offset = 360.0 * np.round((target_deg - phase_deg[idx]) / 360.0)
    return phase_deg + offset


def rolling_median(y: np.ndarray, window: int) -> np.ndarray:
    if window <= 1:
        return y
    window = max(3, int(window))
    if window % 2 == 0:
        window += 1
    half = window // 2
    out = np.empty_like(y, dtype=float)
    for i in range(len(y)):
        lo, hi = max(0, i - half), min(len(y), i + half + 1)
        out[i] = np.nanmedian(y[lo:hi])
    return out


def rolling_phase_median_deg(y: np.ndarray, window: int, principal_output: bool) -> np.ndarray:
    if window <= 1:
        return y
    window = max(3, int(window))
    if window % 2 == 0:
        window += 1
    half = window // 2
    out = np.empty_like(y, dtype=float)
    for i in range(len(y)):
        lo, hi = max(0, i - half), min(len(y), i + half + 1)
        vals = y[lo:hi]
        vals_near_center = vals + 360.0 * np.round((y[i] - vals) / 360.0)
        out[i] = np.nanmedian(vals_near_center)
    if principal_output:
        out = principal_phase_deg(out)
    return out


def break_phase_branch_cuts(y: np.ndarray, jump_deg: float = 60.0) -> np.ndarray:
    out = y.copy()
    jumps = np.abs(np.diff(out)) > jump_deg
    out[1:][jumps] = np.nan
    return out


def load_trace(path: Path, args: argparse.Namespace) -> BodeTrace:
    freq_hz, mag_db, phase_raw_deg = read_bode_csv(path, args.freq_unit)

    mask = np.ones_like(freq_hz, dtype=bool)
    if args.freq_min is not None:
        mask &= freq_hz >= args.freq_min
    if args.freq_max is not None:
        mask &= freq_hz <= args.freq_max
    freq_hz, mag_db, phase_raw_deg = freq_hz[mask], mag_db[mask], phase_raw_deg[mask]
    if len(freq_hz) == 0:
        raise ValueError(f"{path}: no points remain after frequency filtering")

    if args.phase_mode == "raw" or args.no_unwrap:
        phase_plot_deg = phase_raw_deg.copy()
    elif args.phase_mode == "principal":
        phase_plot_deg = principal_phase_deg(phase_raw_deg)
    elif args.phase_mode == "unwrap":
        phase_plot_deg = unwrap_phase_deg(phase_raw_deg, args.unwrap_jump)
    else:
        raise ValueError(f"unhandled phase mode: {args.phase_mode}")
    phase_plot_deg = shift_phase_branch(
        phase_plot_deg,
        freq_hz,
        args.phase_anchor,
        args.phase_anchor_at,
    )
    if args.fit_phase_mode == "raw":
        phase_fit_deg = phase_raw_deg.copy()
    elif args.fit_phase_mode == "principal":
        phase_fit_deg = principal_phase_deg(phase_raw_deg)
    elif args.fit_phase_mode == "unwrap":
        phase_fit_deg = continuous_principal_phase_deg(phase_raw_deg, args.unwrap_jump)
    else:
        raise ValueError(f"unhandled fit phase mode: {args.fit_phase_mode}")
    phase_fit_deg = shift_phase_branch(
        phase_fit_deg,
        freq_hz,
        args.phase_anchor,
        args.phase_anchor_at,
    )
    return BodeTrace(path, freq_hz, mag_db, phase_raw_deg, phase_plot_deg, phase_fit_deg)


def default_label(path: Path) -> str:
    stem = path.stem
    if len(stem) > 42:
        return stem[:18] + "..." + stem[-18:]
    return stem


def write_fit_csv(trace: BodeTrace, out_path: Path) -> None:
    with out_path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["freq_hz", "mag_db", "phase_deg"])
        for row in zip(trace.freq_hz, trace.mag_db, trace.phase_fit_deg):
            writer.writerow([f"{row[0]:.12g}", f"{row[1]:.12g}", f"{row[2]:.12g}"])


def output_fit_paths(inputs: list[Path], requested: str) -> list[Path]:
    out = Path(requested)
    if len(inputs) == 1:
        return [out]
    if out.suffix.lower() == ".csv":
        stem = out.with_suffix("")
        return [stem.parent / f"{stem.name}_{path.stem}.csv" for path in inputs]
    out.mkdir(parents=True, exist_ok=True)
    return [out / f"{path.stem}_for_bode_fit.csv" for path in inputs]


def plot_traces(traces: list[BodeTrace], args: argparse.Namespace) -> None:
    plt.style.use("seaborn-v0_8-whitegrid")
    fig, (ax_mag, ax_phase) = plt.subplots(2, 1, figsize=(11, 7.5), sharex=True)
    fig.suptitle(args.title or "Open-loop Bode measurement", fontsize=14)

    for i, trace in enumerate(traces):
        color = f"C{i % 10}"
        label = default_label(trace.path)
        show_raw_in_legend = args.median_window <= 1 and args.line_width <= 0
        raw_label = label if show_raw_in_legend else "_nolegend_"
        ax_mag.semilogx(
            trace.freq_hz,
            trace.mag_db,
            ".",
            ms=args.marker_size,
            alpha=args.raw_alpha,
            color=color,
            label=raw_label,
        )
        ax_phase.semilogx(
            trace.freq_hz,
            trace.phase_plot_deg,
            ".",
            ms=args.marker_size,
            alpha=args.raw_alpha,
            color=color,
            label=raw_label,
        )

        if args.line_width > 0:
            phase_line = trace.phase_plot_deg
            if args.phase_mode == "principal":
                phase_line = break_phase_branch_cuts(phase_line, args.phase_break)
            ax_mag.semilogx(
                trace.freq_hz,
                trace.mag_db,
                lw=args.line_width,
                alpha=args.line_alpha,
                color=color,
                label=label,
            )
            ax_phase.semilogx(
                trace.freq_hz,
                phase_line,
                lw=args.line_width,
                alpha=args.line_alpha,
                color=color,
                label=label,
            )

        if args.median_window > 1:
            mag_trend = rolling_median(trace.mag_db, args.median_window)
            phase_trend = rolling_phase_median_deg(
                trace.phase_plot_deg,
                args.median_window,
                principal_output=args.phase_mode == "principal",
            )
            ax_mag.semilogx(trace.freq_hz, mag_trend, lw=2.0, color=color, label=label)
            if args.phase_mode == "principal":
                phase_trend = break_phase_branch_cuts(phase_trend)
            ax_phase.semilogx(trace.freq_hz, phase_trend, lw=2.0, color=color, label=label)

    ax_mag.axhline(0.0, color="0.25", lw=0.9, ls="--", alpha=0.7)
    ax_mag.set_ylabel("Magnitude (dB)")
    ax_mag.grid(True, which="both", alpha=0.35)
    ax_mag.legend(fontsize=8, loc="best")

    ax_phase.axhline(-180.0, color="0.35", lw=0.8, ls=":", alpha=0.7)
    ax_phase.axhline(-90.0, color="0.35", lw=0.8, ls=":", alpha=0.6)
    ax_phase.set_xlabel("Frequency (Hz)")
    ax_phase.set_ylabel("Phase (deg)")
    ax_phase.grid(True, which="both", alpha=0.35)
    if args.phase_mode == "principal" and args.phase_ylim and args.phase_ylim > 0:
        ax_phase.set_ylim(-args.phase_ylim, args.phase_ylim)
    ax_phase.legend(fontsize=8, loc="best")

    if args.phase_anchor is not None:
        ax_phase.text(
            0.01,
            0.02,
            f"phase branch shifted near {args.phase_anchor:g} deg at {args.phase_anchor_at}",
            transform=ax_phase.transAxes,
            fontsize=8,
            color="0.35",
        )

    fig.tight_layout()
    out = Path(args.out)
    fig.savefig(out, dpi=args.dpi)
    print(f"Saved {out}")
    if not args.no_show:
        plt.show()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot measured transfer-function Bode CSVs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("csv", nargs="+", type=Path, help="input Bode CSV file(s)")
    parser.add_argument("--out", default="bode_data.png", help="output plot path")
    parser.add_argument("--title", default=None, help="figure title")
    parser.add_argument(
        "--freq-unit",
        choices=("auto", "hz", "log10hz"),
        default="auto",
        help="input frequency unit",
    )
    parser.add_argument("--freq-min", type=float, default=None, help="minimum plotted/exported Hz")
    parser.add_argument("--freq-max", type=float, default=None, help="maximum plotted/exported Hz")
    parser.add_argument(
        "--phase-mode",
        choices=("principal", "unwrap", "raw"),
        default="principal",
        help="phase display branch handling",
    )
    parser.add_argument(
        "--fit-phase-mode",
        choices=("unwrap", "principal", "raw"),
        default="unwrap",
        help="phase branch written by --fit-csv; smoothing is never applied",
    )
    parser.add_argument(
        "--no-unwrap",
        action="store_true",
        help="deprecated alias for --phase-mode raw",
    )
    parser.add_argument(
        "--unwrap-jump",
        type=float,
        default=180.0,
        help="phase jump, in degrees, treated as a branch cut",
    )
    parser.add_argument(
        "--phase-anchor",
        type=float,
        default=None,
        help="shift phase by an integer multiple of 360 deg near this target",
    )
    parser.add_argument(
        "--phase-anchor-at",
        default="first",
        help="'first', 'last', or a frequency in Hz used by --phase-anchor",
    )
    parser.add_argument(
        "--median-window",
        type=int,
        default=1,
        help="optional odd-point rolling median trend window; 1 means no smoothing",
    )
    parser.add_argument(
        "--phase-ylim",
        type=float,
        default=190.0,
        help="symmetric phase y-limit for --phase-mode principal; use 0 for autoscale",
    )
    parser.add_argument(
        "--phase-break",
        type=float,
        default=120.0,
        help="do not connect principal-phase line segments across jumps larger than this many degrees",
    )
    parser.add_argument("--line-width", type=float, default=0.8, help="raw unsmoothed line width; use 0 for points only")
    parser.add_argument("--line-alpha", type=float, default=0.85, help="raw unsmoothed line alpha")
    parser.add_argument("--marker-size", type=float, default=2.8, help="raw data marker size")
    parser.add_argument("--raw-alpha", type=float, default=0.28, help="raw data marker alpha")
    parser.add_argument("--dpi", type=int, default=160, help="saved figure DPI")
    parser.add_argument(
        "--fit-csv",
        default=None,
        help="write unsmoothed bode_fit.py-ready CSV(s). With multiple inputs, pass a directory or a base CSV name.",
    )
    parser.add_argument("--no-show", action="store_true", help="save without opening a GUI window")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    traces = [load_trace(path, args) for path in args.csv]

    print("Loaded:")
    for trace in traces:
        print(
            f"  {trace.path}: {len(trace.freq_hz)} points, "
            f"{trace.freq_hz[0]:.3g} Hz to {trace.freq_hz[-1]:.3g} Hz"
        )

    if args.fit_csv:
        for trace, out_path in zip(traces, output_fit_paths([t.path for t in traces], args.fit_csv)):
            write_fit_csv(trace, out_path)
            print(f"Wrote {out_path}")

    plot_traces(traces, args)


if __name__ == "__main__":
    main()
