"""
Closed-loop simulation using a fitted plant and the actual fixed-point controller.

This models the controller at the decimated HLS controller rate, not at the
125 MHz ADC/DAC fabric rate.  That keeps long runs practical while preserving
the important fixed-point controller behavior:

    ADC counts -> HinfFilter fixed-point SOS -> DAC counts -> fitted plant -> ADC counts

The plant NPZ is interpreted as continuous-time rad/s ZPK, then discretized at
the controller update rate.  The controller is read from src/hinf_coeffs.h and
uses the same AP_TRN_ZERO integer arithmetic as scripts/controller_fixedpoint_sim.py.

Usage:
    python3 scripts/closed_loop_fixedpoint_sim.py plant_fit.npz
    python3 scripts/closed_loop_fixedpoint_sim.py plant_fit.npz --disturbance-counts 100 --samples 200000
    python3 scripts/closed_loop_fixedpoint_sim.py plant_fit.npz --controller-sign -1
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as sig

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(SCRIPT_DIR))

from controller_fixedpoint_sim import (  # noqa: E402
    ACC_FRAC,
    COEF_FRAC,
    PIPE_FRAC,
    PIPE_SHIFT,
    SIG_MAX,
    SIG_MIN,
    parse_hinf_coeffs,
    qshift_zero,
    sat_signed,
)


DAC_MIN = -8192
DAC_MAX = 8191
TWO_PI = 2.0 * math.pi


def sat_dac14(x: int) -> int:
    return max(DAC_MIN, min(DAC_MAX, int(x)))


def wrap_signed(x: int, bits: int) -> int:
    mask = (1 << bits) - 1
    x = int(x) & mask
    sign = 1 << (bits - 1)
    return x - (1 << bits) if x & sign else x


def cast_adc16(x: float, mode: str) -> int:
    xi = int(round(float(x)))
    if mode == "wrap":
        return wrap_signed(xi, 16)
    return sat_signed(xi, 16)


class FixedPointController:
    def __init__(self, sos_int: list[tuple[int, int, int, int, int]]):
        self.sos_int = sos_int
        self.x_hist = [[0, 0] for _ in sos_int]
        self.y_hist = [[0, 0] for _ in sos_int]
        self.product_shift = COEF_FRAC + PIPE_FRAC - ACC_FRAC
        self.pipe_shift = ACC_FRAC - PIPE_FRAC

    def process(self, x_sig: int) -> int:
        pipe = sat_signed(x_sig, 16) << PIPE_SHIFT
        for i, (b0, b1, b2, a1, a2) in enumerate(self.sos_int):
            terms = (
                b0 * pipe,
                b1 * self.x_hist[i][0],
                b2 * self.x_hist[i][1],
                -a1 * self.y_hist[i][0],
                -a2 * self.y_hist[i][1],
            )
            acc = sum(qshift_zero(term, self.product_shift) for term in terms)
            pipe_out = sat_signed(qshift_zero(acc, self.pipe_shift), 32)
            self.x_hist[i][1] = self.x_hist[i][0]
            self.x_hist[i][0] = pipe
            self.y_hist[i][1] = self.y_hist[i][0]
            self.y_hist[i][0] = pipe_out
            pipe = pipe_out
        return sat_signed(qshift_zero(pipe, PIPE_SHIFT), 16)


def controller_float_dc_gain(sos_int: list[tuple[int, int, int, int, int]], scale: int) -> float:
    gain = 1.0
    for b0, b1, b2, a1, a2 in sos_int:
        num = b0 + b1 + b2
        den = scale + a1 + a2
        if den == 0:
            return math.inf if num else math.nan
        gain *= num / den
    return float(gain)


def print_controller_dc_diagnostics(sos_int: list[tuple[int, int, int, int, int]], scale: int, samples: int) -> None:
    print("\nController DC diagnostics from quantized SOS:")
    gain = controller_float_dc_gain(sos_int, scale)
    if np.isfinite(gain):
        print(f"  Quantized-float C(z=1): {gain:+.8g} ({20*np.log10(abs(gain)+1e-300):+.3f} dB)")
    else:
        print(f"  Quantized-float C(z=1): {gain}")

    fragile = []
    for i, (b0, b1, b2, a1, a2) in enumerate(sos_int):
        b_sum = b0 + b1 + b2
        den_sum = scale + a1 + a2
        print(f"  Sec {i}: b0+b1+b2={b_sum:+d}, 1+a1+a2={den_sum:+d}")
        if b_sum == 0 or abs(den_sum) <= 4096:
            fragile.append(i)
    if fragile:
        print(f"  WARNING: sections {fragile} have exact/near DC cancellation.")
        print("           This can make the implemented controller unable to correct a static offset.")

    for step in (1, 10, 100):
        ctrl = FixedPointController(sos_int)
        y = 0
        for _ in range(samples):
            y = ctrl.process(step)
        print(f"  Fixed-point constant input {step:+d} counts -> final output {y:+d} counts")


def reflect_lhp(roots: np.ndarray) -> tuple[np.ndarray, int]:
    roots = np.asarray(roots, dtype=complex).copy()
    count = 0
    for i, root in enumerate(roots):
        if root.real > 0:
            roots[i] = complex(-root.real, root.imag)
            count += 1
    return roots, count


def load_continuous_plant(path: Path, reflect_unstable: bool):
    data = np.load(path, allow_pickle=False)
    z = np.asarray(data["z"], dtype=complex).flatten()
    p = np.asarray(data["p"], dtype=complex).flatten()
    k = float(np.asarray(data["k"]).item())
    reflected_z = reflected_p = 0
    if reflect_unstable:
        z, reflected_z = reflect_lhp(z)
        p, reflected_p = reflect_lhp(p)
    return z, p, k, data.files, reflected_z, reflected_p


def continuous_dc_gain(z, p, k: float) -> float:
    if len(p) and np.any(np.isclose(p, 0.0)):
        return float("nan")
    num = np.prod(-z) if len(z) else 1.0
    den = np.prod(-p) if len(p) else 1.0
    return float(np.real(k * num / den))


def bilinear_zpk(z, p, k: float, fs: float):
    c = 2.0 * fs
    zd = (c + z) / (c - z)
    pd = (c + p) / (c - p)
    relative_degree = max(0, len(p) - len(z))
    if relative_degree:
        zd = np.concatenate([zd, -np.ones(relative_degree)])

    dc = continuous_dc_gain(z, p, k)
    if np.isfinite(dc):
        den = np.prod(1.0 - zd) if len(zd) else 1.0
        num = np.prod(1.0 - pd) if len(pd) else 1.0
        kd = float(np.real(dc * num / den))
    else:
        kd = float(k)
    return zd, pd, kd


def discretize_plant(z, p, k, fs: float, method: str):
    if method == "bilinear":
        zd, pd, kd = bilinear_zpk(z, p, k, fs)
        method_used = "bilinear"
    else:
        try:
            zd, pd, kd, _ = sig.cont2discrete((z, p, k), 1.0 / fs, method=method)
            method_used = method
        except Exception as exc:
            print(f"WARNING: {method} discretization failed ({exc}); falling back to bilinear ZPK mapping.")
            zd, pd, kd = bilinear_zpk(z, p, k, fs)
            method_used = "bilinear-fallback"
    sos = sig.zpk2sos(zd, pd, float(np.real(kd)), pairing="nearest")
    return np.asarray(sos, dtype=float), method_used


class FloatSosPlant:
    def __init__(self, sos: np.ndarray):
        self.sos = np.asarray(sos, dtype=float)
        self.state = np.zeros((len(self.sos), 2), dtype=float)

    def process(self, x: float) -> float:
        y = float(x)
        for i, sec in enumerate(self.sos):
            b0, b1, b2, a0, a1, a2 = sec
            if a0 != 1.0:
                b0, b1, b2, a1, a2 = b0 / a0, b1 / a0, b2 / a0, a1 / a0, a2 / a0
            out = b0 * y + self.state[i, 0]
            self.state[i, 0] = b1 * y - a1 * out + self.state[i, 1]
            self.state[i, 1] = b2 * y - a2 * out
            y = out
        return float(y)

    def set_initial_output(self, target_y: float, first_input: float = 0.0) -> None:
        if len(self.sos) == 0:
            return
        b0 = self.sos[0, 0] / self.sos[0, 3]
        self.state[0, 0] = float(target_y) - b0 * float(first_input)


def normalize_rms(x: np.ndarray, target_rms: float) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    x -= np.mean(x)
    rms = float(np.std(x))
    if rms > 0:
        x *= target_rms / rms
    return x


def generate_acoustic_noise(args) -> np.ndarray:
    if args.no_noise or args.acoustic_rms_counts <= 0:
        return np.zeros(args.samples, dtype=float)

    rng = np.random.default_rng(args.noise_seed)
    freqs = np.fft.rfftfreq(args.samples, d=1.0 / args.fs)
    white = rng.normal(size=len(freqs)) + 1j * rng.normal(size=len(freqs))
    white[0] = rng.normal()
    if len(white) > 1 and args.samples % 2 == 0:
        white[-1] = rng.normal()

    f_floor = max(args.acoustic_low_hz, args.fs / args.samples)
    acoustic_asd = np.sqrt(args.acoustic_corner_hz / np.maximum(freqs, f_floor))
    asd = np.where(freqs <= args.acoustic_corner_hz, acoustic_asd, 0.0)

    spectrum = white * asd
    noise = np.fft.irfft(spectrum, n=args.samples)
    return normalize_rms(noise, args.acoustic_rms_counts)


def generate_sensor_noise(args) -> np.ndarray:
    if args.no_noise or args.sensor_rms_counts <= 0:
        return np.zeros(args.samples, dtype=float)
    rng = np.random.default_rng(args.noise_seed + 1)
    return normalize_rms(rng.normal(size=args.samples), args.sensor_rms_counts)


def initialize_state_for_output(cd: np.ndarray, dd: np.ndarray, u_prev: float, target_y: float) -> np.ndarray:
    n = cd.shape[1]
    x = np.zeros(n, dtype=float)
    c_row = cd.reshape(-1)
    c_norm2 = float(np.dot(c_row, c_row))
    if c_norm2 > 0:
        target = float(target_y) - float(dd.reshape(-1)[0]) * float(u_prev)
        x = c_row * (target / c_norm2)
    return x


def simulate(args, sos_int):
    z, p, k, keys, reflected_z, reflected_p = load_continuous_plant(
        args.plant_npz, args.reflect_unstable_plant
    )
    plant_sos, method_used = discretize_plant(z, p, k, args.fs, args.discretization)
    plant = FloatSosPlant(plant_sos)

    ctrl = FixedPointController(sos_int)
    u_prev = float(args.initial_dac_counts)
    plant.set_initial_output(args.initial_plant_output_counts, u_prev)
    acoustic = generate_acoustic_noise(args)
    sensor = generate_sensor_noise(args)

    adc = np.zeros(args.samples, dtype=np.int16)
    plant_y = np.zeros(args.samples, dtype=float)
    measured_y = np.zeros(args.samples, dtype=float)
    ctrl_y = np.zeros(args.samples, dtype=np.int16)
    dac = np.zeros(args.samples, dtype=np.int16)
    diverged_at = None

    for n in range(args.samples):
        y = float(plant.process(u_prev) + args.disturbance_counts + acoustic[n])
        y_meas = y + sensor[n]
        if (
            not np.isfinite(y_meas)
            or abs(y_meas) > args.diverge_limit_counts
            or not np.all(np.isfinite(plant.state))
            or np.max(np.abs(plant.state)) > args.diverge_limit_counts
        ):
            diverged_at = n
            break
        adc_in = cast_adc16(y_meas + args.adc_offset_counts, args.adc_cast)
        c_out = ctrl.process(adc_in)
        dac_cmd = sat_dac14(args.controller_sign * c_out + args.dac_offset_counts)

        adc[n] = adc_in
        plant_y[n] = y
        measured_y[n] = y_meas
        ctrl_y[n] = c_out
        dac[n] = dac_cmd

        u_prev = float(dac_cmd)

    end = args.samples if diverged_at is None else max(diverged_at, 0)

    return {
        "plant_keys": keys,
        "plant_order": len(p),
        "plant_zeros": len(z),
        "unstable_poles": int(np.count_nonzero(np.real(p) > 0)),
        "unstable_zeros": int(np.count_nonzero(np.real(z) > 0)),
        "reflected_poles": reflected_p,
        "reflected_zeros": reflected_z,
        "plant_discretization": method_used,
        "diverged_at": diverged_at,
        "adc": adc[:end],
        "plant_y": plant_y[:end],
        "measured_y": measured_y[:end],
        "acoustic": acoustic[:end],
        "sensor": sensor[:end],
        "ctrl": ctrl_y[:end],
        "dac": dac[:end],
    }


def summarize_trace(name: str, x: np.ndarray, tail: int) -> None:
    if len(x) == 0:
        print(f"  {name:8s}: no finite samples")
        return
    tail_x = x[-min(tail, len(x)) :]
    print(
        f"  {name:8s}: final={x[-1]:+.6g}, tail mean={np.mean(tail_x):+.6g}, "
        f"tail std={np.std(tail_x):.6g}, min={np.min(x):+.6g}, max={np.max(x):+.6g}"
    )


def write_csv(path: Path, fs: float, result: dict[str, np.ndarray]) -> None:
    with path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "n",
            "t_s",
            "adc_counts",
            "physical_plant_output_counts",
            "measured_before_adc_offset_counts",
            "acoustic_disturbance_counts",
            "sensor_noise_counts",
            "controller_counts",
            "dac_counts",
        ])
        for n, (adc, plant_y, measured_y, acoustic, sensor, ctrl, dac) in enumerate(
            zip(
                result["adc"],
                result["plant_y"],
                result["measured_y"],
                result["acoustic"],
                result["sensor"],
                result["ctrl"],
                result["dac"],
            )
        ):
            w.writerow([
                n,
                n / fs,
                int(adc),
                float(plant_y),
                float(measured_y),
                float(acoustic),
                float(sensor),
                int(ctrl),
                int(dac),
            ])


def plot_result(path: Path, fs: float, result: dict[str, np.ndarray], args) -> None:
    if len(result["adc"]) == 0:
        print("Skipping plot: no finite samples")
        return
    n = np.arange(len(result["adc"]))
    t_ms = 1e3 * n / fs
    fig, axes = plt.subplots(3, 1, figsize=(11, 8), sharex=True)
    axes[0].plot(t_ms, result["adc"], lw=1.2, label="ADC input to controller")
    axes[0].plot(t_ms, result["plant_y"], lw=0.9, alpha=0.75, label="physical plant output + acoustic")
    axes[0].plot(t_ms, result["measured_y"], lw=0.7, alpha=0.55, label="measured before ADC offset")
    axes[0].plot(t_ms, result["acoustic"], lw=0.6, alpha=0.45, label="acoustic disturbance")
    axes[0].plot(t_ms, result["sensor"], lw=0.5, alpha=0.35, label="sensor white noise")
    axes[1].plot(t_ms, result["ctrl"], lw=1.2, color="C1", label="controller output")
    axes[2].plot(t_ms, result["dac"], lw=1.2, color="C3", label="DAC after offset/saturation")
    for ax in axes:
        ax.grid(True, alpha=0.3)
        ax.legend(loc="best")
        ax.set_ylabel("counts")
    axes[-1].set_xlabel("time (ms)")
    fig.suptitle(
        f"Closed-loop fixed-point sim | fs={fs/1e3:.3f} kHz | "
        f"adc_offset={args.adc_offset_counts}, dac_offset={args.dac_offset_counts}, sign={args.controller_sign:+d}"
    )
    fig.tight_layout()
    fig.savefig(path, dpi=130)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("plant_npz", type=Path, help="continuous-time fitted plant NPZ")
    parser.add_argument("--coeffs", type=Path, default=ROOT / "src" / "hinf_coeffs.h")
    parser.add_argument("--fs", type=float, default=125e6 / 128.0, help="controller update rate in Hz")
    parser.add_argument("--samples", type=int, default=100_000)
    parser.add_argument("--dc-test-samples", type=int, default=50_000)
    parser.add_argument("--discretization", choices=("zoh", "bilinear", "euler", "backward_diff"), default="bilinear")
    parser.add_argument("--reflect-unstable-plant", action="store_true",
                        help="mirror fitted RHP plant roots into the LHP before time-domain simulation")
    parser.add_argument("--controller-sign", type=int, choices=(-1, 1), default=1)
    parser.add_argument("--adc-offset-counts", type=float, default=0.0)
    parser.add_argument("--dac-offset-counts", type=float, default=0.0)
    parser.add_argument("--disturbance-counts", type=float, default=0.0, help="constant plant-output bias before ADC offset")
    parser.add_argument("--acoustic-rms-counts", type=float, default=None,
                        help="RMS of physical 1/f acoustic output disturbance in ADC counts")
    parser.add_argument("--noise-rms-counts", type=float, default=None,
                        help="deprecated alias for --acoustic-rms-counts")
    parser.add_argument("--sensor-rms-counts", type=float, default=1.0,
                        help="RMS of white sensor/readout noise in ADC counts")
    parser.add_argument("--no-noise", action="store_true", help="disable injected acoustic/readout noise")
    parser.add_argument("--noise-seed", type=int, default=1)
    parser.add_argument("--acoustic-corner-hz", type=float, default=5e3,
                        help="1/f acoustic ASD extends up to this frequency; white above")
    parser.add_argument("--acoustic-low-hz", type=float, default=1.0,
                        help="low-frequency floor for finite 1/f noise generation")
    parser.add_argument("--initial-plant-output-counts", type=float, default=0.0)
    parser.add_argument("--initial-dac-counts", type=float, default=0.0)
    parser.add_argument("--adc-cast", choices=("saturate", "wrap"), default="saturate")
    parser.add_argument("--tail", type=int, default=10_000)
    parser.add_argument("--diverge-limit-counts", type=float, default=1e12,
                        help="stop simulation when plant output/state exceeds this count scale")
    parser.add_argument("--csv", type=Path)
    parser.add_argument("--plot", type=Path, default=Path("closed_loop_fixedpoint_sim.png"))
    parser.add_argument("--no-plot", action="store_true")
    args = parser.parse_args()
    if args.acoustic_rms_counts is None:
        args.acoustic_rms_counts = 5.0 if args.noise_rms_counts is None else args.noise_rms_counts

    sos_int, scale, int_bits, frac_bits = parse_hinf_coeffs(args.coeffs)
    print(f"Loaded controller: {len(sos_int)} SOS sections from {args.coeffs} (Q{int_bits}.{frac_bits})")
    print_controller_dc_diagnostics(sos_int, scale, args.dc_test_samples)

    result = simulate(args, sos_int)
    print(f"\nLoaded plant: {args.plant_npz} ({result['plant_zeros']} zeros, {result['plant_order']} poles)")
    print(f"Plant discretization: {result['plant_discretization']}")
    if result["reflected_poles"] or result["reflected_zeros"]:
        print(
            f"Reflected fitted RHP roots to LHP: {result['reflected_zeros']} zeros, "
            f"{result['reflected_poles']} poles"
        )
    elif result["unstable_poles"]:
        print(
            f"WARNING: fitted plant has {result['unstable_poles']} RHP pole(s). "
            "Time-domain simulation may diverge; try --reflect-unstable-plant for a physical-stable approximation."
        )
    if result["diverged_at"] is None:
        print(f"Simulated {args.samples} controller samples ({args.samples / args.fs:.6g} s)")
    else:
        print(
            f"Simulation diverged at sample {result['diverged_at']} "
            f"({result['diverged_at'] / args.fs:.6g} s); reporting finite prefix."
        )
    if args.no_noise:
        print("Noise model: disabled")
    else:
        print(
            f"Noise model: acoustic/output disturbance {args.acoustic_rms_counts:g} counts RMS "
            f"with 1/f ASD to {args.acoustic_corner_hz:g} Hz; "
            f"sensor/readout white {args.sensor_rms_counts:g} counts RMS; seed={args.noise_seed}"
        )
    print("\nClosed-loop trace summary:")
    summarize_trace("ADC", result["adc"], args.tail)
    summarize_trace("plant", result["plant_y"], args.tail)
    summarize_trace("measured", result["measured_y"], args.tail)
    summarize_trace("acoustic", result["acoustic"], args.tail)
    summarize_trace("sensor", result["sensor"], args.tail)
    summarize_trace("ctrl", result["ctrl"], args.tail)
    summarize_trace("DAC", result["dac"], args.tail)

    if args.csv:
        write_csv(args.csv, args.fs, result)
        print(f"Wrote CSV: {args.csv}")
    if not args.no_plot:
        plot_result(args.plot, args.fs, result, args)
        print(f"Saved plot: {args.plot}")


if __name__ == "__main__":
    main()
