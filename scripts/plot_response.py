#!/usr/bin/env python3
"""
Plot testbench CSV outputs (step / impulse / frequency sweep).

Usage: python3 scripts/plot_response.py [sim_output_dir]
"""
import os
import sys
import csv
import matplotlib.pyplot as plt

DEFAULT_SIM_DIR = "build/sim_workspace/dy_cavity_relocker_2_csim/solution1/csim/build"


def read_csv(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def plot_frequency(sim_dir):
    path = os.path.join(sim_dir, "freq_response.csv")
    if not os.path.exists(path):
        return
    rows = read_csv(path)
    f_hz   = [float(r["freq_hz"])   for r in rows]
    gain   = [float(r["gain_db"])   for r in rows]
    phase  = [float(r["phase_deg"]) for r in rows]

    dc = gain[0]
    bw = next((f_hz[i] for i, g in enumerate(gain) if g < dc - 3.0), None)
    print(f"Frequency response: DC={dc:.1f} dB, -3dB BW ≈ "
          f"{bw/1e3 if bw else float('nan'):.1f} kHz")

    fig, (a1, a2) = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    a1.semilogx(f_hz, gain, lw=1.5)
    a1.axhline(dc - 3, color="red", ls="--", alpha=0.5, label="-3 dB")
    if bw: a1.axvline(bw, color="red", ls="--", alpha=0.5)
    a1.set_ylabel("Gain (dB)")
    a1.set_title("Controller frequency response")
    a1.grid(True, which="both", alpha=0.3); a1.legend()

    a2.semilogx(f_hz, phase, color="C3", lw=1.5)
    a2.set_xlabel("Frequency (Hz)"); a2.set_ylabel("Phase (deg)")
    a2.grid(True, which="both", alpha=0.3)

    plt.tight_layout()
    out = os.path.join(sim_dir, "freq_response.png")
    plt.savefig(out, dpi=150)
    print(f"  → {out}")


def plot_time(sim_dir, csv_name, title, ref_line=None):
    path = os.path.join(sim_dir, csv_name)
    if not os.path.exists(path):
        return
    rows = read_csv(path)
    n = [int(r["frame"])         for r in rows]
    y = [int(r["output_ch1"])    for r in rows]

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(n, y, lw=1.0)
    if ref_line is not None:
        ax.axhline(ref_line, color="red", ls="--", alpha=0.5, label=f"ref={ref_line}")
        ax.legend()
    ax.set_xlabel("Frame"); ax.set_ylabel("Output ch1")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)

    out = os.path.join(sim_dir, csv_name.replace(".csv", ".png"))
    plt.savefig(out, dpi=150)
    print(f"  → {out}")


def main():
    sim_dir = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_SIM_DIR
    if not os.path.isdir(sim_dir):
        print(f"ERROR: {sim_dir} does not exist", file=sys.stderr)
        sys.exit(1)

    plot_frequency(sim_dir)
    plot_time(sim_dir, "step_response.csv",    "Step response (input=50)")
    plot_time(sim_dir, "impulse_response.csv", "Impulse response")


if __name__ == "__main__":
    main()
