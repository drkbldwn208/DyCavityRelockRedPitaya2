#!/usr/bin/env python3
"""
Plot the closed-loop acquisition from the HLS C-simulation.
Overlays the float-controller and fixed-point-controller runs if both exist.

Usage: python3 scripts/plot_closed_loop.py [sim_output_dir]
"""
import os
import sys
import csv
import matplotlib.pyplot as plt

DEFAULT_SIM_DIR = "build/sim_workspace/dy_cavity_relocker_2_csim/solution1/csim/build"


def read_run(path):
    if not os.path.exists(path):
        return None
    samples, errors, efforts = [], [], []
    with open(path) as f:
        for r in csv.DictReader(f):
            samples.append(int(r["sample"]))
            errors.append(float(r["error"]))
            efforts.append(float(r["effort"]))
    return samples, errors, efforts


def main():
    sim_dir = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_SIM_DIR

    runs = {
        "float":      read_run(os.path.join(sim_dir, "closed_loop_float.csv")),
        "fixed-point": read_run(os.path.join(sim_dir, "closed_loop_fixed.csv")),
    }
    runs = {k: v for k, v in runs.items() if v is not None}
    if not runs:
        print(f"No closed_loop_{{float,fixed}}.csv found in {sim_dir}", file=sys.stderr)
        sys.exit(1)

    fig, (a1, a2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    fig.suptitle("H∞ controller: closed-loop acquisition")

    colors = {"float": "#1f77b4", "fixed-point": "#d62728"}
    for label, (s, e, u) in runs.items():
        a1.plot(s, e, color=colors[label], lw=1.4, label=label, alpha=0.85)
        a2.plot(s, u, color=colors[label], lw=1.4, label=label, alpha=0.85)

    a1.axhline(0, color="black", ls=":", alpha=0.5)
    a1.set_ylabel("Residual error")
    a1.set_title("Error  (target: 0)")
    a1.grid(True, alpha=0.3); a1.legend()

    a2.set_xlabel("Filter sample  (~1 µs each)")
    a2.set_ylabel("Control effort")
    a2.set_title("DAC drive")
    a2.grid(True, alpha=0.3); a2.legend()

    plt.tight_layout()
    out = os.path.join(sim_dir, "closed_loop.png")
    plt.savefig(out, dpi=150)
    print(f"Saved {out}")


if __name__ == "__main__":
    main()
