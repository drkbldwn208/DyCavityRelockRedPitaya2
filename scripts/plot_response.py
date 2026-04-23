#!/usr/bin/env python3
"""
Plot frequency response, impulse response, and step response
from the HLS C-simulation CSV outputs.

Usage: python3 scripts/plot_response.py [sim_output_dir]
  Default sim_output_dir: build/sim_workspace/dy_cavity_relocker_2_csim/solution1/csim/build
"""
import sys
import os
import csv
import math

def read_csv(path):
    """Read CSV into list of dicts."""
    with open(path, 'r') as f:
        reader = csv.DictReader(f)
        return list(reader)

def main():
    if len(sys.argv) > 1:
        sim_dir = sys.argv[1]
    else:
        sim_dir = "build/sim_workspace/dy_cavity_relocker_2_csim/solution1/csim/build"

    freq_csv = os.path.join(sim_dir, "freq_response.csv")
    imp_csv = os.path.join(sim_dir, "impulse_response.csv")
    step_csv = os.path.join(sim_dir, "step_response.csv")

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        HAS_MPL = True
    except ImportError:
        HAS_MPL = False
        print("matplotlib not found — printing text summary only")

    # ---- Frequency Response ----
    if os.path.exists(freq_csv):
        data = read_csv(freq_csv)
        freqs_mhz = [float(r['freq_mhz']) for r in data]
        gains_db   = [float(r['gain_db'])  for r in data]
        phases_deg = [float(r['phase_deg']) for r in data]

        print("\n=== Frequency Response ===")
        print(f"{'Freq (MHz)':>12}  {'Gain (dB)':>10}  {'Phase (deg)':>12}")
        for i in range(0, len(data), max(1, len(data)//20)):
            print(f"{freqs_mhz[i]:12.4f}  {gains_db[i]:10.2f}  {phases_deg[i]:12.1f}")

        # Find -3dB point
        dc_gain = gains_db[0] if gains_db else 0
        bw_freq = None
        for i, g in enumerate(gains_db):
            if g < dc_gain - 3.0:
                bw_freq = freqs_mhz[i]
                break

        if bw_freq:
            print(f"\nEstimated -3dB bandwidth: {bw_freq:.4f} MHz ({bw_freq*1e3:.1f} kHz)")
        print(f"DC gain: {dc_gain:.2f} dB")
        print(f"Gain at Nyquist ({freqs_mhz[-1]:.2f} MHz): {gains_db[-1]:.2f} dB")

        if HAS_MPL:
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

            ax1.semilogx(freqs_mhz, gains_db, 'b-', linewidth=1.5)
            ax1.set_ylabel('Gain (dB)')
            ax1.set_title('H-infinity Filter Frequency Response (HLS C-Sim)')
            ax1.grid(True, which='both', alpha=0.3)
            ax1.axhline(y=dc_gain - 3, color='r', linestyle='--', alpha=0.5, label='-3dB')
            if bw_freq:
                ax1.axvline(x=bw_freq, color='r', linestyle='--', alpha=0.5)
            ax1.legend()

            ax2.semilogx(freqs_mhz, phases_deg, 'r-', linewidth=1.5)
            ax2.set_ylabel('Phase (deg)')
            ax2.set_xlabel('Frequency (MHz)')
            ax2.grid(True, which='both', alpha=0.3)

            plt.tight_layout()
            out_path = os.path.join(sim_dir, "freq_response.png")
            plt.savefig(out_path, dpi=150)
            print(f"\nBode plot saved to {out_path}")
    else:
        print(f"WARNING: {freq_csv} not found")

    # ---- Impulse Response ----
    if os.path.exists(imp_csv):
        data = read_csv(imp_csv)
        samples = [int(r['frame']) for r in data]
        inp     = [int(r['input']) for r in data]
        out_ch1 = [int(r['output_ch1']) for r in data]

        print("\n=== Impulse Response (first 40 nonzero-region samples) ===")
        for i in range(min(40, len(data))):
            print(f"  n={samples[i]:4d}  in={inp[i]:6d}  out_ch1={out_ch1[i]:6d}")

        if HAS_MPL:
            fig, ax = plt.subplots(figsize=(10, 4))
            ax.plot(samples[:200], out_ch1[:200], 'b-', linewidth=1)
            ax.set_xlabel('Sample')
            ax.set_ylabel('Output (ch1)')
            ax.set_title('Impulse Response')
            ax.grid(True, alpha=0.3)
            out_path = os.path.join(sim_dir, "impulse_response.png")
            plt.savefig(out_path, dpi=150)
            print(f"Impulse response plot saved to {out_path}")

    # ---- Step Response ----
    if os.path.exists(step_csv):
        data = read_csv(step_csv)
        samples = [int(r['frame']) for r in data]
        out_ch1 = [int(r['output_ch1']) for r in data]

        print("\n=== Step Response (first 20 and last 20 samples) ===")
        for i in range(min(20, len(data))):
            print(f"  n={samples[i]:4d}  out_ch1={out_ch1[i]:6d}")
        print("  ...")
        for i in range(max(0, len(data)-20), len(data)):
            print(f"  n={samples[i]:4d}  out_ch1={out_ch1[i]:6d}")

        if HAS_MPL:
            fig, ax = plt.subplots(figsize=(10, 4))
            ax.plot(samples, out_ch1, 'b-', linewidth=1)
            ax.axhline(y=2000, color='r', linestyle='--', alpha=0.5, label='Input (2000)')
            ax.set_xlabel('Sample')
            ax.set_ylabel('Output (ch1)')
            ax.set_title('Step Response')
            ax.grid(True, alpha=0.3)
            ax.legend()
            out_path = os.path.join(sim_dir, "step_response.png")
            plt.savefig(out_path, dpi=150)
            print(f"Step response plot saved to {out_path}")

    print("\nDone.")

if __name__ == '__main__':
    main()