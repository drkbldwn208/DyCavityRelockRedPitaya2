#!/usr/bin/env python3
"""
Plot the closed-loop acquisition from the HLS C-simulation.

Usage: python3 scripts/plot_closed_loop.py [sim_output_dir]
  Default sim_output_dir: build/sim_workspace/dy_cavity_relocker_2_csim/solution1/csim/build
"""
import sys
import os
import csv

try:
    import matplotlib
    # matplotlib.use('Agg') # Uncomment if you are running entirely headless and plt.show() crashes
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False
    print("Error: matplotlib is required to generate the plots.")
    sys.exit(1)

def main():
    if len(sys.argv) > 1:
        sim_dir = sys.argv[1]
    else:
        sim_dir = "build/sim_workspace/dy_cavity_relocker_2_csim/solution1/csim/build"

    csv_path = os.path.join(sim_dir, "closed_loop.csv")

    if not os.path.exists(csv_path):
        print(f"ERROR: Could not find {csv_path}")
        print("Make sure you ran the C-simulation with Test 5 enabled.")
        sys.exit(1)

    samples = []
    errors = []
    efforts = []

    print(f"Reading data from {csv_path}...")
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            samples.append(int(row['sample']))
            errors.append(float(row['adc_error']))
            efforts.append(int(row['dac_effort']))

    if not samples:
        print("ERROR: CSV file is empty.")
        sys.exit(1)

    print("Generating plot...")
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    fig.suptitle('H∞ Controller: Closed-Loop Lock Acquisition', fontsize=14)

    # Top Plot: Error Signal
    ax1.plot(samples, errors, color='#d62728', linewidth=2.0, label='Residual Error (ADC)')
    ax1.axhline(0, color='black', linestyle=':', alpha=0.6)
    ax1.set_ylabel('Cavity Error (ADU)')
    ax1.set_title('Error Signal (Target: 0.0)')
    ax1.grid(True, which='both', alpha=0.3)
    ax1.legend(loc='upper right')

    # Bottom Plot: Control Effort
    ax2.plot(samples, efforts, color='#1f77b4', linewidth=2.0, label='Drive Signal (DAC)')
    ax2.set_xlabel('Time (Samples / µs)')
    ax2.set_ylabel('Control Effort (ADU)')
    ax2.set_title('Piezo Actuator Drive')
    ax2.grid(True, which='both', alpha=0.3)
    ax2.legend(loc='lower right')

    plt.tight_layout()
    
    out_path = os.path.join(sim_dir, "closed_loop_acquisition.png")
    plt.savefig(out_path, dpi=150)
    print(f"\nSuccess! Plot saved to: {out_path}")

    # Display the plot if running in an interactive GUI environment
    try:
        plt.show()
    except Exception as e:
        print("Note: Could not display interactive window (likely headless environment).")

if __name__ == '__main__':
    main()