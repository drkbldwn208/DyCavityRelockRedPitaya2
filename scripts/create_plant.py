"""
create_plant.py

Generates a continuous-time ZPK (Zero-Pole-Gain) .npz file 
for use with ctrl.py via the --plant-npz flag.
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig

# =========================================================================
# HELPER FUNCTIONS
# =========================================================================

def real_root(f_hz):
    """Returns a single real root at the specified frequency in Hz."""
    return [-abs(f_hz)]

def complex_pair(fn_hz, Q):
    """
    Returns a pair of roots for a 2nd-order system given natural frequency and Q.
    Handles both underdamped (complex conjugate) and overdamped (two real) cases.
    """
    zeta = 1.0 / (2.0 * Q)
    if zeta < 1.0:
        # Underdamped: Complex conjugate pair
        wd = fn_hz * np.sqrt(1.0 - zeta**2)
        sigma = -zeta * fn_hz
        return [sigma + 1j*wd, sigma - 1j*wd]
    else:
        # Overdamped: Two real roots
        w1 = fn_hz * (zeta - np.sqrt(zeta**2 - 1.0))
        w2 = fn_hz * (zeta + np.sqrt(zeta**2 - 1.0))
        return [-w1, -w2]

# =========================================================================
# USER CONFIGURATION
# Construct your plant by adding lists of roots together.
# =========================================================================

ZEROS_HZ = sum([
    # Add your zeros here
    # complex_pair(2.6e3, Q=15.0),
    # complex_pair(4.8e3, Q=2),
    # complex_pair(70e3, Q=0.4),
    # complex_pair(13e3, Q=25),
    complex_pair(42134.4, Q=1)
], [])

POLES_HZ = sum([
    # Add your poles here
    complex_pair(2744.1, Q=3.11)   # Example: Overdamped 2nd order roll-off
], [])

# Set the desired DC gain of the plant
TARGET_DC_GAIN = 10.0
OUTPUT_FILE = "custom_plant.npz"

# =========================================================================
# GENERATION MATH
# =========================================================================

def main():
    print(f"Building custom plant: {len(ZEROS_HZ)} zeros, {len(POLES_HZ)} poles.")
    
    # 1. Convert from Hz to rad/s
    z_rad = np.array(ZEROS_HZ) * 2 * np.pi
    p_rad = np.array(POLES_HZ) * 2 * np.pi
    
    # 2. Calculate the raw gain 'k' required to hit the TARGET_DC_GAIN
    prod_z = np.prod(-z_rad) if len(z_rad) > 0 else 1.0
    prod_p = np.prod(-p_rad) if len(p_rad) > 0 else 1.0
    
    if np.isclose(np.abs(prod_p), 0):
        print("Warning: Pole at origin detected. DC gain is infinite. Setting k=1.0.")
        k_raw = 1.0
    else:
        k_raw = np.real(TARGET_DC_GAIN * prod_p / prod_z)
        
    print(f"Target DC Gain: {TARGET_DC_GAIN}")
    print(f"Calculated raw k: {k_raw:.5e}")

    # 3. Save to .npz
    np.savez(OUTPUT_FILE, z=z_rad, p=p_rad, k=k_raw)
    print(f"Saved to {OUTPUT_FILE}\n")

    # 4. Plot verification Bode
    w = np.logspace(1, 6, 2000) * 2 * np.pi
    w_hz = w / (2 * np.pi)
    
    _, mag, phase = sig.bode((z_rad, p_rad, k_raw), w=w)
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    fig.suptitle(f"Custom Plant Verification\nDC Gain = {TARGET_DC_GAIN:.2f}")
    
    ax1.semilogx(w_hz, mag, lw=2, color='C3')
    ax1.set_ylabel("Magnitude (dB)")
    ax1.grid(True, which="both", alpha=0.3)
    ax1.axhline(20*np.log10(TARGET_DC_GAIN), color='gray', ls=':', lw=1.5, label="Target DC")
    ax1.legend()
    
    ax2.semilogx(w_hz, phase, lw=2, color='C3')
    ax2.set_ylabel("Phase (deg)")
    ax2.set_xlabel("Frequency (Hz)")
    ax2.grid(True, which="both", alpha=0.3)
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()