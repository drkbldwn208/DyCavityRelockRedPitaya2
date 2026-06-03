import numpy as np
import control as ct
import matplotlib.pyplot as plt
from scipy.signal import tf2zpk

def load_plant_zpk(filename):
    """Loads zeros, poles, and gain from the .npz file."""
    try:
        data = np.load(filename)
        z = data['z']
        p = data['p']
        k = data['k']
        return z, p, k
    except Exception as e:
        print(f"Error loading {filename}: {e}")
        return None, None, None

def main():
    # 1. Load the measured plant data
    filename = 'plant_only.npz'
    z_plant, p_plant, k_plant = load_plant_zpk(filename)

    if z_plant is None:
        return

    # Convert the raw ZPK data into a control.TransferFunction object
    plant_sys = ct.zpk(z_plant, p_plant, k_plant)
    print(f"--- Measured Plant ---")
    print(plant_sys)

    # 2. Define the Controller Parameters
    Kp = 0.0       # Set to 0.0 to test the Pure Integrator mode
    Ki = 51.0
    epsilon = 1e-2 # Very slow pole to safely approximate a pure integrator

    # 3. Formulate the ZPK based on Kp
    if Kp == 0:
        print("\n[Mode: Pure Integrator]")
        # C(s) = Ki / (s + epsilon)
        z_ctrl = []            # No zeros
        p_ctrl = [-epsilon]    # Leaky pole at origin
        k_ctrl = Ki            # Gain is exactly Ki
    else:
        print("\n[Mode: Proportional-Integral]")
        # C(s) = Kp * (s + Ki/Kp) / (s + epsilon)
        z_ctrl = [-Ki / Kp]    # Zero at -Ki/Kp
        p_ctrl = [-epsilon]    # Leaky pole at origin
        k_ctrl = Kp            # High frequency gain is Kp

    # Build the controller system object
    ctrl_sys = ct.zpk(z_ctrl, p_ctrl, k_ctrl)
    print(f"--- Controller Model ---")
    print(ctrl_sys)

    # 4. Create the Effective Plant for Synthesis
    # Multiply the transfer functions: P_eff(s) = C(s) * P(s)
    effective_sys = ct.series(ctrl_sys, plant_sys)
    print(f"\n--- Effective Plant (Controller * Plant) ---")
    print(effective_sys)

    num = effective_sys.num[0][0]
    den = effective_sys.den[0][0]

    z_eff, p_eff, k_eff = tf2zpk(num, den)
    zeros = np.asarray(z_eff)
    poles = np.asarray(p_eff)
    gain = np.asarray(float(np.real(k_eff)))

    output_filename = "plant_with_controller.npz"
    np.savez(output_filename, z=zeros, p=poles, k=gain)
    print(f"\nSaved effective plant ZPK to {output_filename}")

    # 5. Visualize the results
    plt.figure(figsize=(10, 8))

    # Plot Bode for both the raw plant and the effective plant
    ct.bode_plot([plant_sys, effective_sys],
                 label=['Bare Plant', f'Effective Plant (Kp={Kp}, Ki={Ki})'],
                 Hz=True, deg=True, margins=False)

    plt.suptitle("Bode Plot: Bare Plant vs. Effective Plant")
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()