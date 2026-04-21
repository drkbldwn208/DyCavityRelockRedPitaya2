"""
H-infinity synthesis — v12

Goal: modest but real improvement over existing loop.
  - 20-40 dB better suppression from DC to 5 kHz (factor 10-100x)
  - Crossover stays near 30 kHz
  - No S+T=1 conflict (W3 DC is tiny, W1_MID is small)
  - gamma should land 0.5-2.0 — feasible, well-posed

The -100 dB result in v11 came from W1_DC=500 being wildly over-specified.
W1_DC=50 gives 34 dB improvement over existing which is already excellent
for a cavity intensity lock.

Primary knob: W1_DC
  10  → 20 dB better than existing
  30  → 30 dB better
  100 → 40 dB better
  500 → overconstrained, gamma explodes (what we had in v11)
"""

import control as ct
import matplotlib.pyplot as plt
import numpy as np
import time

# ---------------------------------------------------------------------------
# NORMALIZATION
# ---------------------------------------------------------------------------
f_xover = 15e3
w_norm  = 2 * np.pi * f_xover

def w(f_hz):
    return 2 * np.pi * f_hz / w_norm

def makeweight(dcgain, wc_norm, hfgain):
    if dcgain > hfgain:
        M, A = dcgain, hfgain / dcgain
        return ct.tf([1.0/M, wc_norm], [1.0, wc_norm * A])
    else:
        M, A = hfgain, dcgain / hfgain
        return ct.tf([1.0, wc_norm * A], [1.0/M, wc_norm])

# ===========================================================================
# TUNING KNOBS
# ===========================================================================
F_POLE1_HZ = 30e3
F_POLE2_HZ = 50e3
PLANT_HAS_INTEGRATOR = False

# W1_DC is your primary knob — this is all you should need to change.
# Existing loop has ~-35 dB at 100-500 Hz.
# W1_DC=50 demands |S| < 1/50 below 5 kHz.
# Combined with existing ~-35 dB that gives ~-69 dB total — factor ~30x better.
# Raise W1_DC until gamma > 2, then you're near the plant's limit.
W1_DC         = 50.0   # <-- primary knob. Try 10, 30, 50, 100.
W1_CORNER     = 8e3    # Hz — edge of suppression band
W1_MID        = 3.0    # gentle anchor: keeps crossover near 25-30 kHz
                       # without conflicting with W3. Keep at 2-5.
W1_MID_CORNER = 15e3   # Hz
W1_HF         = 0.5

# W3: HF guardrail ONLY. Tiny DC to avoid S+T=1 conflict with W1.
W3_DC     = 0.05
W3_CORNER = 20e3
W3_HF     = 0.3

# W2: suppress K above actuator bandwidth
W2_DC     = 1e-3
W2_CORNER = 80e3
W2_HF     = 1.0

TAU_S     = 0.5e-6
# ===========================================================================

# ---------------------------------------------------------------------------
# 1. PLANT
# ---------------------------------------------------------------------------
# wp1 = w(F_POLE1_HZ)
# wp2 = w(F_POLE2_HZ)
# den_poles = np.polymul([1.0/wp1, 1.0], [1.0/wp2, 1.0])
# den_G = np.polymul(den_poles, [1.0, 0.0]) if PLANT_HAS_INTEGRATOR else den_poles
# G_raw = ct.tf([1.0], den_G)

# G_mag = float(np.abs(ct.frequency_response(G_raw,
#               np.array([w(f_xover)])).frdata.squeeze()))
# G = G_raw * (1.0 / G_mag)

# K_plant = 3533.3
# w_leak = 2*np.pi*1.0
# norm_K = K_plant / w_norm
# norm_w_leak = w_leak / w_norm
# G = ct.tf([norm_K], [1.0, norm_w_leak])

K_PLANT = 3533.3 # rad/s

# 1a. Build the Core Leaky Integrator
w_leak = 2 * np.pi * 1.0
norm_K = K_PLANT / w_norm
norm_w_leak = w_leak / w_norm
G_core = ct.tf([norm_K], [1.0, norm_w_leak])

# 1b. Helper Function for Complex Pole/Zero Resonance Pairs
def make_complex_pz_pair(p_real_hz, p_imag_hz, z_real_hz, z_imag_hz):
    """
    Creates a normalized second-order pole/zero pair from complex roots.
    - p_real_hz: Real part of the pole (Must be NEGATIVE for stability, e.g., -80)
    - p_imag_hz: Imaginary part of the pole (The peak frequency, e.g., 2500)
    - z_real_hz: Real part of the zero (e.g., -80)
    - z_imag_hz: Imaginary part of the zero (The notch frequency, e.g., 2600)
    """
    # Convert inputs to normalized angular frequencies
    pr_n = 2 * np.pi * p_real_hz / w_norm
    pi_n = 2 * np.pi * p_imag_hz / w_norm
    zr_n = 2 * np.pi * z_real_hz / w_norm
    zi_n = 2 * np.pi * z_imag_hz / w_norm

    # Calculate the squared magnitudes: |p|^2 = (Real^2 + Imag^2)
    p_mag2 = pr_n**2 + pi_n**2
    z_mag2 = zr_n**2 + zi_n**2

    # DC Correction Factor: Evaluated at s=0, the gain is z_mag2 / p_mag2.
    # We multiply the numerator by the inverse to force DC gain to 1.0 (0 dB).
    dc_correction = p_mag2 / z_mag2

    # Polynomial expansion: (s - p)(s - p*) = s^2 - 2*Real(p)*s + |p|^2
    # Numerator: DC_corr * (s^2 - 2*Real(z)*s + |z|^2)
    num = [dc_correction, dc_correction * (-2 * zr_n), dc_correction * z_mag2]

    # Denominator: s^2 - 2*Real(p)*s + |p|^2
    den = [1.0, -2 * pr_n, p_mag2]

    return ct.tf(num, den)

# 1c. Define Your Physical Table Resonances via Complex Roots
# Example:
# Pole at 2500 Hz. Real part of -80 Hz means a resonance width of ~160 Hz.
# Zero at 2600 Hz. Real part of -80 Hz gives the anti-resonance notch a similar width.

R1 = make_complex_pz_pair(p_real_hz=-6.64, p_imag_hz=2549.7, z_real_hz=-9.05, z_imag_hz=2558.4)
R2 = make_complex_pz_pair(p_real_hz=-19.5, p_imag_hz=5641.1, z_real_hz=-20.9, z_imag_hz=5672.7)
R3 = make_complex_pz_pair(p_real_hz=-17.2, p_imag_hz=4058.1, z_real_hz=-10.5, z_imag_hz=4089.8)
R4 = make_complex_pz_pair(p_real_hz=-1e4, p_imag_hz=3400, z_real_hz=-1e4, z_imag_hz=2e4)
R4 = make_complex_pz_pair(p_real_hz=-1e4, p_imag_hz=9e2, z_real_hz=-5e2, z_imag_hz=1.6e2)

# 1d. Assemble the Fully Dressed Plant
# Cascade the integrator with the complex resonances
G = G_core * R1 * R2 * R3 * R4

try:
    _, pm0, _, wpc0 = ct.margin(G)
    print(f"Bare plant:  xover={float(wpc0)*w_norm/(2*np.pi)/1e3:.1f} kHz  PM={pm0:.1f} deg")
except Exception:
    print("Bare plant: no finite crossover")

# ---------------------------------------------------------------------------
# 2. DELAY — first-order lag
# ---------------------------------------------------------------------------
Delay = ct.tf([1.0], [TAU_S * w_norm, 1.0])
G_aug = G * Delay

try:
    _, pm_d, _, wpc_d = ct.margin(G_aug)
    print(f"Plant+delay: xover={float(wpc_d)*w_norm/(2*np.pi)/1e3:.1f} kHz  "
          f"PM={pm_d:.1f} deg  (target ~10 deg)\n")
except Exception:
    print("Plant+delay: no finite crossover\n")

# ---------------------------------------------------------------------------
# 3. WEIGHTS
# ---------------------------------------------------------------------------
# W1 = two first-order weights multiplied:
#   W1_low: heavy suppression below 5 kHz
#   W1_mid: light anchor up to 25 kHz to prevent crossover dropping
# Product has no S+T=1 conflict because W1_MID=3 is small and W3_DC=0.05
W1_low = makeweight(W1_DC,  w(W1_CORNER),     W1_HF)
W1_mid = makeweight(W1_MID, w(W1_MID_CORNER), W1_HF)
W1     = W1_low * W1_mid

W2 = makeweight(W2_DC, w(W2_CORNER), W2_HF)
W3 = makeweight(W3_DC, w(W3_CORNER), W3_HF)

print(f"W1 demands: |S| < 1/{W1_DC:.0f} = {-20*np.log10(W1_DC):.0f} dB  below {W1_CORNER/1e3:.0f} kHz")
print(f"           |S| < 1/{W1_MID:.0f} = {-20*np.log10(W1_MID):.1f} dB  up to {W1_MID_CORNER/1e3:.0f} kHz  (xover anchor)")
print(f"W3 DC={W3_DC} — pure HF guardrail, no conflict with W1\n")

# ---------------------------------------------------------------------------
# 4. SYNTHESIS
# ---------------------------------------------------------------------------
print("Running mixsyn...")
t0 = time.time()
K, CL, info = ct.mixsyn(G_aug, w1=W1, w2=W2, w3=W3)
elapsed = time.time() - t0

gamma = float(np.atleast_1d(info[0]).flat[0])
rcond = float(np.atleast_1d(info[1]).flat[0])
print(f"  Done {elapsed:.1f}s  gamma={gamma:.4f}  rcond={rcond:.2e}  "
      f"order={ct.ss(K).nstates}")

if gamma < 0.8:
    print(f"  gamma={gamma:.2f}: lots of headroom — raise W1_DC for more suppression")
elif gamma < 1.5:
    print(f"  gamma={gamma:.2f}: near-optimal — plant is close to its limit here")
elif gamma < 3.0:
    print(f"  gamma={gamma:.2f}: acceptable — could push W1_DC a little harder")
elif gamma < 5.0:
    print(f"  gamma={gamma:.2f}: tight — lower W1_DC by ~30%")
else:
    print(f"  gamma={gamma:.2f} > 5: overconstrained — lower W1_DC significantly")

# ---------------------------------------------------------------------------
# 5. LOOP ANALYSIS
# ---------------------------------------------------------------------------
w_dense = np.logspace(-3, 3, 4000)
f_plot  = w_dense * w_norm / (2*np.pi)

def freq_mag(sys):
    return np.abs(np.asarray(
        ct.frequency_response(sys, w_dense).frdata).squeeze())

def freq_phase_deg(sys):
    return np.degrees(np.unwrap(np.angle(np.asarray(
        ct.frequency_response(sys, w_dense).frdata).squeeze())))

L_new = G_aug * K;  S_new = ct.feedback(1, L_new);  T_new = ct.feedback(L_new, 1)
L_old = G_aug;      S_old = ct.feedback(1, L_old);   T_old = ct.feedback(L_old, 1)

def margins_report(L, label):
    try:
        gm, pm, _, wpc = ct.margin(L)
        fc = float(wpc) * w_norm / (2*np.pi)
        if not np.isfinite(pm): raise ValueError
    except Exception:
        print(f"  [{label}]  no finite crossover"); return None, None
    Ms = np.max(freq_mag(ct.feedback(1, L)))
    print(f"  [{label:18s}]  xover={fc/1e3:5.1f} kHz  PM={pm:5.1f} deg  "
          f"GM={20*np.log10(max(gm,1e-9)):5.1f} dB  peak|S|={Ms:.2f}")
    return fc, pm

print()
fc_new, pm_new = margins_report(L_new, "H∞ outer K")
fc_old, pm_old = margins_report(L_old, "existing  ")

S_new_mag = freq_mag(S_new)
S_old_mag = freq_mag(S_old)

print()
print("  Suppression vs existing loop:")
any_worse = False
for f_check in [100, 500, 1000, 2000, 5000, 10000]:
    idx = np.argmin(np.abs(f_plot - f_check))
    s_n  = 20*np.log10(max(S_new_mag[idx], 1e-12))
    s_o  = 20*np.log10(max(S_old_mag[idx], 1e-12))
    imp  = s_o - s_n
    flag = "  *** WORSE" if imp < 0 else f"  ({10**(imp/20):.1f}x better)"
    if imp < 0: any_worse = True
    print(f"    {f_check:6d} Hz:  existing={s_o:+.1f} dB  "
          f"H∞={s_n:+.1f} dB  Δ={imp:+.1f} dB{flag}")

if any_worse:
    print("\n  *** Worse at some freqs — raise W1_MID or W1_DC")
else:
    print("\n  H∞ better than existing at all checked frequencies.")

# ---------------------------------------------------------------------------
# 6. PLOTS
# ---------------------------------------------------------------------------
cn, co = "#1f77b4", "#d62728"
fig, axs = plt.subplots(3, 2, figsize=(14, 12))
ttl = f"Intensity lock v12  |  γ={gamma:.3f}"
if pm_new: ttl += f"  |  H∞: PM={pm_new:.1f}°@{fc_new/1e3:.1f}kHz"
if pm_old: ttl += f"  |  existing: PM={pm_old:.1f}°@{fc_old/1e3:.1f}kHz"
fig.suptitle(ttl, fontsize=10)

# PRIMARY: Sensitivity
ax = axs[0, 0]
ax.loglog(f_plot, S_new_mag, color=cn, lw=2.0, label="|S| H∞")
ax.loglog(f_plot, S_old_mag, color=co, lw=1.5, ls="--", label="|S| existing")
ax.loglog(f_plot, 1/freq_mag(W1), "k--", lw=0.8, label="1/W1 target")
ax.axhline(1.0, color="gray", ls=":", lw=0.8)
ax.axvline(W1_CORNER,     color="green",  ls=":", lw=1.0, label=f"{W1_CORNER/1e3:.0f} kHz band edge")
ax.axvline(W1_MID_CORNER, color="orange", ls=":", lw=1.0, label=f"{W1_MID_CORNER/1e3:.0f} kHz xover anchor")
ax.set_ylim([1e-4, 10])
ax.set_title(f"★ Sensitivity  |  target: {20*np.log10(W1_DC):.0f} dB better below {W1_CORNER/1e3:.0f} kHz")
ax.set_xlabel("Hz"); ax.legend(fontsize=7); ax.grid(True, which="both", alpha=0.3)
for f_ann in [100, 1000, 5000]:
    idx = np.argmin(np.abs(f_plot - f_ann))
    val = S_new_mag[idx]
    if val > 1e-4:
        ax.annotate(f"{20*np.log10(max(val,1e-12)):.0f}dB",
                    xy=(f_ann, val), xytext=(f_ann, min(val*6, 1.5)),
                    fontsize=7, color=cn, ha='center',
                    arrowprops=dict(arrowstyle='->', color=cn, lw=0.8))

# Phase
ax = axs[0, 1]
ph_new = freq_phase_deg(L_new)
ph_old = freq_phase_deg(L_old)
ax.semilogx(f_plot, ph_new, color=cn, lw=1.8, label="∠L H∞")
ax.semilogx(f_plot, ph_old, color=co, lw=1.5, ls="--", label="∠L existing")
ax.axhline(-180, color="gray", ls=":", lw=0.8)
ax.axvline(f_xover, color="gray", ls=":", lw=0.8)
if fc_new and pm_new:
    idx = np.argmin(np.abs(f_plot - fc_new))
    ph_at = ph_new[idx]
    ax.annotate("", xy=(fc_new, -180), xytext=(fc_new, ph_at),
                arrowprops=dict(arrowstyle="<->", color=cn, lw=1.2))
    ax.text(fc_new*1.12, ph_at + pm_new*0.35, f"{pm_new:.0f}°", color=cn, fontsize=9)
if fc_old and pm_old:
    idx = np.argmin(np.abs(f_plot - fc_old))
    ph_at = ph_old[idx]
    ax.annotate("", xy=(fc_old*0.65, -180), xytext=(fc_old*0.65, ph_at),
                arrowprops=dict(arrowstyle="<->", color=co, lw=1.2))
    ax.text(fc_old*0.35, ph_at + pm_old*0.35, f"{pm_old:.0f}°", color=co, fontsize=9)
ax.set_ylim([-360, 0])
ax.set_title("Open-loop phase")
ax.set_xlabel("Hz"); ax.legend(fontsize=8); ax.grid(True, which="both", alpha=0.3)

# Open-loop magnitude
ax = axs[1, 0]
ax.loglog(f_plot, freq_mag(L_new), color=cn, lw=1.8, label="|L| H∞")
ax.loglog(f_plot, freq_mag(L_old), color=co, lw=1.5, ls="--", label="|L| existing")
ax.axvline(f_xover, color="gray", ls=":", lw=0.8)
ax.axhline(1.0,     color="gray", ls=":", lw=0.8)
ax.set_title("Open-loop magnitude")
ax.set_xlabel("Hz"); ax.legend(fontsize=8); ax.grid(True, which="both", alpha=0.3)

# Complementary sensitivity
ax = axs[1, 1]
ax.loglog(f_plot, freq_mag(T_new), color=cn, lw=1.8, label="|T| H∞")
ax.loglog(f_plot, freq_mag(T_old), color=co, lw=1.5, ls="--", label="|T| existing")
ax.loglog(f_plot, 1/freq_mag(W3),  "k--", lw=0.8, label="1/W3 bound")
ax.set_ylim([1e-4, 10])
ax.set_title("Complementary sensitivity")
ax.set_xlabel("Hz"); ax.legend(fontsize=8); ax.grid(True, which="both", alpha=0.3)

# Controller
ax = axs[2, 0]
ax.loglog(f_plot, freq_mag(K),     color=cn, lw=1.8, label="|K| H∞")
ax.loglog(f_plot, 1/freq_mag(W2),  "g--",    lw=0.8, label="1/W2 effort bound")
ax.axvline(W1_CORNER,     color="green",  ls=":", lw=1.0, label="5 kHz")
ax.axvline(W1_MID_CORNER, color="orange", ls=":", lw=1.0, label="25 kHz")
ax.axvline(f_xover,       color="gray",   ls=":", lw=0.8, label="30 kHz")
ax.axhline(10.0, color="red", ls=":", lw=0.8, label="|K|=10 danger")
ax.set_title("|K| shape  (high LF, rolloff above 30 kHz)")
ax.set_xlabel("Hz"); ax.legend(fontsize=7); ax.grid(True, which="both", alpha=0.3)

# Step response
ax = axs[2, 1]
t_n, y_n = ct.step_response(T_new)
t_o, y_o = ct.step_response(T_old)
ax.plot(1e6*t_n/w_norm, y_n, color=cn, lw=1.8, label="H∞")
ax.plot(1e6*t_o/w_norm, y_o, color=co, lw=1.5, ls="--", label="existing")
ax.set_xlabel("Time (µs)"); ax.set_title("Step response")
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("hinf_comparison.png", dpi=120)
print("\nSaved hinf_comparison.png")
plt.show()

# ---------------------------------------------------------------------------
# 7. DISCRETIZE  (31.25 MHz = 125 MHz / 4)
# ---------------------------------------------------------------------------
fs_ctrl = 125e6 / 4
Ts_n    = (1.0 / fs_ctrl) * w_norm
K_d     = ct.c2d(K, Ts_n, method="tustin")

try:
    K_bal = ct.balred(K_d, orders=K_d.nstates)
except Exception as e:
    print(f"  balred failed ({e}), using unbalanced realization")
    K_bal = ct.ss(K_d)

A, B, C, D = K_bal.A, K_bal.B, K_bal.C, K_bal.D
max_r = np.max(np.abs(np.linalg.eigvals(A)))
print(f"\n  Discrete poles: max|z|={max_r:.8f}  "
      f"({'OK' if max_r < 1.0 else 'UNSTABLE'})")

np.set_printoptions(precision=9, suppress=False, linewidth=120)
print(f"\n=== Discrete K @ {fs_ctrl/1e6:.4f} MHz ===")
print("A =\n", A)
print("B =\n", B)
print("C =\n", C)
print("D =\n", D)

np.savez("K_coeffs.npz", A=A, B=B, C=C, D=D, fs=np.array(fs_ctrl))
print("\nSaved K_coeffs.npz")

import scipy.signal as sig

# Extract continuous ZPK directly from the optimal controller K
num, den = ct.tfdata(ct.tf(K))
z_c, p_c, k_c = sig.tf2zpk(num[0][0], den[0][0])

fs_norm = fs_ctrl / w_norm

# Discretize via Bilinear Transform (mathematically identical to c2d tustin)
z_d, p_d, k_d = sig.bilinear_zpk(z_c, p_c, k_c, fs_norm)

# Save the exact discrete ZPK to a dedicated file for the HLS analyzer
np.savez("K_zpk.npz", z=z_d, p=p_d, k=k_d, fs=np.array(fs_ctrl))
print("\nSaved K_zpk.npz (Exact Discrete Poles and Zeros)")