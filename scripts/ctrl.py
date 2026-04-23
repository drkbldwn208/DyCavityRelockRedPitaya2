"""
H-infinity controller synthesis for a simple first-order plant.

Plant (continuous):
    G(s) = K_DC / (1 + s / wp)        wp = 2π · PLANT_CORNER_HZ

This deliberately-simple 50 kHz low-pass is the verification plant. Once the
pipeline is happy (synth → discretize → quantize → HLS csim), swap in a
fitted plant (see bode_fit.py) to synthesize the real controller.

Target bandwidth (loop-gain crossover): 20 kHz.

Pipeline outputs:
    K_zpk.npz           — discrete controller ZPK  (consumed by coeffs_analyze.py)
    hinf_comparison.png — open-loop Bode, plant alone vs. plant × K

===========================================================================
How to choose W1, W2, W3 for stable performance
===========================================================================

Mixed-sensitivity design bounds three closed-loop transfers:

    |W1 · S|  ≤ 1      S = 1/(1+L)         (tracking / disturbance rejection)
    |W2 · K·S| ≤ 1     K·S                  (control effort)
    |W3 · T|  ≤ 1      T = L/(1+L)         (noise rejection / robustness)

S and T satisfy S + T = 1, so W1 and W3 CANNOT both be tight at the same
frequency. The weights trade loop shape across frequency:

  W1 (performance, low frequencies)
      DC gain:   large (e.g. 10–100) → forces |S| ≤ 1/W1_DC below W1_CORNER
      Corner:    the "bandwidth of good rejection" — put it a factor ~2
                 BELOW the desired crossover. For a 20 kHz bandwidth, use
                 W1_CORNER ≈ 8–10 kHz.
      HF gain:   small (< 1) so W1 does NOT fight W3 at high freq

  W3 (robustness, high frequencies)
      DC gain:   small (e.g. 0.01–0.1) — does nothing at low freq
      Corner:    put it a factor ~1.5 ABOVE the desired crossover so |T|
                 must roll off there. For 20 kHz, use W3_CORNER ≈ 30 kHz.
      HF gain:   large (> 5) → forces |T| to decay at HF

  W2 (effort)
      Usually very loose for a simple plant. Set |W2| small at low freq
      (allowing high K), and rising above actuator bandwidth so the
      controller can't command impossible things. |K| > 10 is often a
      warning the actuator model is missing.

Diagnostics:
  γ (gamma)     = attained H∞ norm. < 1.5 is near-optimal; > 3 means the
                  weights ask for something the plant cannot deliver
                  (lower W1_DC, move corners closer, or accept a lower BW).
  crossover     = where |L| = 1. Should land near W1_CORNER · √(W1_DC).
  phase margin  = should be > 30° (60° is comfortable).
  peak |S|      = < 2 (6 dB) means good robustness to modelling error.

Rule of thumb for picking crossover frequency:
  At the desired crossover, |G| · |K| must equal 1. A single-pole plant
  with corner f_p and DC gain K_DC has |G(f_c)| ≈ K_DC · f_p / f_c if
  f_c > f_p, or ≈ K_DC · (1 - (f_c/f_p)^2)^(1/2) if f_c < f_p. If the
  loop needs > ~20 dB gain at DC (W1_DC ≈ 10) and the plant delivers
  0 dB, the controller must supply the boost — so its |K(DC)| ≈ W1_DC.

This file's default weights produce a ~20 kHz crossover with γ ≈ 1–2 on
the 50 kHz plant; adjust the knobs below to retune.
"""

import time
import numpy as np
import scipy.signal as sig
import control as ct
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Normalization
# ---------------------------------------------------------------------------
F_XOVER_HZ = 20e3
W_NORM     = 2 * np.pi * F_XOVER_HZ

def wn(f_hz):
    return 2 * np.pi * f_hz / W_NORM

def makeweight(dcgain, wc_norm, hfgain):
    """First-order weight with specified DC gain, corner, and HF gain.
    Shape:  W(s) = (s/M + wc) / (s + wc·A)   where M,A depend on dcgain/hfgain.
    """
    if dcgain > hfgain:
        M, A = dcgain, hfgain / dcgain
        return ct.tf([1.0 / M, wc_norm], [1.0, wc_norm * A])
    M, A = hfgain, dcgain / hfgain
    return ct.tf([1.0, wc_norm * A], [1.0 / M, wc_norm])

# ===========================================================================
# Tuning knobs
# ===========================================================================
PLANT_CORNER_HZ = 50e3            # single-pole LPF corner
PLANT_DC_GAIN   = 1.0
TAU_DELAY_S     = 0.5e-6          # first-order electronics lag

# W1 — performance (low freq, |S| small)
W1_DC         = 30.0              # |S| < 1/30 = -30 dB at DC
W1_CORNER     = 15e3              # Hz — pushes xover toward 20 kHz
W1_HF         = 0.5               # keep <1 to avoid fighting W3

# W3 — robustness / noise (high freq, |T| small)
W3_DC         = 0.05              # harmless at DC
W3_CORNER     = 30e3              # Hz — ~1.5× F_XOVER_HZ
W3_HF         = 20.0              # aggressive HF rolloff request

# W2 — effort (small → more freedom)
W2_DC         = 1e-3
W2_CORNER     = 80e3
W2_HF         = 1.0

FS_FILTER     = 125e6 / 128       # 976.5625 kHz

# ---------------------------------------------------------------------------
# Plant
# ---------------------------------------------------------------------------
wp_n  = wn(PLANT_CORNER_HZ)
G     = ct.tf([PLANT_DC_GAIN * wp_n], [1.0, wp_n])
Delay = ct.tf([1.0], [TAU_DELAY_S * W_NORM, 1.0])
G_aug = G * Delay

try:
    _, pm0, _, wpc0 = ct.margin(G_aug)
    fc0 = float(wpc0) * W_NORM / (2 * np.pi)
    if np.isfinite(fc0) and np.isfinite(pm0):
        print(f"Plant+delay: xover={fc0/1e3:.1f} kHz  PM={pm0:.1f} deg")
    else:
        print("Plant+delay: |G| never crosses 0 dB  (no natural loop gain)")
except Exception:
    print("Plant+delay: no finite crossover")

# ---------------------------------------------------------------------------
# Weights
# ---------------------------------------------------------------------------
W1 = makeweight(W1_DC, wn(W1_CORNER), W1_HF)
W2 = makeweight(W2_DC, wn(W2_CORNER), W2_HF)
W3 = makeweight(W3_DC, wn(W3_CORNER), W3_HF)

print(f"W1 target: |S| < 1/{W1_DC:.0f} ({-20*np.log10(W1_DC):.0f} dB) below {W1_CORNER/1e3:.0f} kHz")
print(f"W3 target: |T| rolls off above {W3_CORNER/1e3:.0f} kHz")

# ---------------------------------------------------------------------------
# Synthesis
# ---------------------------------------------------------------------------
t0 = time.time()
K, CL, info = ct.mixsyn(G_aug, w1=W1, w2=W2, w3=W3)
gamma = float(np.atleast_1d(info[0]).flat[0])
rcond = float(np.atleast_1d(info[1]).flat[0])
print(f"mixsyn: {time.time()-t0:.1f}s  γ={gamma:.3f}  rcond={rcond:.2e}  "
      f"order={ct.ss(K).nstates}")

if gamma > 3.0:
    print("  γ > 3: weights are overconstrained — reduce W1_DC or move corners apart")
elif gamma < 0.5:
    print("  γ < 0.5: lots of headroom — raise W1_DC for tighter suppression")

# ---------------------------------------------------------------------------
# Loop analysis
# ---------------------------------------------------------------------------
w_dense = np.logspace(-3, 3, 4000)
f_plot  = w_dense * W_NORM / (2 * np.pi)

def freq_mag(sys):
    return np.abs(np.asarray(ct.frequency_response(sys, w_dense).frdata).squeeze())

def freq_phase_deg(sys):
    return np.degrees(np.unwrap(np.angle(np.asarray(
        ct.frequency_response(sys, w_dense).frdata).squeeze())))

L_open = G_aug
L_ctrl = G_aug * K

def margins(L, label):
    try:
        gm, pm, _, wpc = ct.margin(L)
        if not np.isfinite(pm):
            raise ValueError
        fc = float(wpc) * W_NORM / (2 * np.pi)
        print(f"  [{label:14s}] xover={fc/1e3:5.1f} kHz  PM={pm:5.1f} deg  "
              f"GM={20*np.log10(max(gm,1e-9)):5.1f} dB")
        return fc, pm
    except Exception:
        print(f"  [{label:14s}] no finite crossover")
        return None, None

print()
fc_c, pm_c = margins(L_ctrl, "plant × K")
fc_o, pm_o = margins(L_open, "plant only")
peak_S = np.max(freq_mag(ct.feedback(1, L_ctrl)))
print(f"  peak |S| = {peak_S:.2f} ({20*np.log10(peak_S):.1f} dB)")

# ---------------------------------------------------------------------------
# Plot: open-loop Bode, plant alone vs. plant × K
# ---------------------------------------------------------------------------
c_ctrl, c_open = "#1f77b4", "#d62728"
fig, axs = plt.subplots(2, 1, figsize=(10, 7), sharex=True)

ttl = f"Open-loop Bode  |  γ={gamma:.2f}"
if pm_c: ttl += f"  |  with K: PM={pm_c:.0f}°@{fc_c/1e3:.1f} kHz"
if pm_o: ttl += f"  |  plant only: PM={pm_o:.0f}°@{fc_o/1e3:.1f} kHz"
fig.suptitle(ttl)

ax = axs[0]
ax.loglog(f_plot, freq_mag(L_ctrl), color=c_ctrl, lw=2.0, label="|L| = |G·K|  (with controller)")
ax.loglog(f_plot, freq_mag(L_open), color=c_open, lw=1.5, ls="--", label="|G|  (plant only)")
ax.axhline(1.0, color="gray", ls=":", lw=0.8)
ax.axvline(PLANT_CORNER_HZ, color="green", ls=":", lw=0.8, label=f"{PLANT_CORNER_HZ/1e3:.0f} kHz plant corner")
ax.axvline(F_XOVER_HZ, color="purple", ls=":", lw=0.8, label=f"{F_XOVER_HZ/1e3:.0f} kHz xover target")
ax.set_ylabel("Magnitude")
ax.grid(True, which="both", alpha=0.3)
ax.legend(fontsize=9, loc="lower left")

ax = axs[1]
ax.semilogx(f_plot, freq_phase_deg(L_ctrl), color=c_ctrl, lw=2.0, label="∠L  (with controller)")
ax.semilogx(f_plot, freq_phase_deg(L_open), color=c_open, lw=1.5, ls="--", label="∠G  (plant only)")
ax.axhline(-180, color="gray", ls=":", lw=0.8)
if fc_c: ax.axvline(fc_c, color=c_ctrl, ls=":", lw=0.8)
if fc_o: ax.axvline(fc_o, color=c_open, ls=":", lw=0.8)
ax.set_ylabel("Phase (deg)")
ax.set_xlabel("Frequency (Hz)")
ax.set_ylim([-360, 90])
ax.grid(True, which="both", alpha=0.3)
ax.legend(fontsize=9, loc="lower left")

plt.tight_layout()
plt.savefig("hinf_comparison.png", dpi=120)
print("\nSaved hinf_comparison.png")

# ---------------------------------------------------------------------------
# Discretize controller (Tustin at fs = 125 MHz / 128)
# ---------------------------------------------------------------------------
Ts_n = (1.0 / FS_FILTER) * W_NORM
K_d  = ct.c2d(K, Ts_n, method="tustin")

A_d, B_d = np.array(K_d.A), np.array(K_d.B)
C_d, D_d = np.array(K_d.C), np.array(K_d.D)
n        = A_d.shape[0]

p_d  = np.linalg.eigvals(A_d)
z_d  = np.array(ct.zeros(K_d)).flatten()
K_dc = (D_d + C_d @ np.linalg.solve(np.eye(n) - A_d, B_d)).item().real
k_d  = np.real(K_dc * np.prod(1.0 - p_d) / np.prod(1.0 - z_d))

z_list, p_list, cancelled = list(z_d), list(p_d), 0
for p in list(p_list):
    for z in list(z_list):
        if abs(p - z) / (abs(p) + 1e-30) < 1e-6:
            p_list.remove(p); z_list.remove(z); cancelled += 1; break
if cancelled:
    k_d = np.real(K_dc * np.prod(1.0 - np.array(p_list))
                       / np.prod(1.0 - np.array(z_list)))
    print(f"Cancelled {cancelled} near-exact pole-zero pair(s)")

z_d, p_d = np.array(z_list), np.array(p_list)
print(f"Discrete K: {len(p_d)} poles, {len(z_d)} zeros, |p|_max={max(abs(p_d)):.6f}  "
      f"DC gain={K_dc:.3f} ({20*np.log10(abs(K_dc)):.1f} dB)")

np.savez("K_zpk.npz", z=z_d, p=p_d, k=k_d, fs=np.array(FS_FILTER))
print("Saved K_zpk.npz")

# ---------------------------------------------------------------------------
# Export plant SOS for C++ testbench closed-loop sim
# ---------------------------------------------------------------------------
G_aug_d = ct.c2d(G_aug, Ts_n, method="tustin")
zp, pp  = np.array(ct.zeros(G_aug_d)), np.array(ct.poles(G_aug_d))
kp      = float(G_aug_d.dcgain().real)
k_adj_p = kp * np.prod(1.0 - pp) / np.prod(1.0 - zp)
sos_p   = sig.zpk2sos(zp, pp, np.real(k_adj_p), pairing='nearest')

print("\n=== C++ Plant Simulator SOS (paste into tb_freq_response.cpp) ===")
print(f"const int PLANT_N_SEC = {sos_p.shape[0]};")
print("const double PLANT_SOS[][5] = {")
for s in sos_p:
    print(f"    {{{s[0]:.8e}, {s[1]:.8e}, {s[2]:.8e}, {s[4]:.8e}, {s[5]:.8e}}},")
print("};")

plt.show()
