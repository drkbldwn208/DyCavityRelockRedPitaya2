"""
H-infinity controller synthesis.

Default plant is a single-pole 50 kHz low-pass with a 500 ns electronics lag.
Target loop-gain crossover is 20 kHz. Swap in a fitted plant from bode_fit.py
via the --plant-npz flag when you have real measurement data.

  python3 scripts/ctrl.py                      # defaults
  python3 scripts/ctrl.py --xover 30           # retarget crossover
  python3 scripts/ctrl.py --plant-npz fit.npz  # use measured plant

Outputs:
  K_zpk.npz            — discrete controller ZPK (consumed by coeffs_analyze.py)
  hinf_comparison.png  — 3-panel Bode: Loop Gain, Phase, and Sensitivity
"""

import argparse
import time
import numpy as np
import scipy.signal as sig
import control as ct
import matplotlib.pyplot as plt

# ---------- CLI ----------
def parse_args():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--xover",        type=float, default=20e3, help="target crossover, Hz")
    p.add_argument("--plant-corner", type=float, default=50e3, help="plant LPF corner, Hz")
    p.add_argument("--plant-dc",     type=float, default=1.0,  help="plant DC gain")
    p.add_argument("--tau-delay",    type=float, default=0.5e-6, help="first-order delay, s")
    p.add_argument("--plant-npz",    type=str,   default=None, help="override plant from bode_fit ZPK npz")
    p.add_argument("--w1-dc",        type=float, default=30.0)
    p.add_argument("--w1-corner",    type=float, default=15e3)
    p.add_argument("--w2-dc",        type=float, default=1e-3)
    p.add_argument("--w2-corner",    type=float, default=80e3)
    p.add_argument("--w3-dc",        type=float, default=0.05)
    p.add_argument("--w3-corner",    type=float, default=40e3)
    p.add_argument("--fs",           type=float, default=125e6/128, help="controller sample rate, Hz")
    p.add_argument("--out",          type=str,   default="K_zpk.npz")
    p.add_argument("--no-show",      action="store_true")
    return p.parse_args()


# def makeweight(dcgain, wc_norm, hfgain):
#     """ Corrected H-infinity bounding weight formula """
#     # W(s) = (HF * s + DC * wc) / (s + wc)
#     return ct.tf([hfgain, wc_norm * dcgain], [1.0, wc_norm])

def makeweight(dcgain, wc_norm, hfgain):
    if dcgain > hfgain:
        M, A = dcgain, hfgain / dcgain
        return ct.tf([1.0/M, wc_norm], [1.0, wc_norm * A])
    M, A = hfgain, dcgain / hfgain
    return ct.tf([1.0, wc_norm * A], [1.0/M, wc_norm])




def build_plant(args, w_norm):
    if args.plant_npz:
        d = np.load(args.plant_npz)
        zc = d["z"] / w_norm
        pc = d["p"] / w_norm
        kc = float(d["k"]) * w_norm**(len(d["z"]) - len(d["p"]))
        num = np.poly(zc) * kc
        den = np.poly(pc)
        return ct.tf(num, den)
    wp_n = 2*np.pi*args.plant_corner / w_norm
    G    = ct.tf([args.plant_dc * wp_n], [1.0, wp_n])
    if args.tau_delay > 0:
        G *= ct.tf([1.0], [args.tau_delay * w_norm, 1.0])
    return G


def margins(L, label, w_norm):
    try:
        gm, pm, _, wpc = ct.margin(L)
        if not np.isfinite(pm):
            raise ValueError
        fc = float(wpc) * w_norm / (2*np.pi)
        print(f"  [{label:14s}] xover={fc/1e3:5.1f} kHz  PM={pm:5.1f} deg  "
              f"GM={20*np.log10(max(gm,1e-9)):5.1f} dB")
        return fc, pm
    except Exception:
        print(f"  [{label:14s}] no finite crossover")
        return None, None


def main():
    args    = parse_args()
    w_norm  = 2*np.pi*args.xover

    wn = lambda f: 2*np.pi*f / w_norm
    G_aug = build_plant(args, w_norm)

    W1 = makeweight(args.w1_dc, wn(args.w1_corner), 0.5)
    W2 = makeweight(args.w2_dc, wn(args.w2_corner), 1.0)
    W3 = makeweight(args.w3_dc, wn(args.w3_corner), 20.0)
    print(f"Target: xover ≈ {args.xover/1e3:.0f} kHz, |S|(DC) ≤ 1/{args.w1_dc:.0f}")

    t0 = time.time()
    K, _, info = ct.mixsyn(G_aug, w1=W1, w2=W2, w3=W3)
    gamma = float(np.atleast_1d(info[0]).flat[0])
    print(f"mixsyn: {time.time()-t0:.1f}s  γ={gamma:.3f}  order={ct.ss(K).nstates}")
    if gamma > 3.0:
        print("  γ > 3: over-asked — lower --w1-dc or widen corners")

    L_ctrl = G_aug * K
    S_ctrl = ct.feedback(1, L_ctrl) # Sensitivity S = 1 / (1 + L)
    
    print()
    fc_c, pm_c = margins(L_ctrl, "plant × K", w_norm)
    fc_o, pm_o = margins(G_aug,   "plant only", w_norm)
    
    peak_S = float(np.max(np.abs(np.asarray(
        ct.frequency_response(S_ctrl, np.logspace(-3, 3, 2000)).frdata).squeeze())))
    print(f"  peak |S| = {peak_S:.2f} ({20*np.log10(peak_S):.1f} dB)")

    # --- 3-Panel Bode Plot: Magnitude, Phase, and Sensitivity -----------------
    w_dense = np.logspace(-3, 3, 4000)
    f_plot  = w_dense * w_norm / (2*np.pi)
    mag = lambda s: np.abs(np.asarray(ct.frequency_response(s, w_dense).frdata).squeeze())
    phd = lambda s: np.degrees(np.unwrap(np.angle(np.asarray(
        ct.frequency_response(s, w_dense).frdata).squeeze())))

    # Make the figure taller to comfortably fit 3 panels
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 11), sharex=True)
    title = f"H-Infinity Synthesis  |  γ={gamma:.2f}"
    if pm_c: title += f"  |  with K: PM={pm_c:.0f}°@{fc_c/1e3:.1f} kHz"
    fig.suptitle(title, fontsize=14)
    
    # --- 1. Magnitude Plot ---
    ax1.loglog(f_plot, mag(L_ctrl), color="C0", lw=2.5, label="|L| (Loop Gain: G·K)")
    ax1.loglog(f_plot, mag(G_aug),  color="C3", lw=1.5, ls="--", label="|G| (Plant)")
    # Made controller orange, thicker, and dashed for high visibility
    ax1.loglog(f_plot, mag(K),      color="C1", lw=2.5, ls="--", label="|K| (Controller)") 
    ax1.axhline(1.0, color="gray", ls=":", lw=1.0)
    ax1.axvline(args.xover, color="purple", ls=":", lw=1.0, label=f"{args.xover/1e3:.0f} kHz target")
    ax1.set_ylabel("Magnitude", fontsize=11)
    ax1.grid(True, which="both", alpha=0.3)
    ax1.legend(fontsize=10, loc="lower left")

    # --- 2. Phase Plot ---
    ax2.semilogx(f_plot, phd(L_ctrl), color="C0", lw=2.5, label="∠L (Loop Gain)")
    ax2.semilogx(f_plot, phd(G_aug),  color="C3", lw=1.5, ls="--", label="∠G (Plant)")
    ax2.semilogx(f_plot, phd(K),      color="C1", lw=2.5, ls="--", label="∠K (Controller)")
    ax2.axhline(-180, color="gray", ls=":", lw=1.0)
    ax2.axvline(args.xover, color="purple", ls=":", lw=1.0)
    
    # Auto-scale phase y-axis to handle large integrator wrap-around gracefully
    current_ymin, current_ymax = ax2.get_ylim()
    ax2.set_ylim([min(-360, current_ymin), max(90, current_ymax)])
    
    ax2.set_ylabel("Phase (deg)", fontsize=11)
    ax2.grid(True, which="both", alpha=0.3)
    ax2.legend(fontsize=10, loc="lower left")

    # --- 3. Sensitivity Plot ---
    mag_S_db = 20 * np.log10(mag(S_ctrl) + 1e-12) # Convert Sensitivity to dB
    ax3.semilogx(f_plot, mag_S_db, color="C4", lw=2.5, label="|S| (Sensitivity)")
    
    # Prominent safety limit line
    ax3.axhline(6.0, color="red", ls="-.", lw=1.5, label="+6.0 dB Safety Limit")
    ax3.axvline(args.xover, color="purple", ls=":", lw=1.0)
    
    ax3.set_xlabel("Frequency (Hz)", fontsize=11)
    ax3.set_ylabel("Magnitude (dB)", fontsize=11)
    ax3.grid(True, which="both", alpha=0.3)
    ax3.legend(fontsize=10, loc="upper left")

    plt.tight_layout()
    plt.savefig("hinf_comparison.png", dpi=120)
    print("Saved hinf_comparison.png")

    # --- Discretize (Tustin) + save ZPK --------------------------------------
    Ts_n = (1.0 / args.fs) * w_norm
    K_d  = ct.c2d(K, Ts_n, method="tustin")
    A, B, C, D = np.array(K_d.A), np.array(K_d.B), np.array(K_d.C), np.array(K_d.D)
    n     = A.shape[0]
    p_d   = np.linalg.eigvals(A)
    z_d   = np.array(ct.zeros(K_d)).flatten()
    K_dc  = (D + C @ np.linalg.solve(np.eye(n) - A, B)).item().real
    k_d   = np.real(K_dc * np.prod(1.0 - p_d) / np.prod(1.0 - z_d))

    # cancel numerically-exact pole-zero pairs
    zs, ps = list(z_d), list(p_d)
    for p in list(ps):
        for z in list(zs):
            if abs(p - z) / (abs(p) + 1e-30) < 1e-6:
                ps.remove(p); zs.remove(z); break
    if len(ps) != len(p_d):
        k_d = np.real(K_dc * np.prod(1.0 - np.array(ps)) / np.prod(1.0 - np.array(zs)))
    z_d, p_d = np.array(zs), np.array(ps)

    print(f"Discrete K: {len(p_d)}p/{len(z_d)}z  |p|_max={max(abs(p_d)):.6f}  "
          f"DC gain={K_dc:.2f} ({20*np.log10(abs(K_dc)):.1f} dB)")
    np.savez(args.out, z=z_d, p=p_d, k=k_d, fs=np.array(args.fs))
    print(f"Saved {args.out}")

    # --- Plant SOS for the C++ testbench -------------------------------------
    G_aug_d = ct.c2d(G_aug, Ts_n, method="tustin")
    zp, pp  = np.array(ct.zeros(G_aug_d)), np.array(ct.poles(G_aug_d))
    kp      = float(G_aug_d.dcgain().real) * np.prod(1.0 - pp) / np.prod(1.0 - zp)
    sos_p   = sig.zpk2sos(zp, pp, np.real(kp), pairing='nearest')
    print(f"\nPlant SOS for tb_freq_response.cpp (PLANT_N_SEC = {sos_p.shape[0]}):")
    for s in sos_p:
        print(f"    {{{s[0]:.8e}, {s[1]:.8e}, {s[2]:.8e}, {s[4]:.8e}, {s[5]:.8e}}},")

    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()