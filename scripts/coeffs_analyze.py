"""
hinf_coeffs_analyze.py  v3

Key fix: ground truth is now state-space direct evaluation, not polyval on
the monolithic TF. polyval(num, z) loses precision catastrophically when z
is near a pole — that's what caused the noisy blue line in v2's plots.

Error breakdown now tells you three separate things:
  1. SOS decomposition error (should be ~machine epsilon — float64 limit)
  2. Pure quantization error (what the fixed-point implementation adds)
  3. Total error vs ground truth
"""

import numpy as np
import scipy.signal as sig
import matplotlib.pyplot as plt
import control as ct

TOTAL_BITS        = 32
INT_BITS_OVERRIDE = 3
CANCEL_TOL        = 0

data = np.load("K_zpk.npz")
z_d, p_d, k_d = data["z"], data["p"], data["k"]
fs = float(data["fs"])
print(f"Controller order: {len(p_d)}   fs: {fs/1e6:.4f} MHz\n")

# Create SOS directly from the exact discrete ZPK
sos = sig.zpk2sos(z_d, p_d, k_d, pairing='nearest')
# After zpk2sos, remove sections with poles near z = -1
# Replace the "drop sections with poles near z=-1" block with:
# Stabilize any near-unit-circle poles by contracting inward
# MAX_POLE_MAG = 0.998  # Target maximum pole radius

# for i, s in enumerate(sos):
#     poles = np.roots([1.0, s[4], s[5]])
#     max_mag = max(abs(p) for p in poles)
#     if max_mag > MAX_POLE_MAG:
#         contract = MAX_POLE_MAG / max_mag
#         # Contracting poles: if poles are p1,p2 then new polynomial is
#         # (z - p1*c)(z - p2*c) = z^2 - (p1+p2)*c*z + p1*p2*c^2
#         # = z^2 + a1*contract*z + a2*contract^2
#         new_a1 = s[4] * contract
#         new_a2 = s[5] * contract**2
        
#         # Preserve DC gain: adjust numerator to compensate
#         old_dc_den = 1.0 + s[4] + s[5]
#         new_dc_den = 1.0 + new_a1 + new_a2
#         num_scale = new_dc_den / old_dc_den if abs(old_dc_den) > 1e-12 else 1.0
        
#         sos[i] = [s[0]*num_scale, s[1]*num_scale, s[2]*num_scale,
#                   1.0, new_a1, new_a2]
        
#         new_poles = np.roots([1.0, new_a1, new_a2])
#         print(f"  Section {i}: pole |z|={max_mag:.6f} → {max(abs(p) for p in new_poles):.6f} "
#               f"(contracted by {(1-contract)*100:.2f}%)")

n_sec = sos.shape[0]

print(f"\n=== SOS: {n_sec} sections ===")
for i, s in enumerate(sos):
    b0,b1,b2,a0,a1,a2 = s
    poles_s = np.roots([1.0, a1, a2])
    fp = sorted([abs(np.angle(p_))*fs/(2*np.pi) for p_ in poles_s])
    coefs = [abs(c) for c in [b0,b1,b2,a1,a2] if abs(c) > 1e-15]
    dr = 20*np.log10(max(coefs)/min(coefs)) if len(coefs) > 1 else 0
    print(f"  Sec {i}: poles≈{fp[0]:.0f}/{fp[1]:.0f} Hz  DR={dr:.1f} dB")
print()

# Auto format
all_coefs = []
for s in sos:
    all_coefs.extend([abs(s[0]),abs(s[1]),abs(s[2]),abs(s[4]),abs(s[5])])
max_coef = max(c for c in all_coefs if c > 1e-15)
min_coef = min(c for c in all_coefs if c > 1e-15)

if INT_BITS_OVERRIDE is not None:
    INT_BITS = INT_BITS_OVERRIDE
else:
    INT_BITS = max(2, int(np.ceil(np.log2(max_coef + 1))) + 1)
FRAC_BITS = TOTAL_BITS - INT_BITS
COEF_SCALE = 2**FRAC_BITS
COEF_MAX =  2**(TOTAL_BITS-1) - 1
COEF_MIN = -2**(TOTAL_BITS-1)

print(f"=== Q{INT_BITS}.{FRAC_BITS} ===")
print(f"  Max |coef|={max_coef:.4f}   Min |coef|={min_coef:.4e}")
print(f"  Headroom: {2**(INT_BITS-1)/max_coef:.2f}x")
print(f"  Smallest coef: {min_coef*COEF_SCALE:.1f} LSBs ")
print()

def Q(x):
    q = int(round(x * COEF_SCALE))
    return max(COEF_MIN, min(COEF_MAX, q))

scaled_sos = [(Q(s[0]),Q(s[1]),Q(s[2]),Q(s[4]),Q(s[5])) for s in sos]

# Three frequency responses
w_hz  = np.logspace(1, np.log10(fs*0.49), 3000)
w_rad = w_hz / fs * 2*np.pi
z_e   = np.exp(1j * w_rad)

print("Computing factored ground truth...")
# Factored evaluation prevents Float64 polynomial precision loss near DC
H_truth = np.ones(len(z_e), dtype=complex) * k_d
for z in z_d:
    H_truth *= (z_e - z)
for p in p_d:
    H_truth /= (z_e - p)

H_sos_f = np.ones(len(w_hz), dtype=complex)
for s in sos:
    b = np.array([s[0],s[1],s[2]])
    a = np.array([1.0, s[4], s[5]])
    H_sos_f *= np.polyval(b, z_e) / np.polyval(a, z_e)

# data = np.load("K_coeffs.npz")
# A, B, C, D = data["A"], data["B"], data["C"], data["D"]
# fs = float(data["fs"])
# N  = A.shape[0]
# print(f"Controller order: {N}   fs: {fs/1e6:.4f} MHz\n")

# # zpk
# num, den = sig.ss2tf(A, B, C, D)
# num = num.flatten()
# z_k, p_k, k_k = sig.tf2zpk(num, den)


# # Cancel near pole/zero pairs
# z_list = list(z_k); p_list = list(p_k); cancelled = []
# for p in list(p_list):
#     for z in list(z_list):
#         if abs(p - z)/(abs(p)+1e-30) < CANCEL_TOL:
#             cancelled.append((p,z)); p_list.remove(p); z_list.remove(z); break
# if cancelled:
#     print(f"Cancelled {len(cancelled)} pole/zero pair(s) (Pade artifacts)")

# k_adj = k_k
# for p,z in cancelled:
#     d = abs(1.0 - p)
#     if d > 1e-10: k_adj *= abs(1.0 - z)/d

# # SOS
# if len(p_list) > 0:
#     sos = sig.zpk2sos(z_list, p_list, k_adj, pairing='nearest')
# else:
#     sos = np.array([[k_adj, 0, 0, 1, 0, 0]])
# n_sec = sos.shape[0]


# print(f"\n=== SOS: {n_sec} sections ===")
# for i, s in enumerate(sos):
#     b0,b1,b2,a0,a1,a2 = s
#     poles_s = np.roots([1.0, a1, a2])
#     fp = sorted([abs(np.angle(p_))*fs/(2*np.pi) for p_ in poles_s])
#     coefs = [abs(c) for c in [b0,b1,b2,a1,a2] if abs(c) > 1e-15]
#     dr = 20*np.log10(max(coefs)/min(coefs)) if len(coefs) > 1 else 0
#     print(f"  Sec {i}: poles≈{fp[0]:.0f}/{fp[1]:.0f} Hz  DR={dr:.1f} dB")
# print()

# # Auto format
# all_coefs = []
# for s in sos:
#     all_coefs.extend([abs(s[0]),abs(s[1]),abs(s[2]),abs(s[4]),abs(s[5])])
# max_coef = max(c for c in all_coefs if c > 1e-15)
# min_coef = min(c for c in all_coefs if c > 1e-15)

# if INT_BITS_OVERRIDE is not None:
#     INT_BITS = INT_BITS_OVERRIDE
# else:
#     INT_BITS = max(2, int(np.ceil(np.log2(max_coef + 1))) + 1)
# FRAC_BITS = TOTAL_BITS - INT_BITS
# COEF_SCALE = 2**FRAC_BITS
# COEF_MAX =  2**(TOTAL_BITS-1) - 1
# COEF_MIN = -2**(TOTAL_BITS-1)

# print(f"=== Q{INT_BITS}.{FRAC_BITS} ===")
# print(f"  Max |coef|={max_coef:.4f}   Min |coef|={min_coef:.4e}")
# print(f"  Headroom: {2**(INT_BITS-1)/max_coef:.2f}x")
# print(f"  Smallest coef: {min_coef*COEF_SCALE:.1f} LSBs "
#       f"({'adequate' if min_coef*COEF_SCALE > 10 else 'PRECISION LOSS'})")
# if max_coef < 2**(INT_BITS-1)/2 and INT_BITS > 2:
#     print(f"  HINT: could tighten to Q{INT_BITS-1}.{FRAC_BITS+1} "
#           f"(set INT_BITS_OVERRIDE={INT_BITS-1})")
# print()

# def Q(x):
#     q = int(round(x * COEF_SCALE))
#     return max(COEF_MIN, min(COEF_MAX, q))

# scaled_sos = [(Q(s[0]),Q(s[1]),Q(s[2]),Q(s[4]),Q(s[5])) for s in sos]

# # Three frequency responses
# w_hz  = np.logspace(1, np.log10(fs*0.49), 3000)
# w_rad = w_hz / fs * 2*np.pi
# z_e   = np.exp(1j * w_rad)

# print("Computing state-space ground truth (stable everywhere)...")
# nx = A.shape[0]
# I  = np.eye(nx)
# H_truth = np.zeros(len(z_e), dtype=complex)
# for j, zj in enumerate(z_e):
#     H_truth[j] = (C @ np.linalg.solve(zj*I - A, B) + D).item()

# H_sos_f = np.ones(len(w_hz), dtype=complex)
# for s in sos:
#     b = np.array([s[0],s[1],s[2]])
#     a = np.array([1.0, s[4], s[5]])
#     H_sos_f *= np.polyval(b, z_e) / np.polyval(a, z_e)

H_quant = np.ones(len(w_hz), dtype=complex)
H_secs = []
for b0q,b1q,b2q,a1q,a2q in scaled_sos:
    b = np.array([b0q,b1q,b2q]) / COEF_SCALE
    a = np.array([1.0, a1q/COEF_SCALE, a2q/COEF_SCALE])
    H_s = np.polyval(b, z_e) / np.polyval(a, z_e)
    H_quant *= H_s
    H_secs.append(H_s)

err_sos  = np.abs(H_sos_f - H_truth) / (np.abs(H_truth) + 1e-30)
err_quan = np.abs(H_quant - H_truth) / (np.abs(H_truth) + 1e-30)
err_q    = np.abs(H_quant - H_sos_f) / (np.abs(H_sos_f) + 1e-30)

mask = w_hz < fs*0.4
print(f"\nErrors (< 0.4*fs):")
print(f"  SOS decomp (float) vs truth: {20*np.log10(np.max(err_sos[mask])+1e-30):+.1f} dB")
print(f"  Quantized vs truth:          {20*np.log10(np.max(err_quan[mask])+1e-30):+.1f} dB")
print(f"  Quantized vs SOS (pure Q):   {20*np.log10(np.max(err_q[mask])+1e-30):+.1f} dB")
print()

# Plots
fig, axs = plt.subplots(2, 2, figsize=(13, 9))
fig.suptitle(f"SOS verification v3  |  {n_sec} sections  |  Q{INT_BITS}.{FRAC_BITS}",
             fontsize=10)

ax = axs[0,0]
ax.semilogx(w_hz, 20*np.log10(np.abs(H_truth)+1e-30),
            "#1f77b4", lw=2.5, label="Truth (SS direct)")
ax.semilogx(w_hz, 20*np.log10(np.abs(H_sos_f)+1e-30),
            "#2ca02c", lw=1.2, alpha=0.8, label="SOS float")
ax.semilogx(w_hz, 20*np.log10(np.abs(H_quant)+1e-30),
            "#d62728", lw=1.0, ls="--", label="SOS quantized")
ax.set_xlabel("Hz"); ax.set_ylabel("dB"); ax.set_title("Magnitude")
ax.legend(fontsize=8); ax.grid(True, which="both", alpha=0.3)

ax = axs[0,1]
ax.semilogx(w_hz, np.degrees(np.unwrap(np.angle(H_truth))),
            "#1f77b4", lw=2.5, label="Truth")
ax.semilogx(w_hz, np.degrees(np.unwrap(np.angle(H_sos_f))),
            "#2ca02c", lw=1.2, alpha=0.8, label="SOS float")
ax.semilogx(w_hz, np.degrees(np.unwrap(np.angle(H_quant))),
            "#d62728", lw=1.0, ls="--", label="SOS quantized")
ax.set_xlabel("Hz"); ax.set_ylabel("degrees"); ax.set_title("Phase")
ax.legend(fontsize=8); ax.grid(True, which="both", alpha=0.3)

ax = axs[1,0]
ax.semilogx(w_hz, 20*np.log10(err_sos+1e-30),  "#2ca02c", lw=1.2,
            label="SOS vs truth (should be ~-300 dB)")
ax.semilogx(w_hz, 20*np.log10(err_quan+1e-30), "#d62728", lw=1.5,
            label="Quantized vs truth")
ax.semilogx(w_hz, 20*np.log10(err_q+1e-30),    "#ff7f0e", lw=1.0, ls="--",
            label="Pure quantization (vs float SOS)")
ax.axhline(-40, color="gray", ls=":", lw=0.5)
ax.axhline(-60, color="gray", ls=":", lw=0.5)
ax.set_ylim([-160, 10])
ax.set_xlabel("Hz"); ax.set_ylabel("dB relative")
ax.set_title("Error breakdown  (orange = pure FP effect)")
ax.legend(fontsize=8); ax.grid(True, which="both", alpha=0.3)

ax = axs[1,1]
colors = ["#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd"]
for i, H_s in enumerate(H_secs):
    poles_s = np.roots([1.0, sos[i,4], sos[i,5]])
    fp = sorted([abs(np.angle(p_))*fs/(2*np.pi) for p_ in poles_s])
    ax.semilogx(w_hz, 20*np.log10(np.abs(H_s)+1e-30),
                color=colors[i%len(colors)], lw=1.5,
                label=f"Sec {i} ({fp[0]:.0f}/{fp[1]:.0f} Hz)")
ax.set_xlabel("Hz"); ax.set_ylabel("dB"); ax.set_title("Per-section responses")
ax.legend(fontsize=8); ax.grid(True, which="both", alpha=0.3)

plt.tight_layout()
plt.savefig("sos_verification.png", dpi=120)
print("Saved sos_verification.png")
plt.show()

# Write header
max_qerr_db = 20*np.log10(np.max(err_q[mask])+1e-30)
if max_qerr_db > -10:
    print(f"Pure quantization error {max_qerr_db:.1f} dB — header NOT written")
else:
    lines = [
        "// hinf_coeffs.h — auto-generated",
        f"// Q{INT_BITS}.{FRAC_BITS}, {n_sec} SOS sections, fs={fs/1e6:.4f} MHz",
        f"// Quantization error: {max_qerr_db:.1f} dB",
        "#ifndef HINF_COEFFS_H", "#define HINF_COEFFS_H",
        "#include <stdint.h>", "",
        f"#define HINF_N_SECTIONS     {n_sec}",
        f"#define HINF_COEF_INT_BITS  {INT_BITS}",
        f"#define HINF_COEF_FRAC_BITS {FRAC_BITS}",
        f"#define HINF_COEF_SCALE     {COEF_SCALE}L", "",
        "typedef struct { int32_t b0,b1,b2,a1,a2; } hinf_sos_t;", "",
        f"static const hinf_sos_t HINF_SOS[{n_sec}] = {{",
    ]
    for i, (b0q,b1q,b2q,a1q,a2q) in enumerate(scaled_sos):
        comma = "," if i < n_sec-1 else ""
        lines.append(f"    {{{b0q:11d},{b1q:11d},{b2q:11d},"
                     f"{a1q:11d},{a2q:11d}}}{comma}  // section {i}")
    lines += ["};", "", "#endif"]
    with open("src/hinf_coeffs.h","w") as f: f.write("\n".join(lines))
    print(f"Wrote hinf_coeffs.h  (Q err={max_qerr_db:.1f} dB)")