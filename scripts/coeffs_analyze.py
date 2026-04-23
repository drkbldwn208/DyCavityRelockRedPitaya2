"""
Analyze and quantize the discrete controller from K_zpk.npz.

Pipeline:
  K_zpk.npz  →  SOS (scipy.zpk2sos)  →  Q{INT_BITS}.{FRAC_BITS} integers
             →  src/hinf_coeffs.h    (written if quant error < -10 dB)

Ground truth is factored ZPK evaluation (not polyval on monolithic TF, which
loses precision near poles). Reports three separate error budgets:
  1. SOS decomposition vs. ZPK truth  (should be ≈ float64 machine epsilon)
  2. Quantized SOS vs. ZPK truth       (total error hardware will see)
  3. Quantized SOS vs. float SOS       (pure quantization contribution)
"""

import numpy as np
import scipy.signal as sig
import matplotlib.pyplot as plt

TOTAL_BITS   = 32
INT_BITS     = 3
FRAC_BITS    = TOTAL_BITS - INT_BITS
COEF_SCALE   = 2 ** FRAC_BITS
COEF_MAX     =  2 ** (TOTAL_BITS - 1) - 1
COEF_MIN     = -2 ** (TOTAL_BITS - 1)

data = np.load("K_zpk.npz")
z_d, p_d, k_d = data["z"], data["p"], data["k"]
fs = float(data["fs"])
print(f"Controller: {len(p_d)} poles  fs={fs/1e6:.4f} MHz")

sos = sig.zpk2sos(z_d, p_d, k_d, pairing='nearest')
n_sec = sos.shape[0]

print(f"\nSOS: {n_sec} sections")
for i, s in enumerate(sos):
    poles_s = np.roots([1.0, s[4], s[5]])
    fp = sorted([abs(np.angle(q)) * fs / (2 * np.pi) for q in poles_s])
    coefs = [abs(c) for c in [s[0], s[1], s[2], s[4], s[5]] if abs(c) > 1e-15]
    dr = 20 * np.log10(max(coefs) / min(coefs)) if len(coefs) > 1 else 0
    print(f"  Sec {i}: pole freqs ≈ {fp[0]:.0f} / {fp[1]:.0f} Hz  DR={dr:.1f} dB")

max_coef = max(abs(c) for s in sos for c in (s[0], s[1], s[2], s[4], s[5]) if abs(c) > 1e-15)
headroom = 2 ** (INT_BITS - 1) / max_coef
print(f"\nQ{INT_BITS}.{FRAC_BITS}  (max |coef|={max_coef:.3f}, headroom={headroom:.2f}x)")
if headroom < 1.0:
    print(f"  WARNING: coefficient > integer range; increase INT_BITS")

def Q(x):
    return max(COEF_MIN, min(COEF_MAX, int(round(x * COEF_SCALE))))

scaled_sos = [(Q(s[0]), Q(s[1]), Q(s[2]), Q(s[4]), Q(s[5])) for s in sos]

# ---------------------------------------------------------------------------
# Frequency responses — factored ZPK truth, float SOS, quantized SOS
# ---------------------------------------------------------------------------
w_hz  = np.logspace(1, np.log10(fs * 0.49), 3000)
z_e   = np.exp(1j * w_hz / fs * 2 * np.pi)

H_truth = np.full(len(z_e), k_d, dtype=complex)
for z in z_d: H_truth *= (z_e - z)
for p in p_d: H_truth /= (z_e - p)

H_sos_f = np.ones(len(z_e), dtype=complex)
for s in sos:
    H_sos_f *= np.polyval([s[0], s[1], s[2]], z_e) / np.polyval([1.0, s[4], s[5]], z_e)

H_quant = np.ones(len(z_e), dtype=complex)
H_secs  = []
for b0q, b1q, b2q, a1q, a2q in scaled_sos:
    H_s = (np.polyval([b0q, b1q, b2q], z_e) / COEF_SCALE) / \
          np.polyval([1.0, a1q / COEF_SCALE, a2q / COEF_SCALE], z_e)
    H_quant *= H_s
    H_secs.append(H_s)

mask = w_hz < fs * 0.4
err_sos  = np.abs(H_sos_f - H_truth) / (np.abs(H_truth) + 1e-30)
err_quan = np.abs(H_quant - H_truth) / (np.abs(H_truth) + 1e-30)
err_q    = np.abs(H_quant - H_sos_f) / (np.abs(H_sos_f) + 1e-30)

print(f"\nErrors (f < 0.4·fs):")
print(f"  SOS decomp (float) vs truth: {20*np.log10(np.max(err_sos[mask])+1e-30):+.1f} dB")
print(f"  Quantized vs truth:          {20*np.log10(np.max(err_quan[mask])+1e-30):+.1f} dB")
print(f"  Quantized vs SOS (pure Q):   {20*np.log10(np.max(err_q[mask])+1e-30):+.1f} dB")

# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------
fig, axs = plt.subplots(2, 2, figsize=(13, 9))
fig.suptitle(f"SOS verification  |  {n_sec} sections  |  Q{INT_BITS}.{FRAC_BITS}")

ax = axs[0, 0]
ax.semilogx(w_hz, 20 * np.log10(np.abs(H_truth) + 1e-30), "#1f77b4", lw=2.5, label="Truth (ZPK)")
ax.semilogx(w_hz, 20 * np.log10(np.abs(H_sos_f) + 1e-30), "#2ca02c", lw=1.2, label="SOS float")
ax.semilogx(w_hz, 20 * np.log10(np.abs(H_quant) + 1e-30), "#d62728", lw=1.0, ls="--", label="SOS quantized")
ax.set_title("Magnitude"); ax.set_xlabel("Hz"); ax.set_ylabel("dB")
ax.grid(True, which="both", alpha=0.3); ax.legend(fontsize=8)

ax = axs[0, 1]
ax.semilogx(w_hz, np.degrees(np.unwrap(np.angle(H_truth))), "#1f77b4", lw=2.5, label="Truth")
ax.semilogx(w_hz, np.degrees(np.unwrap(np.angle(H_sos_f))), "#2ca02c", lw=1.2, label="SOS float")
ax.semilogx(w_hz, np.degrees(np.unwrap(np.angle(H_quant))), "#d62728", lw=1.0, ls="--", label="SOS quantized")
ax.set_title("Phase"); ax.set_xlabel("Hz"); ax.set_ylabel("degrees")
ax.grid(True, which="both", alpha=0.3); ax.legend(fontsize=8)

ax = axs[1, 0]
ax.semilogx(w_hz, 20 * np.log10(err_sos + 1e-30),  "#2ca02c", lw=1.2, label="SOS vs truth")
ax.semilogx(w_hz, 20 * np.log10(err_quan + 1e-30), "#d62728", lw=1.5, label="Quantized vs truth")
ax.semilogx(w_hz, 20 * np.log10(err_q + 1e-30),    "#ff7f0e", lw=1.0, ls="--", label="Pure quantization")
ax.set_title("Error breakdown"); ax.set_xlabel("Hz"); ax.set_ylabel("dB")
ax.set_ylim([-160, 10])
ax.grid(True, which="both", alpha=0.3); ax.legend(fontsize=8)

ax = axs[1, 1]
colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2"]
for i, H_s in enumerate(H_secs):
    poles_s = np.roots([1.0, sos[i, 4], sos[i, 5]])
    fp = sorted([abs(np.angle(q)) * fs / (2 * np.pi) for q in poles_s])
    ax.semilogx(w_hz, 20 * np.log10(np.abs(H_s) + 1e-30),
                color=colors[i % len(colors)], lw=1.5,
                label=f"Sec {i} ({fp[0]:.0f}/{fp[1]:.0f} Hz)")
ax.set_title("Per-section responses"); ax.set_xlabel("Hz"); ax.set_ylabel("dB")
ax.grid(True, which="both", alpha=0.3); ax.legend(fontsize=8)

plt.tight_layout()
plt.savefig("sos_verification.png", dpi=120)
print("Saved sos_verification.png")

# ---------------------------------------------------------------------------
# Write src/hinf_coeffs.h (only if quantization error is small)
# ---------------------------------------------------------------------------
max_qerr_db = 20 * np.log10(np.max(err_q[mask]) + 1e-30)
if max_qerr_db > -10:
    print(f"\nPure quantization error {max_qerr_db:.1f} dB — hinf_coeffs.h NOT written")
else:
    lines = [
        "// hinf_coeffs.h — auto-generated by coeffs_analyze.py",
        f"// Q{INT_BITS}.{FRAC_BITS}, {n_sec} SOS sections, fs={fs/1e6:.4f} MHz",
        f"// Quantization error: {max_qerr_db:.1f} dB",
        "#ifndef HINF_COEFFS_H",
        "#define HINF_COEFFS_H",
        "#include <stdint.h>",
        "",
        f"#define HINF_N_SECTIONS     {n_sec}",
        f"#define HINF_COEF_INT_BITS  {INT_BITS}",
        f"#define HINF_COEF_FRAC_BITS {FRAC_BITS}",
        f"#define HINF_COEF_SCALE     {COEF_SCALE}L",
        "",
        "typedef struct { int32_t b0,b1,b2,a1,a2; } hinf_sos_t;",
        "",
        f"static const hinf_sos_t HINF_SOS[{n_sec}] = {{",
    ]
    for i, (b0q, b1q, b2q, a1q, a2q) in enumerate(scaled_sos):
        comma = "," if i < n_sec - 1 else ""
        lines.append(f"    {{{b0q:11d},{b1q:11d},{b2q:11d},{a1q:11d},{a2q:11d}}}{comma}  // section {i}")
    lines += ["};", "", "#endif", ""]
    with open("src/hinf_coeffs.h", "w") as f:
        f.write("\n".join(lines))
    print(f"Wrote src/hinf_coeffs.h  (quant err = {max_qerr_db:.1f} dB)")

plt.show()
