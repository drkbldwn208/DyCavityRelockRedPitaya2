# DyCavityRelockRedPitaya2 — H∞ cavity-relock controller

Red Pitaya (Zynq xc7z010) servo that runs an H∞ controller on a decimated ADC
stream and drives the DAC. This README covers the **one workflow that matters
now**: take an already-fitted plant (`plant_only.npz`), synthesize a controller,
and build the bitstream that implements it.

> Plant *fitting* has moved to another workspace. The fitting/analysis scripts
> and stray `.npz`/`.png` files left here are bloat — see
> [Legacy / bloat](#legacy--bloat). You start from a finished `plant_only.npz`.

---

## TL;DR

```bash
# 1. Synthesize the H∞ controller from the fitted plant
python3 scripts/ctrl.py --plant-npz plant_only.npz --no-show
#    -> writes K_zpk.npz, hinf_comparison.png
#    -> PRINTS a "Plant SOS for tb_freq_response.cpp" block  <-- copy this

# 2. Quantize the controller to Q3.29 SOS (this is what the FPGA runs)
python3 scripts/coeffs_analyze.py --no-show
#    -> writes src/hinf_coeffs.h, sos_verification.png

# 3. *** MANUAL, EASY TO FORGET ***  paste the plant SOS from step 1 into
#    tb/tb_freq_response.cpp (PLANT_N_SEC + PLANT_SOS). Sim-only, but the
#    sim is meaningless without it. (See "The step you were missing".)

# 4. (Optional but recommended) verify in C-sim before the ~hour Vivado run
make csim_plot
python3 scripts/plot_closed_loop.py

# 5. Build the bitstream
make all
#    -> bitstream + .hwh land in pynq_overlays/latest/ via the bit monitor
```

---

## The step you were missing

There are **two independent consumers of the synthesis output**, and only one of
them updates itself:

| Output of `ctrl.py`                         | Consumer                          | Updated by |
|---------------------------------------------|-----------------------------------|------------|
| `K_zpk.npz` (controller ZPK)                | `coeffs_analyze.py` → `src/hinf_coeffs.h` → **bitstream** | automatic |
| printed *"Plant SOS for tb_freq_response.cpp"* | `tb/tb_freq_response.cpp` (sim plant model) | **you, by hand** |

`tb/tb_freq_response.cpp` ships with a hardcoded **single-pole default plant**:

```cpp
const int PLANT_N_SEC = 1;
const double PLANT_SOS[][5] = {
    {7.01024743e-02, 1.40204949e-01, 7.01024743e-02, -7.11018449e-01, -8.57165403e-03},
};
```

If you synthesize against `plant_only.npz` but never replace this block, the
C-sim runs your new fitted-plant controller against the **old 50 kHz LPF** —
which is exactly why the measured transfer function "looks like one pole" and
the closed-loop sim looks strange. `ctrl.py` prints the correct replacement; copy
the whole block (and update `PLANT_N_SEC`) into the testbench.

**Important:** this only affects the **C-simulation**. The bitstream from
`make all` is built from `src/` (i.e. `src/hinf_coeffs.h`) and never sees the
testbench, so the deployed controller is correct as long as steps 1–2 ran. The
stale plant only corrupts what the sim *shows you*.

---

## Workflow in detail

### 1. H∞ synthesis — `scripts/ctrl.py`
```bash
python3 scripts/ctrl.py --plant-npz plant_only.npz --no-show
```
- Loads the fitted plant ZPK, runs `control.mixsyn` with the W1/W2/W3 loop-shaping
  weights, discretizes the controller via **Tustin** at `fs = 125e6/128 =
  976.5625 kHz`.
- Writes **`K_zpk.npz`** (discrete controller z/p/k + fs) and
  **`hinf_comparison.png`** (loop gain / phase / sensitivity).
- Prints discrete controller poles/zeros and the **plant SOS block for the
  testbench**.
- Useful knobs: `--xover` (crossover, Hz), `--w1-dc` (DC sensitivity bound),
  `--out`. Watch the printed `γ`; `γ > 3` means you over-asked.

### 2. Quantize coefficients — `scripts/coeffs_analyze.py`
```bash
python3 scripts/coeffs_analyze.py --no-show          # reads K_zpk.npz
```
- Converts the controller to **Q3.29 second-order sections** and writes
  **`src/hinf_coeffs.h`** (the header compiled into the FPGA IP) plus
  **`sos_verification.png`**.
- Refuses to write the header if pure-quantization error exceeds
  `--max-qerr-db` (default −10 dB). Check the printed "Pure quantization" line.

### 3. Update the testbench plant (manual)
Paste the `ctrl.py` plant-SOS block into `tb/tb_freq_response.cpp`, setting both
`PLANT_N_SEC` and the `PLANT_SOS[][5]` rows. (No-op for the bitstream; required
for a meaningful sim.)

### 4. Verify in C-simulation (optional, recommended)
```bash
make csim_plot                          # csim + step/impulse/freq plots
python3 scripts/plot_closed_loop.py     # float vs fixed-point acquisition
python3 scripts/controller_fixedpoint_sim.py   # bit-accurate biquad model
```
The testbench drives `HinfFilter` directly at the filter rate. It does **not**
exercise the streaming top-level (`while(true)` + DATAFLOW can't run in csim);
use `cosim` or hardware for the full pipeline.

### 5. Build the bitstream — `make all`
`make all` = `external_cores` → `hls` → `vivado`:
- **`external_cores`** — `git submodule update --init --recursive` (pulls
  `extern/` cores).
- **`hls`** (`scripts/hls_build.tcl`) — `csynth_design` + `export_design` of
  `src/*` into `build/ip_repo/` (consumes `src/hinf_coeffs.h`).
- **`vivado`** (`scripts/vivado_build.tcl`) — block design, synth, impl,
  bitstream. First starts the **bit monitor** (`scripts/move_bit.sh`), which
  watches `build/` and copies each new `.bit` + matching `.hwh` into
  `pynq_overlays/latest/` (and a timestamped archive). Stop it with
  `make stop_monitor`.

Deploy the contents of `pynq_overlays/latest/` to the Red Pitaya / PYNQ.

---

## Prerequisites
- **Vitis HLS + Vivado 2024.1.** The Makefile sources them from
  `/home/levlabcukomen/tools/...`; edit those `source` lines for your machine.
- **Python 3** with `numpy`, `scipy`, `control`, `matplotlib`.
- `git` submodules for `extern/` and `inotifywait` (inotify-tools) for the bit
  monitor.

## Other make targets
| Target | Does |
|--------|------|
| `make csim` | C-simulation only (CSVs under `build/sim_workspace/.../csim/build`) |
| `make csim_plot` | csim + `plot_response.py` |
| `make hls` / `make vivado` | individual build stages |
| `make clean` | remove `build/`, logs |
| `make start_monitor` / `make stop_monitor` | manage the bitstream monitor |

## Key files
```
scripts/ctrl.py             H∞ synthesis        -> K_zpk.npz (+ plant SOS print)
scripts/coeffs_analyze.py   quantize controller -> src/hinf_coeffs.h
src/hinf_coeffs.h           Q3.29 SOS (auto-generated; the FPGA controller)
src/hinf_filter.hpp         fixed-point DF-I biquad cascade
src/dy_cavity_relocker_2.*  HLS top: decimate -> controller -> servo/DAC
tb/tb_freq_response.cpp     C-sim testbench (hand-pasted plant SOS lives here)
plant_only.npz              fitted plant ZPK (synthesis input)
Makefile                    build orchestration
```

## Legacy / bloat
These supported plant fitting, now done elsewhere — not part of the build:
- `scripts/bode_fit.py`, `scripts/deembed_integrator.py`,
  `scripts/extract_openloop.py`, `scripts/plot_bode_data.py`,
  `scripts/create_plant.py`
- stray data/figures: `plant_fit.npz`, `custom_plant.npz`, `stitched_plant.npz`,
  `*.csv` Bode exports, `fit_csv_exports/`, and the assorted `*.png`.

`plant_only.npz` is the only plant artifact the current workflow needs.
