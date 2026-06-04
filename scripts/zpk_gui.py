"""
Interactive plant/controller ZPK editor.

Loads a plant in the same continuous-time NPZ format used by scripts/ctrl.py:

    z: zeros in rad/s
    p: poles in rad/s
    k: continuous-time gain

The loaded plant is held fixed.  The editable roots are a hand-shaped controller
C(s).  The plot displays:

    G(s)       original loaded plant
    C(s)       editable controller
    G(s) C(s)  loop gain / plant times controller

Complex controller pairs are edited by natural frequency in Hz and Q, which is
usually the most physical way to move broad lead/lag or notch-like features.

Usage:
    python3 scripts/zpk_gui.py plant_fit.npz
    python3 scripts/zpk_gui.py plant_fit.npz --controller-npz controller_zpk.npz
    python3 scripts/zpk_gui.py
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import scipy.signal as sig

import tkinter as tk
from tkinter import filedialog, messagebox, ttk

import matplotlib

matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure


TWO_PI = 2.0 * math.pi


@dataclass
class RootGroup:
    role: str
    shape: str
    freq_hz: float
    q: float = 0.707
    real_hz: float | None = None

    def roots(self) -> list[complex]:
        freq_hz = max(float(self.freq_hz), 1e-12)
        if self.shape == "real":
            return [complex(-TWO_PI * freq_hz, 0.0)]

        wn = TWO_PI * freq_hz
        if self.real_hz is None:
            q = max(float(self.q), 0.500001)
            sigma = -wn / (2.0 * q)
        else:
            sigma = -abs(float(self.real_hz)) * TWO_PI
            if abs(sigma) >= wn:
                sigma = -0.999999 * wn
        wd = math.sqrt(max(0.0, wn * wn - sigma * sigma))
        return [complex(sigma, wd), complex(sigma, -wd)]

    def pair_real_hz(self) -> float:
        if self.shape == "real":
            return -abs(float(self.freq_hz))
        if self.real_hz is not None:
            return -abs(float(self.real_hz))
        return -float(self.freq_hz) / (2.0 * max(float(self.q), 0.500001))

    def pair_q(self) -> float:
        if self.shape == "real":
            return float("nan")
        real_hz = abs(self.pair_real_hz())
        return max(0.500001, float(self.freq_hz) / (2.0 * max(real_hz, 1e-300)))


def pair_to_freq_q(root: complex) -> tuple[float, float]:
    wn = abs(root)
    freq_hz = wn / TWO_PI
    if root.real >= 0:
        q = 1000.0
    else:
        q = max(0.500001, wn / (-2.0 * root.real))
    return freq_hz, q


def group_roots(roots: Iterable[complex], role: str) -> list[RootGroup]:
    roots = list(np.asarray(list(roots), dtype=complex).flatten())
    used = [False] * len(roots)
    groups: list[RootGroup] = []

    for i, root in enumerate(roots):
        if used[i]:
            continue
        used[i] = True

        if abs(root.imag) < max(1e-8, 1e-8 * abs(root)):
            groups.append(RootGroup(role=role, shape="real", freq_hz=abs(root.real) / TWO_PI))
            continue

        conj_idx = None
        for j in range(i + 1, len(roots)):
            if used[j]:
                continue
            if abs(roots[j] - np.conjugate(root)) <= 1e-7 * max(1.0, abs(root)):
                conj_idx = j
                break

        if conj_idx is None:
            groups.append(RootGroup(role=role, shape="pair", freq_hz=pair_to_freq_q(root)[0], q=pair_to_freq_q(root)[1]))
        else:
            used[conj_idx] = True
            freq_hz, q = pair_to_freq_q(root)
            groups.append(RootGroup(role=role, shape="pair", freq_hz=freq_hz, q=q))

    return groups


def flatten_groups(groups: list[RootGroup], role: str) -> np.ndarray:
    roots: list[complex] = []
    for group in groups:
        if group.role == role:
            roots.extend(group.roots())
    return np.asarray(roots, dtype=complex)


def dc_gain(z: np.ndarray, p: np.ndarray, k: float) -> float | None:
    if np.any(np.isclose(p, 0.0)):
        return None
    num = np.prod(-z) if len(z) else 1.0
    den = np.prod(-p) if len(p) else 1.0
    if abs(den) == 0:
        return None
    return float(np.real(k * num / den))


class ZpkEditor(tk.Tk):
    def __init__(self, npz_path: str | None):
        super().__init__()
        self.title("Plant + Controller ZPK Editor")
        self.geometry("1180x780")

        self.groups: list[RootGroup] = []
        self.k = 1.0  # controller gain
        self.plant_z = np.asarray([], dtype=complex)
        self.plant_p = np.asarray([], dtype=complex)
        self.plant_k = 1.0
        self.plant_extra_arrays: dict[str, np.ndarray] = {}
        self.controller_extra_arrays: dict[str, np.ndarray] = {}
        self.plant_path: Path | None = None
        self.controller_path: Path | None = None
        self.selected_iid: str | None = None

        self.f_min = tk.DoubleVar(value=1.0)
        self.f_max = tk.DoubleVar(value=200_000.0)
        self.n_points = tk.IntVar(value=2400)
        self.gain_var = tk.StringVar(value="1.0")
        self.status_var = tk.StringVar(value="Load a plant NPZ to begin.")

        self.role_var = tk.StringVar(value="pole")
        self.shape_var = tk.StringVar(value="pair")
        self.pair_param_var = tk.StringVar(value="Q")
        self.param_label_var = tk.StringVar(value="Q")
        self.freq_var = tk.StringVar(value="2500")
        self.q_var = tk.StringVar(value="10")
        self.freq_slider = tk.DoubleVar(value=math.log10(2500.0))
        self.q_slider = tk.DoubleVar(value=math.log10(10.0))
        self.export_mode = tk.StringVar(value="controller")

        self._build_ui()

        if npz_path:
            self.load_npz(Path(npz_path))
        else:
            self.update_plot()

    def _build_ui(self) -> None:
        self.columnconfigure(0, weight=0)
        self.columnconfigure(1, weight=1)
        self.rowconfigure(0, weight=1)

        side = ttk.Frame(self, padding=8)
        side.grid(row=0, column=0, sticky="nsw")
        side.columnconfigure(0, weight=1)

        buttons = ttk.Frame(side)
        buttons.grid(row=0, column=0, sticky="ew")
        ttk.Button(buttons, text="Load Plant", command=self.open_npz).grid(row=0, column=0, padx=(0, 4))
        ttk.Button(buttons, text="Load C", command=self.open_controller_npz).grid(row=0, column=1, padx=(0, 4))
        ttk.Button(buttons, text="Refresh", command=self.update_plot).grid(row=0, column=2)

        gain_box = ttk.LabelFrame(side, text="Controller Gain", padding=8)
        gain_box.grid(row=1, column=0, sticky="ew", pady=(8, 0))
        ttk.Label(gain_box, text="k").grid(row=0, column=0, sticky="w")
        ttk.Entry(gain_box, textvariable=self.gain_var, width=16).grid(row=0, column=1, sticky="ew", padx=(6, 0))
        ttk.Button(gain_box, text="Apply", command=self.apply_gain).grid(row=0, column=2, padx=(6, 0))
        ttk.Button(gain_box, text="Normalize C DC", command=self.normalize_dc_gain).grid(row=1, column=0, columnspan=3, sticky="ew", pady=(6, 0))

        tree_box = ttk.LabelFrame(side, text="Controller Roots", padding=8)
        tree_box.grid(row=2, column=0, sticky="nsew", pady=(8, 0))
        side.rowconfigure(2, weight=1)

        columns = ("role", "shape", "freq", "q", "real")
        self.tree = ttk.Treeview(tree_box, columns=columns, show="headings", height=15, selectmode="browse")
        for col, label, width in (
            ("role", "Type", 54),
            ("shape", "Shape", 54),
            ("freq", "Freq Hz", 86),
            ("q", "Q", 58),
            ("real", "Re Hz", 70),
        ):
            self.tree.heading(col, text=label)
            self.tree.column(col, width=width, anchor="e" if col in {"freq", "q"} else "center")
        self.tree.grid(row=0, column=0, sticky="nsew")
        scroll = ttk.Scrollbar(tree_box, orient="vertical", command=self.tree.yview)
        scroll.grid(row=0, column=1, sticky="ns")
        self.tree.configure(yscrollcommand=scroll.set)
        self.tree.bind("<<TreeviewSelect>>", self.on_select)
        tree_box.rowconfigure(0, weight=1)
        tree_box.columnconfigure(0, weight=1)

        edit = ttk.LabelFrame(side, text="Edit Selected", padding=8)
        edit.grid(row=3, column=0, sticky="ew", pady=(8, 0))
        ttk.Label(edit, text="Type").grid(row=0, column=0, sticky="w")
        ttk.Combobox(edit, textvariable=self.role_var, values=("pole", "zero"), width=8, state="readonly").grid(row=0, column=1, sticky="ew")
        ttk.Label(edit, text="Shape").grid(row=1, column=0, sticky="w")
        ttk.Combobox(edit, textvariable=self.shape_var, values=("real", "pair"), width=8, state="readonly").grid(row=1, column=1, sticky="ew")
        ttk.Label(edit, text="Pair field").grid(row=2, column=0, sticky="w")
        pair_mode = ttk.Combobox(edit, textvariable=self.pair_param_var, values=("Q", "Re Hz"), width=8, state="readonly")
        pair_mode.grid(row=2, column=1, sticky="ew")
        pair_mode.bind("<<ComboboxSelected>>", self.on_pair_param_mode)
        ttk.Label(edit, text="Freq Hz").grid(row=3, column=0, sticky="w")
        ttk.Entry(edit, textvariable=self.freq_var, width=12).grid(row=3, column=1, sticky="ew")
        ttk.Label(edit, textvariable=self.param_label_var).grid(row=4, column=0, sticky="w")
        ttk.Entry(edit, textvariable=self.q_var, width=12).grid(row=4, column=1, sticky="ew")
        ttk.Scale(edit, variable=self.freq_slider, from_=0.0, to=6.0, orient="horizontal", command=self.on_freq_slider).grid(row=5, column=0, columnspan=2, sticky="ew", pady=(6, 0))
        ttk.Scale(edit, variable=self.q_slider, from_=math.log10(0.51), to=math.log10(200.0), orient="horizontal", command=self.on_q_slider).grid(row=6, column=0, columnspan=2, sticky="ew", pady=(4, 0))
        ttk.Button(edit, text="Apply Root", command=self.apply_root).grid(row=7, column=0, columnspan=2, sticky="ew", pady=(6, 0))

        add = ttk.Frame(side)
        add.grid(row=4, column=0, sticky="ew", pady=(8, 0))
        ttk.Button(add, text="+ Pole Pair", command=lambda: self.add_root("pole", "pair")).grid(row=0, column=0, sticky="ew")
        ttk.Button(add, text="+ Zero Pair", command=lambda: self.add_root("zero", "pair")).grid(row=0, column=1, sticky="ew", padx=(4, 0))
        ttk.Button(add, text="+ Real Pole", command=lambda: self.add_root("pole", "real")).grid(row=1, column=0, sticky="ew", pady=(4, 0))
        ttk.Button(add, text="+ Real Zero", command=lambda: self.add_root("zero", "real")).grid(row=1, column=1, sticky="ew", padx=(4, 0), pady=(4, 0))
        ttk.Button(add, text="Delete Selected", command=self.delete_selected).grid(row=2, column=0, columnspan=2, sticky="ew", pady=(4, 0))

        export_box = ttk.LabelFrame(side, text="Export", padding=8)
        export_box.grid(row=5, column=0, sticky="ew", pady=(8, 0))
        ttk.Combobox(
            export_box,
            textvariable=self.export_mode,
            values=("controller", "plant_times_controller"),
            state="readonly",
            width=24,
        ).grid(row=0, column=0, sticky="ew")
        ttk.Button(export_box, text="Export NPZ", command=self.export_npz).grid(row=1, column=0, sticky="ew", pady=(6, 0))

        plot_opts = ttk.LabelFrame(side, text="Plot Range", padding=8)
        plot_opts.grid(row=6, column=0, sticky="ew", pady=(8, 0))
        ttk.Label(plot_opts, text="Min Hz").grid(row=0, column=0, sticky="w")
        ttk.Entry(plot_opts, textvariable=self.f_min, width=9).grid(row=0, column=1)
        ttk.Label(plot_opts, text="Max Hz").grid(row=1, column=0, sticky="w")
        ttk.Entry(plot_opts, textvariable=self.f_max, width=9).grid(row=1, column=1)
        ttk.Button(plot_opts, text="Apply Range", command=self.update_plot).grid(row=2, column=0, columnspan=2, sticky="ew", pady=(6, 0))

        ttk.Label(side, textvariable=self.status_var, wraplength=310).grid(row=7, column=0, sticky="ew", pady=(8, 0))

        plot_frame = ttk.Frame(self)
        plot_frame.grid(row=0, column=1, sticky="nsew")
        plot_frame.rowconfigure(0, weight=1)
        plot_frame.columnconfigure(0, weight=1)

        self.fig = Figure(figsize=(8.5, 6.5), dpi=100)
        self.ax_mag = self.fig.add_subplot(211)
        self.ax_phase = self.fig.add_subplot(212, sharex=self.ax_mag)
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")
        toolbar = NavigationToolbar2Tk(self.canvas, plot_frame, pack_toolbar=False)
        toolbar.update()
        toolbar.grid(row=1, column=0, sticky="ew")

    def open_npz(self) -> None:
        path = filedialog.askopenfilename(
            title="Load plant NPZ",
            filetypes=(("NPZ files", "*.npz"), ("All files", "*.*")),
        )
        if path:
            self.load_npz(Path(path))

    def open_controller_npz(self) -> None:
        path = filedialog.askopenfilename(
            title="Load controller NPZ",
            filetypes=(("NPZ files", "*.npz"), ("All files", "*.*")),
        )
        if path:
            self.load_controller_npz(Path(path))

    def load_npz(self, path: Path) -> None:
        try:
            data = np.load(path, allow_pickle=False)
            self.plant_z = np.asarray(data["z"], dtype=complex)
            self.plant_p = np.asarray(data["p"], dtype=complex)
            self.plant_k = float(np.asarray(data["k"]).item())
            self.plant_extra_arrays = {name: data[name] for name in data.files if name not in {"z", "p", "k"}}
        except Exception as exc:
            messagebox.showerror("Load failed", str(exc))
            return

        self.plant_path = path
        if not self.groups:
            self.k = 1.0
            self.gain_var.set("1.0")
        self.refresh_tree()
        self.status_var.set(
            f"Loaded plant {path} with {len(self.plant_z)} zeros, {len(self.plant_p)} poles. Editing C(s)."
        )
        self.update_plot()

    def load_controller_npz(self, path: Path) -> None:
        try:
            data = np.load(path, allow_pickle=False)
            z = np.asarray(data["z"], dtype=complex)
            p = np.asarray(data["p"], dtype=complex)
            self.k = float(np.asarray(data["k"]).item())
            self.controller_extra_arrays = {name: data[name] for name in data.files if name not in {"z", "p", "k"}}
        except Exception as exc:
            messagebox.showerror("Load failed", str(exc))
            return

        self.controller_path = path
        self.groups = group_roots(z, "zero") + group_roots(p, "pole")
        self.gain_var.set(f"{self.k:.12g}")
        self.refresh_tree()
        self.status_var.set(f"Loaded controller {path} with {len(z)} zeros, {len(p)} poles.")
        self.update_plot()

    def export_npz(self) -> None:
        z_c, p_c, k_c = self.current_zpk()
        mode = self.export_mode.get()
        if mode == "plant_times_controller":
            z = np.concatenate([self.plant_z, z_c])
            p = np.concatenate([self.plant_p, p_c])
            k = self.plant_k * k_c
            default = "plant_with_controller.npz"
        else:
            z, p, k = z_c, p_c, k_c
            default = "controller_zpk.npz"
            if self.controller_path is not None:
                default = f"{self.controller_path.stem}_edited.npz"

        path = filedialog.asksaveasfilename(
            title="Export ZPK",
            initialfile=default,
            defaultextension=".npz",
            filetypes=(("NPZ files", "*.npz"), ("All files", "*.*")),
        )
        if not path:
            return

        payload = dict(self.controller_extra_arrays if mode == "controller" else {})
        payload.update({"z": z, "p": p, "k": np.asarray(k)})
        try:
            np.savez(path, **payload)
        except Exception as exc:
            messagebox.showerror("Export failed", str(exc))
            return
        self.status_var.set(f"Exported {mode} {path} with {len(z)} zeros, {len(p)} poles.")

    def current_zpk(self) -> tuple[np.ndarray, np.ndarray, float]:
        self.apply_gain(show_errors=False)
        return flatten_groups(self.groups, "zero"), flatten_groups(self.groups, "pole"), float(self.k)

    def refresh_tree(self) -> None:
        self.tree.delete(*self.tree.get_children())
        for idx, group in enumerate(self.groups):
            q_text = "" if group.shape == "real" else f"{group.pair_q():.4g}"
            real_text = "" if group.shape == "real" else f"{group.pair_real_hz():.6g}"
            self.tree.insert(
                "",
                "end",
                iid=str(idx),
                values=(group.role, group.shape, f"{group.freq_hz:.6g}", q_text, real_text),
            )

    def on_select(self, _event=None) -> None:
        sel = self.tree.selection()
        if not sel:
            return
        self.selected_iid = sel[0]
        group = self.groups[int(sel[0])]
        self.role_var.set(group.role)
        self.shape_var.set(group.shape)
        self.pair_param_var.set("Re Hz" if group.real_hz is not None else "Q")
        self.freq_var.set(f"{group.freq_hz:.12g}")
        self.sync_pair_param_field(group)
        self.freq_slider.set(math.log10(max(group.freq_hz, 1.0)))

    def sync_pair_param_field(self, group: RootGroup) -> None:
        if self.pair_param_var.get() == "Re Hz":
            self.param_label_var.set("Re Hz")
            self.q_var.set(f"{group.pair_real_hz():.12g}")
            self.q_slider.configure(from_=math.log10(1e-3), to=math.log10(100_000.0))
            self.q_slider.set(math.log10(max(abs(group.pair_real_hz()), 1e-3)))
        else:
            self.param_label_var.set("Q")
            self.q_var.set(f"{group.pair_q():.12g}")
            self.q_slider.configure(from_=math.log10(0.51), to=math.log10(200.0))
            self.q_slider.set(math.log10(max(group.pair_q(), 0.51)))

    def on_pair_param_mode(self, _event=None) -> None:
        sel = self.tree.selection()
        if not sel:
            self.param_label_var.set(self.pair_param_var.get())
            return
        group = self.groups[int(sel[0])]
        self.sync_pair_param_field(group)

    def on_freq_slider(self, _value=None) -> None:
        freq = 10.0 ** float(self.freq_slider.get())
        self.freq_var.set(f"{freq:.7g}")
        self.apply_root(show_errors=False)

    def on_q_slider(self, _value=None) -> None:
        value = 10.0 ** float(self.q_slider.get())
        if self.pair_param_var.get() == "Re Hz":
            self.q_var.set(f"{-value:.7g}")
        else:
            self.q_var.set(f"{value:.7g}")
        self.apply_root(show_errors=False)

    def apply_gain(self, show_errors: bool = True) -> None:
        try:
            self.k = float(self.gain_var.get())
        except ValueError as exc:
            if show_errors:
                messagebox.showerror("Bad gain", str(exc))
        else:
            self.update_plot()

    def apply_root(self, show_errors: bool = True) -> None:
        if self.selected_iid is None:
            return
        try:
            idx = int(self.selected_iid)
            freq_hz = float(self.freq_var.get())
            pair_value = float(self.q_var.get())
            if freq_hz <= 0:
                raise ValueError("Frequency must be positive")
            if self.shape_var.get() == "pair":
                if self.pair_param_var.get() == "Q" and pair_value <= 0.5:
                    raise ValueError("Complex-pair Q must be > 0.5")
                if self.pair_param_var.get() == "Re Hz" and abs(pair_value) <= 0:
                    raise ValueError("Complex-pair Re Hz must be nonzero")
                if self.pair_param_var.get() == "Re Hz" and abs(pair_value) >= freq_hz:
                    raise ValueError("|Re Hz| must be smaller than Freq Hz for a complex pair")
        except Exception as exc:
            if show_errors:
                messagebox.showerror("Bad root", str(exc))
            return

        if self.shape_var.get() == "pair" and self.pair_param_var.get() == "Re Hz":
            real_hz = -abs(pair_value)
            q = max(0.500001, freq_hz / (2.0 * abs(real_hz)))
        else:
            real_hz = None
            q = max(pair_value, 0.500001)

        self.groups[idx] = RootGroup(
            role=self.role_var.get(),
            shape=self.shape_var.get(),
            freq_hz=freq_hz,
            q=q,
            real_hz=real_hz,
        )
        self.refresh_tree()
        self.tree.selection_set(str(idx))
        self.selected_iid = str(idx)
        self.update_plot()

    def add_root(self, role: str, shape: str) -> None:
        try:
            freq_hz = float(self.freq_var.get())
            pair_value = float(self.q_var.get())
        except ValueError:
            freq_hz = 2500.0
            pair_value = 10.0
        if shape == "pair" and self.pair_param_var.get() == "Re Hz":
            real_hz = -min(abs(pair_value), max(freq_hz * 0.99, 1e-12))
            q = max(0.500001, freq_hz / (2.0 * abs(real_hz)))
        else:
            real_hz = None
            q = max(pair_value, 0.500001)
        self.groups.append(
            RootGroup(role=role, shape=shape, freq_hz=max(freq_hz, 1e-12), q=q, real_hz=real_hz)
        )
        self.refresh_tree()
        iid = str(len(self.groups) - 1)
        self.tree.selection_set(iid)
        self.tree.see(iid)
        self.on_select()
        self.update_plot()

    def delete_selected(self) -> None:
        sel = self.tree.selection()
        if not sel:
            return
        del self.groups[int(sel[0])]
        self.selected_iid = None
        self.refresh_tree()
        self.update_plot()

    def normalize_dc_gain(self) -> None:
        z, p, k = self.current_zpk()
        current = dc_gain(z, p, k)
        if current is None or abs(current) < 1e-30:
            messagebox.showwarning("Cannot normalize", "DC gain is zero or infinite because of a root at/near the origin.")
            return
        self.k = k / current
        self.gain_var.set(f"{self.k:.12g}")
        self.status_var.set(f"Normalized controller DC gain from {current:.6g} to 1.")
        self.update_plot()

    def update_plot(self) -> None:
        z_c, p_c, k_c = self.current_zpk_no_side_effects()
        self.ax_mag.clear()
        self.ax_phase.clear()

        f_min = max(float(self.f_min.get()), 1e-6)
        f_max = max(float(self.f_max.get()), f_min * 1.01)
        n_points = max(int(self.n_points.get()), 100)
        f_hz = np.logspace(math.log10(f_min), math.log10(f_max), n_points)
        w = TWO_PI * f_hz

        _, G = sig.freqs_zpk(self.plant_z, self.plant_p, self.plant_k, worN=w)
        _, C = sig.freqs_zpk(z_c, p_c, k_c, worN=w)
        L = G * C

        traces = (
            ("G plant", G, "C0", "-"),
            ("C controller", C, "C1", "--"),
            ("G*C loop", L, "C3", "-"),
        )

        for label, H, color, style in traces:
            mag_db = 20.0 * np.log10(np.abs(H) + 1e-300)
            phase_deg = np.rad2deg(np.unwrap(np.angle(H)))
            lw = 2.4 if label == "G*C loop" else 1.8
            self.ax_mag.semilogx(f_hz, mag_db, lw=lw, color=color, ls=style, label=label)
            self.ax_phase.semilogx(f_hz, phase_deg, lw=lw, color=color, ls=style, label=label)

        for group in self.groups:
            color = "C1" if group.role == "zero" else "C2"
            ls = "--" if group.role == "zero" else ":"
            self.ax_mag.axvline(group.freq_hz, color=color, ls=ls, lw=0.8, alpha=0.45)
            self.ax_phase.axvline(group.freq_hz, color=color, ls=ls, lw=0.8, alpha=0.35)

        self.ax_mag.set_ylabel("Magnitude (dB)")
        self.ax_phase.set_ylabel("Phase (deg)")
        self.ax_phase.set_xlabel("Frequency (Hz)")
        self.ax_mag.grid(True, which="both", alpha=0.3)
        self.ax_phase.grid(True, which="both", alpha=0.3)
        self.ax_mag.legend(loc="best")
        self.ax_phase.legend(loc="best")

        dc_g = dc_gain(self.plant_z, self.plant_p, self.plant_k)
        dc_c = dc_gain(z_c, p_c, k_c)
        dc_l = None if dc_g is None or dc_c is None else dc_g * dc_c
        dc_text = "undefined" if dc_l is None else f"{dc_l:.5g} ({20 * math.log10(abs(dc_l) + 1e-300):+.2f} dB)"
        self.ax_mag.set_title(
            f"Plant + Controller: G {len(self.plant_z)}z/{len(self.plant_p)}p, "
            f"C {len(z_c)}z/{len(p_c)}p, loop DC={dc_text}"
        )
        self.fig.tight_layout()
        self.canvas.draw_idle()

    def current_zpk_no_side_effects(self) -> tuple[np.ndarray, np.ndarray, float]:
        try:
            k = float(self.gain_var.get())
        except ValueError:
            k = self.k
        return flatten_groups(self.groups, "zero"), flatten_groups(self.groups, "pole"), k


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("npz", nargs="?", help="plant NPZ to load")
    parser.add_argument("--controller-npz", help="optional continuous-time controller NPZ to load")
    args = parser.parse_args()
    app = ZpkEditor(args.npz)
    if args.controller_npz:
        app.load_controller_npz(Path(args.controller_npz))
    app.mainloop()


if __name__ == "__main__":
    main()
