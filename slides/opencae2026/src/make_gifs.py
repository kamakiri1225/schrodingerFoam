"""Animated GIFs for the opencae2026 slides, drawn from the SAME runs as the paper figures.

  tunneling : run/00_1b_tunneling_visible   (paper Fig. tunneling)
  harmonic  : run/00_2_harmonicOscillator_1D (paper Fig. harmonic)
  release   : run/00_3_release2D/VTK         (paper Fig. release; classical row = analytic, as in the paper)
  soliton   : run/04_darkSoliton_whiteNoise/VTK (paper Fig. soliton; density + phase)

Style follows the paper figures (blue |psi|^2 fill, red dashed potential, jet / twilight colormaps).
Run:  python3 slides/opencae2026/src/make_gifs.py   (from the repository root or anywhere)
Output: slides/opencae2026/img/anim_*.gif
"""
import io
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
RUN = os.path.join(REPO, "run")
OUT = os.path.join(HERE, "..", "img")

plt.rcParams.update({"font.family": ["Noto Sans CJK JP", "DejaVu Sans"], "font.size": 17,
                     "axes.titlesize": 20, "axes.labelsize": 18,
                     "xtick.labelsize": 15, "ytick.labelsize": 15})
BLUE, BLUE_FILL, RED = "#2a6fb0", "#9cc3e6", "#e41a1c"
DPI = 100


# ---------------------------------------------------------------- readers
def read_internal(path):
    import re
    data = open(path, "rb").read()
    m = re.search(rb"internalField\s+nonuniform\s+List<scalar>\s*(\d+)\s*\(", data)
    if m is None:
        return float(re.search(rb"internalField\s+uniform\s+([-+0-9.eE]+)", data).group(1))
    n = int(m.group(1))
    if re.search(rb"format\s+binary", data):
        return np.frombuffer(data[m.end():m.end() + 8 * n], dtype="<f8").copy()
    return np.array(data[m.end():].split(b")", 1)[0].split(), dtype=float)[:n]


def time_dirs(case):
    ts = []
    for n in os.listdir(os.path.join(RUN, case)):
        try:
            ts.append((float(n), n))
        except ValueError:
            pass
    return sorted(ts)


def psi_1d(case, name):
    d = os.path.join(RUN, case, name)
    re_, im_ = read_internal(os.path.join(d, "Psire")), read_internal(os.path.join(d, "Psiim"))
    return re_ * re_ + im_ * im_


def load_series(case):
    vdir = os.path.join(RUN, case, "VTK")
    files = json.load(open(os.path.join(vdir, f"{case}.vtm.series"), encoding="utf-8"))["files"]
    return vdir, {float(e["time"]): e["name"] for e in files}


def read_psi_2d(vdir, vtm):
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy
    r = vtk.vtkXMLUnstructuredGridReader()
    r.SetFileName(os.path.join(vdir, vtm.replace(".vtm", ""), "internal.vtu")); r.Update()
    g = r.GetOutput(); c = vtk.vtkCellCenters(); c.SetInputData(g); c.Update()
    pts = vtk_to_numpy(c.GetOutput().GetPoints().GetData())
    re_ = vtk_to_numpy(g.GetCellData().GetArray("Psire")); im_ = vtk_to_numpy(g.GetCellData().GetArray("Psiim"))
    xs = np.unique(np.round(pts[:, 0], 7)); ys = np.unique(np.round(pts[:, 1], 7))
    ix = np.searchsorted(xs, np.round(pts[:, 0], 7)); iy = np.searchsorted(ys, np.round(pts[:, 1], 7))
    psi = np.zeros((len(ys), len(xs)), complex); psi[iy, ix] = re_ + 1j * im_
    return xs, ys, psi


def count_vortices(psi):
    ph = np.angle(psi)
    wrap = lambda a: (a + np.pi) % (2 * np.pi) - np.pi
    d1 = wrap(np.roll(ph, -1, 1) - ph)
    d2 = wrap(np.roll(np.roll(ph, -1, 1), -1, 0) - np.roll(ph, -1, 1))
    d3 = wrap(np.roll(ph, -1, 0) - np.roll(np.roll(ph, -1, 1), -1, 0))
    d4 = wrap(ph - np.roll(ph, -1, 0))
    w = np.rint((d1 + d2 + d3 + d4) / (2 * np.pi)).astype(int)
    return int(np.count_nonzero(w[:-1, :-1]))


# ---------------------------------------------------------------- gif writer
def grab(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=DPI, facecolor="white")
    plt.close(fig)
    buf.seek(0)
    return Image.open(buf).convert("RGB")


def write_gif(frames, name, ms, colors=128):
    pal = frames[len(frames) // 2].quantize(colors=colors, method=Image.Quantize.MEDIANCUT)
    q = [f.quantize(palette=pal, dither=Image.Dither.NONE) for f in frames]
    path = os.path.join(OUT, name)
    durations = [ms] * (len(q) - 1) + [1500]          # hold the last frame
    q[0].save(path, save_all=True, append_images=q[1:], duration=durations, loop=0, optimize=True)
    print(f"wrote {path}  {len(q)} frames  {os.path.getsize(path)/1e6:.2f} MB")


# ---------------------------------------------------------------- 1-D cases
def anim_1d(case, xr, nx, vfun, xlim, ylim, vmax, title, name, step=1, ms=100, barrier=False):
    dx = (xr[1] - xr[0]) / nx
    x = xr[0] + dx * (np.arange(nx) + 0.5)
    frames = []
    for t, n in time_dirs(case)[::step]:
        if t == 0:      # 0/ holds the #codeStream initial condition, not data
            continue
        p = psi_1d(case, n)
        fig, ax = plt.subplots(figsize=(12, 5.0))
        ax.fill_between(x, p, color=BLUE_FILL, lw=0)
        ax.plot(x, p, color=BLUE, lw=2.4)
        ax2 = ax.twinx()
        ax2.plot(x, vfun(x), color=RED, ls="--", lw=2.2)
        ax2.set_ylim(0, vmax); ax2.set_ylabel("V(x)", color=RED); ax2.tick_params(axis="y", colors=RED)
        ax.set_xlim(*xlim); ax.set_ylim(*ylim)
        ax.set_xlabel("x"); ax.set_ylabel(r"$|\psi|^2$")
        ax.set_title(f"{title}    t = {t:4.1f}", loc="left")
        if barrier:
            R = p[x < -0.4].sum() * dx; T = p[x > 0.4].sum() * dx
            ax.text(0.98, 0.93, f"反射 {R:.3f}   透過 {T:.3f}", transform=ax.transAxes, ha="right", va="top",
                    fontsize=18, bbox=dict(fc="white", ec="#cccccc"))
        fig.tight_layout()
        frames.append(grab(fig))
    write_gif(frames, name, ms)


# ---------------------------------------------------------------- release (classical vs quantum)
def anim_release(step=2, ms=90):
    wx, wy, t_off, sv = 0.03125, 2.0, 2.0, 1.0           # same as the paper (make_2d.py)
    sx0, sy0 = 1 / np.sqrt(2 * wx), 1 / np.sqrt(2 * wy)
    vdir, series = load_series("00_3_release2D")
    frames = []
    for t in sorted(series)[::step]:
        xs, ys, psi = read_psi_2d(vdir, series[t])
        q = np.abs(psi) ** 2
        xg, yg = np.meshgrid(xs, ys)
        tau = max(t - t_off, 0.0)
        c = np.exp(-0.5 * ((xg / np.hypot(sx0, sv * tau)) ** 2 + (yg / np.hypot(sy0, sv * tau)) ** 2))
        fig, axes = plt.subplots(1, 2, figsize=(12, 6.0))
        for ax, rho, lab in zip(axes, (c, q), ("古典粒子（熱的）", "量子（schrodingerFoam）")):
            im = ax.imshow(rho / rho.max(), origin="lower", cmap="jet", vmin=0, vmax=1,
                           extent=[xs.min(), xs.max(), ys.min(), ys.max()], interpolation="bilinear")
            ax.set_xlim(-16, 16); ax.set_ylim(-16, 16); ax.set_aspect("equal")
            ax.set_xticks([-10, 0, 10]); ax.set_yticks([-10, 0, 10])
            ax.set_title(lab); ax.set_xlabel("x")
        axes[0].set_ylabel("y")
        state = "トラップ中" if t < t_off else "開放後"
        fig.suptitle(f"t = {t:4.1f}（{state}）   各コマ内で規格化した |ψ|²", fontsize=20, y=0.97)
        fig.tight_layout(rect=(0, 0, 0.9, 0.95))
        cax = fig.add_axes([0.91, 0.14, 0.018, 0.68]); fig.colorbar(im, cax=cax)
        frames.append(grab(fig))
    write_gif(frames, "anim_release.gif", ms)


# ---------------------------------------------------------------- dark soliton (density + phase)
def anim_soliton(ms=220):
    vdir, series = load_series("04_darkSoliton_whiteNoise")
    data = [(t, *read_psi_2d(vdir, series[t])) for t in sorted(series) if t <= 120]
    rmax = max((np.abs(p) ** 2).max() for _, _, _, p in data)
    frames = []
    for t, xs, ys, psi in data:
        ext = [xs.min(), xs.max(), ys.min(), ys.max()]
        fig, axes = plt.subplots(1, 2, figsize=(12, 6.3))
        im0 = axes[0].imshow(np.abs(psi) ** 2, origin="lower", cmap="jet", vmin=0, vmax=rmax, extent=ext,
                             interpolation="bilinear")
        im1 = axes[1].imshow(np.angle(psi), origin="lower", cmap="twilight_shifted", vmin=-np.pi, vmax=np.pi,
                             extent=ext, interpolation="nearest")
        axes[0].set_title(r"密度 $|\psi|^2$"); axes[1].set_title(r"位相 arg $\psi$")
        for ax in axes:
            ax.set_xticks([-8, 0, 8]); ax.set_yticks([-8, 0, 8]); ax.set_xlabel("x")
        axes[0].set_ylabel("y")
        fig.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.03)
        cb = fig.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.03, ticks=[-np.pi, 0, np.pi])
        cb.ax.set_yticklabels(["−π", "0", "π"])
        fig.suptitle(f"ダークソリトンの崩壊と量子渦の生成    t = {t:5.1f}", fontsize=20, y=0.98)
        fig.tight_layout(rect=(0, 0.02, 1, 0.95))
        frames.append(grab(fig))
    write_gif([f.resize((1000, round(f.height * 1000 / f.width)), Image.LANCZOS) for f in frames],
              "anim_soliton.gif", ms, colors=96)


if __name__ == "__main__":
    import sys
    if sys.argv[1:] == ["soliton"]:
        anim_soliton(); sys.exit()
    anim_1d("00_1b_tunneling_visible", (-40, 40), 400, lambda x: np.where(np.abs(x) < 0.4, 2.0, 0.0),
            (-40, 40), (0, 0.39), 2.3, "量子トンネル効果（V₀ = 2、E = 1.125）", "anim_tunneling.gif",
            ms=160, barrier=True)
    anim_1d("00_2_harmonicOscillator_1D", (-12, 12), 240, lambda x: 0.5 * x ** 2,
            (-9, 9), (0, 0.67), 13.4, "調和振動子のコヒーレント状態（ω = 1、x₀ = 4）", "anim_harmonic.gif",
            ms=80)
    anim_release()
    anim_soliton()
