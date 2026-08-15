#!/usr/bin/env python3
"""Render a 1D field (magSqrPsi) vs x from foamToVTK output into an animated GIF.

Usage:
    python3 render1d.py <case_VTK_dir> <out_dir> [field]

Plots |psi|^2 (filled) along x for every output time, and overlays the static
external potential Vext (dashed, read from the first frame) on a second y-axis.
Designed for the 1D tunneling / harmonic-oscillator cases.
"""
import sys, os, glob, json
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

vtk_dir = sys.argv[1]
out_dir = sys.argv[2]
field   = sys.argv[3] if len(sys.argv) > 3 else "magSqrPsi"
os.makedirs(out_dir, exist_ok=True)

series = glob.glob(os.path.join(vtk_dir, "*.vtm.series"))[0]
with open(series) as f:
    entries = json.load(f)["files"]


def read_internal(vtm_name, fname):
    idx_dir = vtm_name.replace(".vtm", "")
    vtu = os.path.join(vtk_dir, idx_dir, "internal.vtu")
    r = vtk.vtkXMLUnstructuredGridReader()
    r.SetFileName(vtu)
    r.Update()
    ug = r.GetOutput()
    cc = vtk.vtkCellCenters()
    cc.SetInputData(ug)
    cc.Update()
    pts = vtk_to_numpy(cc.GetOutput().GetPoints().GetData())
    arr = ug.GetCellData().GetArray(fname)
    vals = vtk_to_numpy(arr) if arr is not None else None
    return pts, vals


# x-axis ordering from the first frame
pts0, _ = read_internal(entries[0]["name"], field)
order = np.argsort(pts0[:, 0])
x = pts0[order, 0]

# static potential from the first frame (Vext is only written at t=0)
_, V = read_internal(entries[0]["name"], "Vext")
V = V[order] if V is not None else None

# global density scale
gmax = 0.0
data = []
for e in entries:
    _, vals = read_internal(e["name"], field)
    if vals is None:
        continue
    vals = vals[order]
    gmax = max(gmax, float(np.nanmax(vals)))
    data.append((e["time"], vals))

frames = []
for i, (t, vals) in enumerate(data):
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    ax.fill_between(x, vals, color="tab:blue", alpha=0.35)
    ax.plot(x, vals, color="tab:blue", lw=1.8, label=r"$|\psi|^2$")
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(0, 1.15*gmax)
    ax.set_xlabel("x")
    ax.set_ylabel(r"$|\psi|^2$", color="tab:blue")
    ax.tick_params(axis="y", labelcolor="tab:blue")
    ax.set_title(f"t = {t:.2f}", fontsize=12)

    if V is not None and np.nanmax(V) > 0:
        ax2 = ax.twinx()
        ax2.plot(x, V, color="tab:red", lw=1.4, ls="--", label=r"$V_{\rm ext}$")
        ax2.set_ylim(0, 1.15*float(np.nanmax(V)))
        ax2.set_ylabel(r"$V_{\rm ext}(x)$", color="tab:red")
        ax2.tick_params(axis="y", labelcolor="tab:red")

    fig.tight_layout()
    fpng = os.path.join(out_dir, f"{field}_{i:04d}.png")
    fig.savefig(fpng, dpi=90)
    plt.close(fig)
    frames.append(fpng)

imgs = [Image.open(f) for f in frames]
gif = os.path.join(out_dir, f"{field}.gif")
imgs[0].save(gif, save_all=True, append_images=imgs[1:], duration=100, loop=0)
print(f"wrote {len(frames)} frames + {gif}")
