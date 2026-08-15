#!/usr/bin/env python3
"""Render density (magSqrPsi) / phase from foamToVTK output into PNG frames + GIF.

Usage:
    python3 render.py <case_VTK_dir> <out_dir> [field] [vmax|auto]

Reads the *.vtm.series (time mapping) and each <block>_<idx>/internal.vtu,
rasterises the uniform structured grid and colour-maps it.
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
vmax_arg = sys.argv[4] if len(sys.argv) > 4 else None
vmax_override = (
    "auto" if vmax_arg == "auto"
    else float(vmax_arg) if vmax_arg is not None
    else None
)
os.makedirs(out_dir, exist_ok=True)

series = glob.glob(os.path.join(vtk_dir, "*.vtm.series"))[0]
with open(series) as f:
    entries = json.load(f)["files"]


def read_internal(vtm_name):
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
    arr = ug.GetCellData().GetArray(field)
    vals = vtk_to_numpy(arr) if arr is not None else None
    return pts, vals


# grid from first entry
pts0, _ = read_internal(entries[0]["name"])
xs = np.unique(np.round(pts0[:, 0], 6))
ys = np.unique(np.round(pts0[:, 1], 6))
nx, ny = len(xs), len(ys)
ix = np.searchsorted(xs, np.round(pts0[:, 0], 6))
iy = np.searchsorted(ys, np.round(pts0[:, 1], 6))

if field == "magSqrPsi":
    vmin, vmax, cmap, clabel = 0.0, 1.3, "jet", r"$|\psi|^2$"
else:
    vmin, vmax, cmap, clabel = -np.pi, np.pi, "twilight", r"$\arg\psi$"

if vmax_override is not None and vmax_override != "auto":
    vmax = vmax_override

frames = []
for frame_index, e in enumerate(entries):
    t = e["time"]
    pts, vals = read_internal(e["name"])
    if vals is None:
        continue
    grid = np.full((ny, nx), np.nan)
    grid[iy, ix] = vals
    frame_vmax = np.nanmax(vals) if vmax_override == "auto" else vmax
    fig, ax = plt.subplots(figsize=(6, 5.4))
    im = ax.imshow(grid, origin="lower",
                   extent=[xs.min(), xs.max(), ys.min(), ys.max()],
                   cmap=cmap, vmin=vmin, vmax=frame_vmax, interpolation="bilinear")
    auto_note = f"   max = {frame_vmax:.3g}" if vmax_override == "auto" else ""
    ax.set_title(f"{field}   t = {t}{auto_note}", fontsize=12)
    ax.set_xlabel("x"); ax.set_ylabel("y")
    cb = fig.colorbar(im, ax=ax, shrink=0.85); cb.set_label(clabel)
    fig.tight_layout()
    # Use the series index, not rounded physical time.  Fractional write
    # intervals (e.g. 0.1) otherwise overwrite several frames with one name.
    fpng = os.path.join(out_dir, f"{field}_{frame_index:04d}.png")
    fig.savefig(fpng, dpi=90)
    plt.close(fig)
    frames.append(fpng)

imgs = [Image.open(f) for f in frames]
gif = os.path.join(out_dir, f"{field}.gif")
imgs[0].save(gif, save_all=True, append_images=imgs[1:], duration=120, loop=0)
print(f"wrote {len(frames)} frames + {gif}")
