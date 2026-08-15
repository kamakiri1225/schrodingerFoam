#!/usr/bin/env python3
"""Compare one field from two foamToVTK directories side by side into one GIF.

Handy for showing the SAME physics rendered two ways, e.g. the phase field with
and without chemical-potential subtraction (flickering vs frozen background).

Usage:
    python3 render_compare.py <VTK_A> <VTK_B> <field> <out_gif> <labelA> <labelB>

Both directories must have matching time series (same number of frames).
"""
import sys
import os
import glob
import json

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

vtk_a, vtk_b, field, out_gif = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
label_a = sys.argv[5] if len(sys.argv) > 5 else "A"
label_b = sys.argv[6] if len(sys.argv) > 6 else "B"

frames_dir = out_gif + "_frames"          # matched by .gitignore (figures/*/*frames/)
os.makedirs(frames_dir, exist_ok=True)

if field == "magSqrPsi":
    vmin, vmax, cmap = 0.0, 1.3, "jet"
else:
    vmin, vmax, cmap = -np.pi, np.pi, "twilight"


def load(vtk_dir):
    series = glob.glob(os.path.join(vtk_dir, "*.vtm.series"))[0]
    with open(series) as f:
        return json.load(f)["files"]


def read(vtk_dir, vtm_name):
    vtu = os.path.join(vtk_dir, vtm_name.replace(".vtm", ""), "internal.vtu")
    r = vtk.vtkXMLUnstructuredGridReader()
    r.SetFileName(vtu)
    r.Update()
    ug = r.GetOutput()
    cc = vtk.vtkCellCenters()
    cc.SetInputData(ug)
    cc.Update()
    pts = vtk_to_numpy(cc.GetOutput().GetPoints().GetData())
    arr = ug.GetCellData().GetArray(field)
    return pts, (vtk_to_numpy(arr) if arr is not None else None)


ent_a, ent_b = load(vtk_a), load(vtk_b)
n = min(len(ent_a), len(ent_b))

pts0, _ = read(vtk_a, ent_a[0]["name"])
xs = np.unique(np.round(pts0[:, 0], 6))
ys = np.unique(np.round(pts0[:, 1], 6))
nx, ny = len(xs), len(ys)
ix = np.searchsorted(xs, np.round(pts0[:, 0], 6))
iy = np.searchsorted(ys, np.round(pts0[:, 1], 6))
extent = [xs.min(), xs.max(), ys.min(), ys.max()]

frames = []
for k in range(n):
    t = ent_a[k]["time"]
    _, va = read(vtk_a, ent_a[k]["name"])
    _, vb = read(vtk_b, ent_b[k]["name"])
    if va is None or vb is None:
        continue
    ga = np.full((ny, nx), np.nan); ga[iy, ix] = va
    gb = np.full((ny, nx), np.nan); gb[iy, ix] = vb

    fig, (axa, axb) = plt.subplots(1, 2, figsize=(11, 5.0))
    for ax, g, lab in ((axa, ga, label_a), (axb, gb, label_b)):
        im = ax.imshow(g, origin="lower", extent=extent, cmap=cmap,
                       vmin=vmin, vmax=vmax, interpolation="bilinear")
        ax.set_title(lab); ax.set_xlabel("x"); ax.set_ylabel("y")
        fig.colorbar(im, ax=ax, shrink=0.85)
    fig.suptitle(f"{field}   t = {t}", fontsize=13)
    fig.tight_layout()
    fpng = os.path.join(frames_dir, f"c_{k:04d}.png")
    fig.savefig(fpng, dpi=85)
    plt.close(fig)
    frames.append(fpng)

imgs = [Image.open(f) for f in frames]
imgs[0].save(out_gif, save_all=True, append_images=imgs[1:], duration=120, loop=0)
print(f"wrote {len(frames)} frames + {out_gif}")
