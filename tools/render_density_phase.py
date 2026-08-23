#!/usr/bin/env python3
"""Render density (magSqrPsi) and phase (arg Psi) side by side into one GIF.

Reads a single foamToVTK output directory (which must contain BOTH magSqrPsi and
phase, i.e. run foamToVTK -fields "(magSqrPsi phase)") and writes a 2-panel
animation: density on the left, phase on the right, sharing the same time.

Usage:
    python3 render_density_phase.py <VTK_dir> <out_gif> [tag]

<tag> is an optional label drawn in the figure title (e.g. "muShift = 1").
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

vtk_dir = sys.argv[1]
out_gif = sys.argv[2]
tag = sys.argv[3] if len(sys.argv) > 3 else ""
rho_vmax = float(sys.argv[4]) if len(sys.argv) > 4 else 1.3   # 密度カラーバー上限

frames_dir = out_gif + "_frames"          # matched by .gitignore (figures/*/*frames/)
os.makedirs(frames_dir, exist_ok=True)

series = glob.glob(os.path.join(vtk_dir, "*.vtm.series"))[0]
with open(series) as f:
    entries = json.load(f)["files"]


def read_fields(vtm_name):
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
    dens = ug.GetCellData().GetArray("magSqrPsi")
    pha = ug.GetCellData().GetArray("phase")
    dens = vtk_to_numpy(dens) if dens is not None else None
    pha = vtk_to_numpy(pha) if pha is not None else None
    return pts, dens, pha


# structured grid from the first frame
pts0, _, _ = read_fields(entries[0]["name"])
xs = np.unique(np.round(pts0[:, 0], 6))
ys = np.unique(np.round(pts0[:, 1], 6))
nx, ny = len(xs), len(ys)
ix = np.searchsorted(xs, np.round(pts0[:, 0], 6))
iy = np.searchsorted(ys, np.round(pts0[:, 1], 6))
extent = [xs.min(), xs.max(), ys.min(), ys.max()]

frames = []
for k, e in enumerate(entries):
    t = e["time"]
    pts, dens, pha = read_fields(e["name"])
    if dens is None or pha is None:
        continue
    dgrid = np.full((ny, nx), np.nan)
    pgrid = np.full((ny, nx), np.nan)
    dgrid[iy, ix] = dens
    pgrid[iy, ix] = pha

    fig, (axd, axp) = plt.subplots(1, 2, figsize=(11, 5.0))
    imd = axd.imshow(dgrid, origin="lower", extent=extent, cmap="jet",
                     vmin=0.0, vmax=rho_vmax, interpolation="bilinear")
    axd.set_title(r"density  $|\psi|^2$"); axd.set_xlabel("x"); axd.set_ylabel("y")
    fig.colorbar(imd, ax=axd, shrink=0.85)

    imp = axp.imshow(pgrid, origin="lower", extent=extent, cmap="twilight",
                     vmin=-np.pi, vmax=np.pi, interpolation="bilinear")
    axp.set_title(r"phase  $\arg\psi$"); axp.set_xlabel("x"); axp.set_ylabel("y")
    fig.colorbar(imp, ax=axp, shrink=0.85)

    suptitle = f"t = {t}"
    if tag:
        suptitle = f"{tag}    " + suptitle
    fig.suptitle(suptitle, fontsize=13)
    fig.tight_layout()
    fpng = os.path.join(frames_dir, f"f_{k:04d}.png")
    fig.savefig(fpng, dpi=85)
    plt.close(fig)
    frames.append(fpng)

imgs = [Image.open(f) for f in frames]
imgs[0].save(out_gif, save_all=True, append_images=imgs[1:], duration=120, loop=0)
print(f"wrote {len(frames)} frames + {out_gif}")
