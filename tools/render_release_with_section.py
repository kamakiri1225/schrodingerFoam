#!/usr/bin/env python3
"""量子・古典密度と、y=0断面での横方向の広がりを比較するGIFを作る。

量子密度はschrodingerFoamのVTK出力から読み込む。古典密度は、量子と同じ
初期位置分布に等方的な速度幅 sigma_vx=sigma_vy=1 を与え、解放後に

    sigma_j(t)^2 = sigma_j0^2 + sigma_vj^2 tau^2

で自由飛行させた解析的ガウス集団である。これは同一トラップ内の正準熱平衡
ではなく、初期形状を揃えた比較用モデルである。

使い方:
    python3 tools/render_release_with_section.py <VTK_dir> <output_dir>
"""

import glob
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
import vtk
from vtk.util.numpy_support import vtk_to_numpy


vtk_dir = sys.argv[1]
out_dir = sys.argv[2]
frames_dir = os.path.join(out_dir, "section_frames")
os.makedirs(frames_dir, exist_ok=True)

omega_x = 0.03125
omega_y = 2.0
release_time = 2.0
thermal_speed = 1.0
display_limit = 16.0
sigma_x0 = 1.0/np.sqrt(2.0*omega_x)
sigma_y0 = 1.0/np.sqrt(2.0*omega_y)

series = glob.glob(os.path.join(vtk_dir, "*.vtm.series"))[0]
with open(series, encoding="utf-8") as stream:
    entries = json.load(stream)["files"]


def read_density(vtm_name):
    """OpenFOAMの実部・虚部から量子密度 |psi|^2 を読み出す。"""
    vtu = os.path.join(vtk_dir, vtm_name.replace(".vtm", ""), "internal.vtu")
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(vtu)
    reader.Update()
    mesh = reader.GetOutput()

    centers = vtk.vtkCellCenters()
    centers.SetInputData(mesh)
    centers.Update()
    points = vtk_to_numpy(centers.GetOutput().GetPoints().GetData())
    re = vtk_to_numpy(mesh.GetCellData().GetArray("Psire"))
    im = vtk_to_numpy(mesh.GetCellData().GetArray("Psiim"))
    return points, re*re + im*im


points0, _ = read_density(entries[0]["name"])
xs = np.unique(np.round(points0[:, 0], 7))
ys = np.unique(np.round(points0[:, 1], 7))
ix = np.searchsorted(xs, np.round(points0[:, 0], 7))
iy = np.searchsorted(ys, np.round(points0[:, 1], 7))
xg, yg = np.meshgrid(xs, ys)
near_y0 = np.isclose(np.abs(ys), np.min(np.abs(ys)))
near_x0 = np.isclose(np.abs(xs), np.min(np.abs(xs)))


def gaussian_density(sx, sy):
    """標準偏差sx、syの規格化された2次元ガウス密度を返す。"""
    return (
        np.exp(-0.5*((xg/sx)**2 + (yg/sy)**2))
        /(2.0*np.pi*sx*sy)
    )


def moments(rho):
    """密度の重み付き2次モーメントからsigma_x、sigma_yを求める。"""
    total = np.sum(rho)
    mx = np.sum(rho*xg)/total
    my = np.sum(rho*yg)/total
    sx = np.sqrt(np.sum(rho*(xg - mx)**2)/total)
    sy = np.sqrt(np.sum(rho*(yg - my)**2)/total)
    return sx, sy


frame_paths = []
for frame_index, entry in enumerate(entries):
    t = float(entry["time"])
    _, values = read_density(entry["name"])
    quantum = np.zeros((len(ys), len(xs)))
    quantum[iy, ix] = values

    # 古典粒子は解放後 x=x0+vx*tau と自由飛行する。独立な位置・速度の
    # 分散を加算して、各時刻のガウス密度を解析的に作る。
    tau = max(t - release_time, 0.0)
    sx_classical = np.sqrt(sigma_x0**2 + (thermal_speed*tau)**2)
    sy_classical = np.sqrt(sigma_y0**2 + (thermal_speed*tau)**2)
    classical = gaussian_density(sx_classical, sy_classical)
    sx_quantum, sy_quantum = moments(quantum)

    # x=0、y=0はセル中心に存在しないため、それぞれ最寄り2列を平均する。
    quantum_cut = np.mean(quantum[near_y0, :], axis=0)
    classical_cut = np.mean(classical[near_y0, :], axis=0)
    quantum_vertical_cut = np.mean(quantum[:, near_x0], axis=1)
    classical_vertical_cut = np.mean(classical[:, near_x0], axis=1)
    quantum_cut /= max(float(np.max(quantum_cut)), np.finfo(float).tiny)
    classical_cut /= max(float(np.max(classical_cut)), np.finfo(float).tiny)
    quantum_vertical_cut /= max(
        float(np.max(quantum_vertical_cut)), np.finfo(float).tiny
    )
    classical_vertical_cut /= max(
        float(np.max(classical_vertical_cut)), np.finfo(float).tiny
    )

    trap_on = t < release_time
    potential_cut = 0.5*np.square(omega_x*xs) if trap_on else np.zeros_like(xs)
    state = "TRAP ON" if trap_on else "RELEASED: V=0"

    fig = plt.figure(figsize=(13.2, 7.6))
    layout = fig.add_gridspec(2, 1, height_ratios=(3.0, 1.45), hspace=0.30)
    top_layout = layout[0].subgridspec(
        1, 4, width_ratios=(0.42, 1.0, 0.42, 1.0), wspace=0.08
    )

    # x=0断面を各密度図の左横へ置く。密度0を画像側にするため横軸を反転し、
    # 曲線のy座標と2次元密度図のy座標がそのまま対応するようにする。
    quantum_vertical_ax = fig.add_subplot(top_layout[0, 0])
    quantum_ax = fig.add_subplot(top_layout[0, 1], sharey=quantum_vertical_ax)
    classical_vertical_ax = fig.add_subplot(
        top_layout[0, 2], sharey=quantum_vertical_ax
    )
    classical_ax = fig.add_subplot(
        top_layout[0, 3], sharex=quantum_ax, sharey=quantum_vertical_ax
    )

    for ax, cut, color, sy_value in (
        (quantum_vertical_ax, quantum_vertical_cut, "#1565c0", sy_quantum),
        (classical_vertical_ax, classical_vertical_cut, "#ef6c00", sy_classical),
    ):
        ax.plot(cut, ys, color=color, lw=2.3)
        ax.fill_betweenx(ys, 0, cut, color=color, alpha=0.14)
        ax.set_xlim(1.08, 0)
        ax.set_ylim(-display_limit, display_limit)
        ax.set_xlabel("line peak")
        ax.set_title(f"x=0 section\nsigma_y={sy_value:.2f}", fontsize=9.5)
        ax.grid(alpha=0.22)
    quantum_vertical_ax.set_ylabel("y")
    classical_vertical_ax.tick_params(labelleft=False)

    for ax, rho, title in (
        (
            quantum_ax,
            quantum,
            f"Quantum: OpenFOAM\nsigma=({sx_quantum:.2f}, {sy_quantum:.2f})",
        ),
        (
            classical_ax,
            classical,
            "Classical: same density + isotropic velocities\n"
            f"sigma=({sx_classical:.2f}, {sy_classical:.2f})",
        ),
    ):
        shown = rho/max(float(np.max(rho)), np.finfo(float).tiny)
        image = ax.imshow(
            shown,
            origin="lower",
            extent=[xs.min(), xs.max(), ys.min(), ys.max()],
            cmap="turbo",
            vmin=0,
            vmax=1,
            interpolation="bilinear",
        )
        ax.axhline(0, color="white", lw=0.8, ls="--", alpha=0.8)
        ax.set_xlim(-display_limit, display_limit)
        ax.set_ylim(-display_limit, display_limit)
        ax.set_aspect("equal")
        ax.set_xlabel("x")
        ax.set_title(title, fontsize=10.5)
    quantum_ax.tick_params(labelleft=False)
    classical_ax.tick_params(labelleft=False)
    colorbar = fig.colorbar(image, ax=[quantum_ax, classical_ax], fraction=0.026, pad=0.025)
    colorbar.set_label("density / instantaneous peak")

    section_ax = fig.add_subplot(layout[1])
    section_ax.plot(
        xs, quantum_cut, color="#1565c0", lw=2.5,
        label=fr"Quantum y=0   $\sigma_x={sx_quantum:.2f}$",
    )
    section_ax.plot(
        xs, classical_cut, color="#ef6c00", lw=2.5,
        label=fr"Classical y=0   $\sigma_x={sx_classical:.2f}$",
    )
    section_ax.fill_between(xs, quantum_cut, color="#64b5f6", alpha=0.15)
    section_ax.set_xlim(-display_limit, display_limit)
    section_ax.set_ylim(0, 1.08)
    section_ax.set_xlabel("x at y=0")
    section_ax.set_ylabel("density / line peak")
    section_ax.grid(alpha=0.22)

    potential_ax = section_ax.twinx()
    potential_ax.plot(xs, potential_cut, "--", color="#c62828", lw=1.6, label=r"$V(x,0,t)$")
    potential_max = float(np.max(0.5*np.square(omega_x*xs)))
    potential_ax.set_ylim(0, max(1.1*potential_max, 1e-4))
    potential_ax.set_ylabel("effective potential V", color="#c62828")
    potential_ax.tick_params(axis="y", labelcolor="#c62828")

    lines = section_ax.get_lines() + potential_ax.get_lines()
    section_ax.legend(lines, [line.get_label() for line in lines], loc="upper right", fontsize=9)

    fig.suptitle(f"Quantum vs classical release   t={t:g}   |   {state}", fontsize=14)
    fig.subplots_adjust(left=0.06, right=0.92, bottom=0.08, top=0.88)

    frame_path = os.path.join(frames_dir, f"section_{frame_index:04d}.png")
    # 全フレームで同じキャンバス寸法を保ち、GIF結合時の寸法不一致を防ぐ。
    fig.savefig(frame_path, dpi=100)
    plt.close(fig)
    frame_paths.append(frame_path)

images = [Image.open(path) for path in frame_paths]
gif_path = os.path.join(out_dir, "density_with_potential_section.gif")
images[0].save(
    gif_path,
    save_all=True,
    append_images=images[1:],
    duration=120,
    loop=0,
    optimize=True,
)
print(f"wrote {len(frame_paths)} frames + {gif_path}")
