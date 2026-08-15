#!/usr/bin/env python3
"""異方的な2次元雲を解放したときの、量子系と古典系の比較GIFを作る。

この図で扱う数式
----------------

量子パネル（左）
    このPythonスクリプト自身は量子方程式を解かない。カスタマイズしたOpenFOAM
    ソルバー ``schrodingerFoam`` が出力した Psire と Psiim を読み込み、

        rho_q(x, y, t) = |psi|^2 = Psire^2 + Psiim^2.

    を表示する。OpenFOAM側では、今回の無次元単位系において

        i d(psi)/dt = [-D nabla^2 + V(x,y,t) + g |psi|^2] psi,

    を D=1/2、g=0 として解く。解放前のポテンシャルは

        V = 1/2 [(omega_x x)^2 + (omega_y y)^2],

    であり、release_time以降は V=0 とする。初期状態には調和トラップの
    ガウス型基底状態を用いる。その密度の標準偏差は

        sigma_x0 = 1/sqrt(2 omega_x),
        sigma_y0 = 1/sqrt(2 omega_y).

古典パネル（中央・右）
    こちらはOpenFOAM計算ではなく、このPython内で評価するガウス集団の解析モデル
    である。解放後の各粒子には力が働かず、運動は

        x(t) = x0 + vx tau,   y(t) = y0 + vy tau,
        tau = max(t - release_time, 0).

    となる。初期位置と速度が互いに独立なガウス変数なら、分散は

        sigma_j(t)^2 = sigma_j0^2 + sigma_vj^2 tau^2.

    と加算される。中央パネルは量子パネルと同じ初期位置密度から出発し、等方的な
    速度幅 sigma_vx=sigma_vy=1 を与える。そのため長時間後には円形へ近づく。
    これは初期形状を揃えた対照比較であり、異方的トラップ内の正準熱平衡そのもの
    ではない。

    右パネルでは全粒子を vx=vy=0 とするので、外部ポテンシャルを除去しても密度は
    膨張しない。

パネル構成:
  1. OpenFOAMで計算した量子密度。
  2. 同じ初期位置密度と等方的な速度分散を持つ古典集団。膨張後は円形へ近づく。
  3. 解放時に全粒子の速度が0である古典集団。

使い方:
    python3 tools/render_release_comparison.py <VTK_dir> <output_dir>
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
frames_dir = os.path.join(out_dir, "comparison_frames")
os.makedirs(frames_dir, exist_ok=True)

omega_x = 0.03125
omega_y = 2.0
release_time = 2.0
sigma_x0 = 1.0/np.sqrt(2.0*omega_x)
sigma_y0 = 1.0/np.sqrt(2.0*omega_y)

# 公平なスナップショット比較として、量子雲と同じ初期位置密度を用い、解放時に
# x、y方向で同じ速度幅を与える。解放前は位置を固定しているため、これは同一
# トラップ内の正準熱平衡ではない。
thermal_speed = 1.0

series = glob.glob(os.path.join(vtk_dir, "*.vtm.series"))[0]
with open(series, encoding="utf-8") as stream:
    entries = json.load(stream)["files"]


def read_internal(vtm_name):
    vtu = os.path.join(vtk_dir, vtm_name.replace(".vtm", ""), "internal.vtu")
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(vtu)
    reader.Update()
    grid = reader.GetOutput()

    centers = vtk.vtkCellCenters()
    centers.SetInputData(grid)
    centers.Update()
    points = vtk_to_numpy(centers.GetOutput().GetPoints().GetData())

    re = vtk_to_numpy(grid.GetCellData().GetArray("Psire"))
    im = vtk_to_numpy(grid.GetCellData().GetArray("Psiim"))
    return points, re*re + im*im


points0, _ = read_internal(entries[0]["name"])
xs = np.unique(np.round(points0[:, 0], 7))
ys = np.unique(np.round(points0[:, 1], 7))
ix = np.searchsorted(xs, np.round(points0[:, 0], 7))
iy = np.searchsorted(ys, np.round(points0[:, 1], 7))
xg, yg = np.meshgrid(xs, ys)


def gaussian_density(sx, sy):
    # 規格化された2次元ガウス確率密度：
    # rho = exp[-x^2/(2 sx^2) - y^2/(2 sy^2)] / (2 pi sx sy)
    return (
        np.exp(-0.5*((xg/sx)**2 + (yg/sy)**2))
        /(2.0*np.pi*sx*sy)
    )


def moments(rho):
    total = np.sum(rho)
    mx = np.sum(rho*xg)/total
    my = np.sum(rho*yg)/total
    sx = np.sqrt(np.sum(rho*(xg - mx)**2)/total)
    sy = np.sqrt(np.sum(rho*(yg - my)**2)/total)
    return sx, sy


frame_paths = []
for frame_index, entry in enumerate(entries):
    t = float(entry["time"])
    _, qvals = read_internal(entry["name"])
    quantum = np.zeros((len(ys), len(xs)))
    quantum[iy, ix] = qvals

    # 解放後の古典自由飛行。x=x0+v*tau で、ガウス変数x0とvが独立なら
    # Var[x] = Var[x0] + Var[v]*tau^2 となる。
    tau = max(t - release_time, 0.0)
    sx_classical = np.sqrt(sigma_x0**2 + (thermal_speed*tau)**2)
    sy_classical = np.sqrt(sigma_y0**2 + (thermal_speed*tau)**2)
    thermal = gaussian_density(sx_classical, sy_classical)

    # 静止状態から解放する比較では sigma_vx=sigma_vy=0 なので、
    # sigma_j(t)=sigma_j0 のまま位置密度は変化しない。
    zero_velocity = gaussian_density(sigma_x0, sigma_y0)

    qsx, qsy = moments(quantum)
    panels = (quantum, thermal, zero_velocity)
    titles = (
        "Quantum: OpenFOAM\n"
        f"max={quantum.max():.3g}, sigma=({qsx:.2f}, {qsy:.2f})",
        "Classical: same density, isotropic velocities\n"
        f"pinned before t=2, sigma_v=(1, 1), sigma=({sx_classical:.2f}, {sy_classical:.2f})",
        "Classical: released from rest\n"
        f"pinned before t=2, sigma=({sigma_x0:.2f}, {sigma_y0:.2f})",
    )

    state = "TRAP ON" if t < release_time else "RELEASED (V=0)"
    fig, axes = plt.subplots(1, 3, figsize=(14.4, 4.8), sharex=True, sharey=True)
    for ax, rho, title in zip(axes, panels, titles):
        # 形状比較のため各パネルを瞬間最大値で規格化する。量子密度の絶対最大値は
        # 左パネルのタイトルに表示する。
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
        ax.set_xlim(-16, 16)
        ax.set_ylim(-16, 16)
        ax.set_aspect("equal")
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("x")
    axes[0].set_ylabel("y")
    fig.suptitle(f"Anisotropic cloud release:  t = {t:g}   |   {state}", fontsize=14)
    colorbar_ax = fig.add_axes([0.945, 0.18, 0.012, 0.58])
    colorbar = fig.colorbar(image, cax=colorbar_ax)
    colorbar.set_label("density / instantaneous peak")
    fig.subplots_adjust(left=0.055, right=0.925, bottom=0.10, top=0.82, wspace=0.12)

    frame_path = os.path.join(frames_dir, f"comparison_{frame_index:04d}.png")
    fig.savefig(frame_path, dpi=100)
    plt.close(fig)
    frame_paths.append(frame_path)

images = [Image.open(path) for path in frame_paths]
gif_path = os.path.join(out_dir, "quantum_vs_classical.gif")
images[0].save(
    gif_path,
    save_all=True,
    append_images=images[1:],
    duration=120,
    loop=0,
    optimize=True,
)
print(f"wrote {len(frame_paths)} frames + {gif_path}")
