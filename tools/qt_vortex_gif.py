#!/usr/bin/env python3
"""量子渦点の分布の時間発展を GIF にする（soliton_hist_ensamble.f90 と同一の巻き数法）。

  python3 tools/qt_vortex_gif.py <VTK_dir> <out.gif> [--rho-frac 0.0]
                                 [--method f90|wrap] [--bg density|phase]

密度（または位相）を背景に、巻き数 +1（赤）/ −1（青）の渦点を重ねる。数え方は
Fortran と同一（プラケット位相巻き）。--rho-frac>0 でトラップ外・低密度ノイズの
偽の巻きを除く（プラケット最大密度で判定）。タイトルに本数を表示。
"""
import sys
import os
import argparse
import glob
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as _fm
from PIL import Image

_jp = os.path.expanduser("~/.fonts/NotoSansCJKjp-Regular.otf")
if os.path.exists(_jp):
    _fm.fontManager.addfont(_jp)
    plt.rcParams["font.family"] = _fm.FontProperties(fname=_jp).get_name()
plt.rcParams["axes.unicode_minus"] = False

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qt_analysis import read_vtk_fields, detect_vortices, detect_vortices_f90


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("vtk_dir")
    ap.add_argument("out")
    ap.add_argument("--rho-frac", type=float, default=0.0,
                    help="低密度ノイズの偽渦を除く密度しきい値（ピーク比。トラップ系は 0.3）")
    ap.add_argument("--method", choices=["f90", "wrap"], default="f90",
                    help="f90=Fortran と同一の dint(.../6) 版（既定）／wrap=丸め版")
    ap.add_argument("--bg", choices=["density", "phase"], default="density")
    ap.add_argument("--pseudo-frac", type=float, default=0.0,
                    help="擬渦度による偽渦除去のしきい値（バルク n0/ξ² 比。0.02〜0.05）")
    ap.add_argument("--compare", action="store_true",
                    help="raw（左）と偽渦除去後（右）を横並びで描く")
    a = ap.parse_args()

    series = glob.glob(os.path.join(a.vtk_dir, "*.vtm.series"))[0]
    nt = len(json.load(open(series))["files"])

    def det(Re, Im, rho_min, dx, dy, pf):
        if a.method == "f90":
            return detect_vortices_f90(Re, Im, rho_min=rho_min,
                                       pseudo_frac=pf, dx=dx, dy=dy)
        return detect_vortices(Re, Im, rho_min=rho_min)

    def panel(ax, rho, Im, Re, ext, xs, ys, w, ttl):
        npl = int(np.sum(w >= 1)); nmi = int(np.sum(w <= -1))
        yy, xx = np.where(np.abs(w) >= 1)
        xv = xs[np.clip(xx, 0, len(xs) - 1)]; yv = ys[np.clip(yy, 0, len(ys) - 1)]
        sgn = w[yy, xx]
        if a.bg == "density":
            ax.imshow(rho, origin="lower", extent=ext, cmap="gray",
                      vmin=0, vmax=1.0, interpolation="bilinear")
        else:
            ax.imshow(np.arctan2(Im, Re), origin="lower", extent=ext,
                      cmap="twilight", vmin=-np.pi, vmax=np.pi)
        ax.scatter(xv[sgn > 0], yv[sgn > 0], s=14, c="red", marker="o",
                   lw=0.4, edgecolors="white", label=f"+1 ({npl})")
        ax.scatter(xv[sgn < 0], yv[sgn < 0], s=14, c="deepskyblue", marker="o",
                   lw=0.4, edgecolors="white", label=f"-1 ({nmi})")
        ax.set_xlabel("x"); ax.set_ylabel("y")
        ax.set_title(f"{ttl}  N_v={npl + nmi}（+{npl}/-{nmi}）", fontsize=11)
        ax.legend(loc="upper right", framealpha=0.6, fontsize=8)

    frames_dir = a.out + "_frames"
    os.makedirs(frames_dir, exist_ok=True)
    imgs = []
    for idx in range(nt):
        f, xs, ys, meta = read_vtk_fields(a.vtk_dir, idx)
        Re, Im = f["Psire"], f["Psiim"]
        rho = Re**2 + Im**2
        dx, dy = meta["dx"], meta["dy"]
        rho_min = a.rho_frac * float(rho.max())
        ext = [xs.min(), xs.max(), ys.min(), ys.max()]
        pf = a.pseudo_frac if a.pseudo_frac > 0 else 0.03   # compare 用の既定
        if a.compare:
            w0, _ = det(Re, Im, rho_min, dx, dy, 0.0)
            w1, _ = det(Re, Im, rho_min, dx, dy, pf)
            fig, ax = plt.subplots(1, 2, figsize=(12.2, 6.0))
            panel(ax[0], rho, Im, Re, ext, xs, ys, w0, "raw（生の巻き数）")
            panel(ax[1], rho, Im, Re, ext, xs, ys, w1,
                  f"偽渦除去後（擬渦度 {pf}）")
            fig.suptitle(f"量子渦の分布  t={meta['time']}", fontsize=13)
        else:
            w0, (npl, nmi, ntot, net) = det(Re, Im, rho_min, dx, dy, a.pseudo_frac)
            fig, ax0 = plt.subplots(figsize=(6.4, 6.0))
            panel(ax0, rho, Im, Re, ext, xs, ys, w0, f"量子渦の分布  t={meta['time']}")
        fig.tight_layout()
        p = os.path.join(frames_dir, f"f{idx:04d}.png")
        fig.savefig(p, dpi=110)
        plt.close(fig)
        imgs.append(Image.open(p).convert("RGB"))
        print(f"t={meta['time']:>6}  done")

    imgs[0].save(a.out, save_all=True, append_images=imgs[1:],
                 duration=150, loop=0)
    print(f"wrote {len(imgs)} frames + {a.out}")


if __name__ == "__main__":
    main()
