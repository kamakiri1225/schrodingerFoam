#!/usr/bin/env python3
"""速度場を Helmholtz 分解して「非圧縮(=回転=渦)成分」と「圧縮(=音)成分」に分け，
それぞれを可視化する後処理ツール。

  python3 tools/qt_velocity.py <VTK_dir> [out.png] [--index N] [--D 0.5] [--weighted]

既定は流体速度 v = 2D grad(theta)。--weighted で密度重み速度 w = sqrt(rho) v
（スペクトルと同じ量。渦芯で有界）。左から：|w|（全体）／|w^i|（非圧縮＝渦）／
|w^c|（圧縮＝音）を並べて出す。渦点(±1)も重ねる。
"""
import sys
import os
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as _fm
_jp = os.path.expanduser("~/.fonts/NotoSansCJKjp-Regular.otf")
if os.path.exists(_jp):
    _fm.fontManager.addfont(_jp)
    plt.rcParams["font.family"] = _fm.FontProperties(fname=_jp).get_name()

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qt_analysis import read_vtk_fields, madelung_fields, helmholtz, detect_vortices


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("vtk_dir")
    ap.add_argument("out", nargs="?", default=None)
    ap.add_argument("--index", type=int, default=None)
    ap.add_argument("--D", type=float, default=0.5)
    ap.add_argument("--weighted", action="store_true",
                    help="密度重み速度 w=sqrt(rho)v を使う（既定は流体速度 v）")
    ap.add_argument("--rho-frac", type=float, default=0.0,
                    help="渦検出の密度しきい値（ピーク密度に対する割合。トラップ系は 0.3）")
    a = ap.parse_args()

    f, xs, ys, meta = read_vtk_fields(a.vtk_dir, a.index)
    Re, Im = f["Psire"], f["Psiim"]
    dx, dy = meta["dx"], meta["dy"]
    rho, (vx, vy), (wx, wy) = madelung_fields(Re, Im, dx, dy, a.D)
    fx, fy = (wx, wy) if a.weighted else (vx, vy)
    label = "w=√ρ·v" if a.weighted else "v"

    # Helmholtz 分解
    (Wix, Wiy), (Wcx, Wcy), _ = helmholtz(fx, fy, dx, dy)
    fi_x = np.real(np.fft.ifft2(Wix)); fi_y = np.real(np.fft.ifft2(Wiy))
    fc_x = np.real(np.fft.ifft2(Wcx)); fc_y = np.real(np.fft.ifft2(Wcy))

    mag_tot = np.hypot(fx, fy)
    mag_i = np.hypot(fi_x, fi_y)          # 非圧縮＝回転＝渦
    mag_c = np.hypot(fc_x, fc_y)          # 圧縮＝音

    rho_min = a.rho_frac * float(rho.max())
    winding, (npl, nmi, ntot, net) = detect_vortices(Re, Im, rho_min=rho_min)
    # 渦点座標（プラケット中心）
    yy, xx = np.where(np.abs(winding) == 1)
    xv = xs[np.clip(xx, 0, len(xs) - 1)]
    yv = ys[np.clip(yy, 0, len(ys) - 1)]
    sgn = winding[yy, xx]

    ext = [xs.min(), xs.max(), ys.min(), ys.max()]
    vmax = np.percentile(mag_tot, 99)
    fig, ax = plt.subplots(1, 3, figsize=(15, 5.2))
    for a_, m, ttl in (
        (ax[0], mag_tot, f"|{label}|  全体"),
        (ax[1], mag_i, f"|{label}$^i$|  非圧縮＝回転（量子渦）"),
        (ax[2], mag_c, f"|{label}$^c$|  圧縮（音波）"),
    ):
        im = a_.imshow(m, origin="lower", extent=ext, cmap="inferno",
                       vmin=0, vmax=vmax)
        a_.set_title(ttl)
        a_.set_xlabel("x")
        fig.colorbar(im, ax=a_, shrink=0.82)
    # 渦点を非圧縮パネルに重ねる（+ 赤 / − 青）
    ax[1].scatter(xv[sgn > 0], yv[sgn > 0], s=8, c="cyan", marker="o", lw=0)
    ax[1].scatter(xv[sgn < 0], yv[sgn < 0], s=8, c="lime", marker="o", lw=0)
    ax[0].set_ylabel("y")
    fig.suptitle(f"速度のHelmholtz分解  t={meta['time']}   "
                 f"渦 +{npl}/-{nmi}（net {net}）", fontsize=13)
    fig.tight_layout()
    out = a.out or os.path.join(os.path.dirname(a.vtk_dir.rstrip("/")),
                                "velocity_decomposition.png")
    fig.savefig(out, dpi=120)
    print(f"time={meta['time']}  vortices +{npl}/-{nmi} (net {net})")
    print("wrote", out)


if __name__ == "__main__":
    main()
