#!/usr/bin/env python3
"""速度場の Helmholtz 分解（非圧縮＝渦／圧縮＝音）の時間発展を GIF にする。

  python3 tools/qt_velocity_gif.py <VTK_dir> <out.gif> [--weighted] [--D 0.5]
                                   [--rho-frac 0.0] [--vmax P]

各フレーム左から：|w^i|（非圧縮＝回転＝量子渦）／|w^c|（圧縮＝音波）。
カラーバー上限は全時刻・両成分をまたいだ共通スケール（既定は 99.5 パーセンタイル）
に固定するので、明るさの時間変化がそのままエネルギーの増減を表す。
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

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qt_analysis import (read_vtk_fields, madelung_fields, helmholtz,
                         detect_vortices, bulk_window)


def decompose(Re, Im, dx, dy, D, weighted, rho_frac=0.0):
    rho, (vx, vy), (wx, wy) = madelung_fields(Re, Im, dx, dy, D)
    fx, fy = (wx, wy) if weighted else (vx, vy)
    if rho_frac > 0.0:                       # トラップ内部だけを残す
        W = bulk_window(rho, rho_frac)
        fx, fy = fx * W, fy * W
    (Wix, Wiy), (Wcx, Wcy), _ = helmholtz(fx, fy, dx, dy)
    mi = np.hypot(np.real(np.fft.ifft2(Wix)), np.real(np.fft.ifft2(Wiy)))
    mc = np.hypot(np.real(np.fft.ifft2(Wcx)), np.real(np.fft.ifft2(Wcy)))
    return rho, mi, mc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("vtk_dir")
    ap.add_argument("out")
    ap.add_argument("--D", type=float, default=0.5)
    ap.add_argument("--weighted", action="store_true",
                    help="密度重み速度 w=sqrt(rho)v（既定は流体速度 v）")
    ap.add_argument("--rho-frac", type=float, default=0.0,
                    help="渦検出の密度しきい値（ピーク密度に対する割合）")
    ap.add_argument("--vmax", type=float, default=99.5,
                    help="共通カラーバー上限のパーセンタイル")
    a = ap.parse_args()

    series = glob.glob(os.path.join(a.vtk_dir, "*.vtm.series"))[0]
    files = json.load(open(series))["files"]
    nt = len(files)
    label = "w" if a.weighted else "v"

    # 1st pass: 共通スケールと全フレームの場を確保
    frames = []
    hi = 0.0
    for idx in range(nt):
        f, xs, ys, meta = read_vtk_fields(a.vtk_dir, idx)
        rho, mi, mc = decompose(f["Psire"], f["Psiim"],
                                meta["dx"], meta["dy"], a.D, a.weighted,
                                a.rho_frac)
        rho_min = a.rho_frac * float(rho.max())
        winding, (npl, nmi, ntot, net) = detect_vortices(
            f["Psire"], f["Psiim"], rho_min=rho_min)
        yy, xx = np.where(np.abs(winding) == 1)
        frames.append(dict(mi=mi, mc=mc, t=meta["time"], xs=xs, ys=ys,
                           xv=xs[np.clip(xx, 0, len(xs) - 1)],
                           yv=ys[np.clip(yy, 0, len(ys) - 1)],
                           sgn=winding[yy, xx], npl=npl, nmi=nmi, net=net))
        hi = max(hi, np.percentile(mi, a.vmax), np.percentile(mc, a.vmax))

    frames_dir = a.out + "_frames"
    os.makedirs(frames_dir, exist_ok=True)
    ext = [frames[0]["xs"].min(), frames[0]["xs"].max(),
           frames[0]["ys"].min(), frames[0]["ys"].max()]
    imgs = []
    for i, fr in enumerate(frames):
        fig, ax = plt.subplots(1, 2, figsize=(10.4, 5.2))
        for a_, m, ttl in (
            (ax[0], fr["mi"], f"|{label}$^i$|  非圧縮＝回転（量子渦）"),
            (ax[1], fr["mc"], f"|{label}$^c$|  圧縮（音波）"),
        ):
            im = a_.imshow(m, origin="lower", extent=ext, cmap="inferno",
                           vmin=0, vmax=hi, interpolation="bilinear")
            a_.set_title(ttl)
            a_.set_xlabel("x")
            fig.colorbar(im, ax=a_, shrink=0.82)
        ax[0].scatter(fr["xv"][fr["sgn"] > 0], fr["yv"][fr["sgn"] > 0],
                      s=7, c="cyan", marker="o", lw=0)
        ax[0].scatter(fr["xv"][fr["sgn"] < 0], fr["yv"][fr["sgn"] < 0],
                      s=7, c="lime", marker="o", lw=0)
        ax[0].set_ylabel("y")
        fig.suptitle(f"速度のHelmholtz分解  t={fr['t']}   "
                     f"渦 +{fr['npl']}/-{fr['nmi']}（net {fr['net']}）",
                     fontsize=13)
        fig.tight_layout()
        p = os.path.join(frames_dir, f"f{i:04d}.png")
        fig.savefig(p, dpi=110)
        plt.close(fig)
        imgs.append(Image.open(p).convert("RGB"))

    imgs[0].save(a.out, save_all=True, append_images=imgs[1:],
                 duration=150, loop=0)
    print(f"wrote {len(imgs)} frames + {a.out}  (common vmax={hi:.4g})")


if __name__ == "__main__":
    main()
