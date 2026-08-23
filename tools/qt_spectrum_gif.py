#!/usr/bin/env python3
"""運動エネルギースペクトル E^i(k), E^c(k) の時間変化を GIF にする。

  python3 tools/qt_spectrum_gif.py <VTK_dir> <out.gif> [--D 0.5]

各フレームで非圧縮 E^i(k)（青）と圧縮 E^c(k)（赤）を両対数で描き、
k^{-5/3}（Kolmogorov–Obukhov）を点線で重ねる。縦軸・横軸は全時刻で固定するので、
慣性領域が立ち上がってカスケードが形成される様子がそのまま見える。
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
from qt_analysis import read_vtk_fields, energy_spectra


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("vtk_dir")
    ap.add_argument("out")
    ap.add_argument("--D", type=float, default=0.5)
    a = ap.parse_args()

    series = glob.glob(os.path.join(a.vtk_dir, "*.vtm.series"))[0]
    files = json.load(open(series))["files"]
    nt = len(files)

    # 1st pass: 全時刻のスペクトルを計算して共通の軸範囲を決める
    snaps = []
    lo, hi = np.inf, 0.0
    for idx in range(nt):
        f, xs, ys, meta = read_vtk_fields(a.vtk_dir, idx)
        sp = energy_spectra(f["Psire"], f["Psiim"], meta, D=a.D)
        snaps.append((meta["time"], sp))
        for key in ("Ei", "Ec"):
            E = sp[key][(sp["k"] > 0) & (sp[key] > 0)]
            if E.size:
                hi = max(hi, E.max()); lo = min(lo, np.percentile(E, 5))
    lo = max(lo, hi * 1e-6)
    k0 = snaps[0][1]["k"]
    kpos = k0[k0 > 0]
    kmin, kmax = kpos.min(), kpos.max()

    frames_dir = a.out + "_frames"
    os.makedirs(frames_dir, exist_ok=True)
    imgs = []
    for i, (t, sp) in enumerate(snaps):
        k = sp["k"]; m = k > 0
        fig, ax = plt.subplots(figsize=(7.2, 5.6))
        ax.loglog(k[m], sp["Ei"][m], "-", color="#1f6feb", lw=2,
                  label=r"$E^i_{kin}(k)$ 非圧縮")
        ax.loglog(k[m], sp["Ec"][m], "-", color="#c0392b", lw=2,
                  label=r"$E^c_{kin}(k)$ 圧縮")
        band = (k[m] > 2 * kmin) & (k[m] < 0.35 * kmax)
        if band.any():
            kref = k[m][band]
            amp = np.median(sp["Ei"][m][band] * kref**(5 / 3))
            if np.isfinite(amp) and amp > 0:
                ax.loglog(kref, amp * kref**(-5 / 3), "k--", lw=1.4,
                          label=r"$k^{-5/3}$")
        ax.set_xlim(kmin * 0.8, kmax * 1.2)
        ax.set_ylim(lo, hi * 1.5)
        ax.set_xlabel("波数  $k$")
        ax.set_ylabel("運動エネルギースペクトル")
        ax.set_title(f"2次元量子乱流スペクトル   t = {t}")
        ax.grid(True, which="both", alpha=0.25)
        ax.legend(loc="lower left")
        fig.tight_layout()
        p = os.path.join(frames_dir, f"f{i:04d}.png")
        fig.savefig(p, dpi=110)
        plt.close(fig)
        imgs.append(Image.open(p).convert("RGB"))

    imgs[0].save(a.out, save_all=True, append_images=imgs[1:],
                 duration=150, loop=0)
    print(f"wrote {len(imgs)} frames + {a.out}")


if __name__ == "__main__":
    main()
