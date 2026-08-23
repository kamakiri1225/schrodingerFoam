#!/usr/bin/env python3
"""非圧縮運動エネルギースペクトル E^i(k) の時間発展を1枚に重ねて描く。
横軸=波数 k、縦軸=非圧縮成分の運動エネルギー、k^{-5/3} を点線で。

  python3 tools/qt_spectrum_time.py <VTK_dir> [out.png] [--n 7] [--comp i|c]

<VTK_dir> は foamToVTK 出力（Psire,Psiim を含む）。--n で重ねる時刻数。
--comp i(非圧縮=既定) / c(圧縮)。
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
_jp = os.path.expanduser("~/.fonts/NotoSansCJKjp-Regular.otf")
if os.path.exists(_jp):
    _fm.fontManager.addfont(_jp)
    plt.rcParams["font.family"] = _fm.FontProperties(fname=_jp).get_name()

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qt_analysis import read_vtk_fields, energy_spectra


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("vtk_dir")
    ap.add_argument("out", nargs="?", default=None)
    ap.add_argument("--n", type=int, default=7, help="重ねる時刻数")
    ap.add_argument("--comp", choices=["i", "c"], default="i")
    ap.add_argument("--D", type=float, default=0.5)
    ap.add_argument("--rho-frac", type=float, default=0.0,
                    help="トラップ内部だけで計算する密度しきい値（ピーク比。トラップ系は 0.3）")
    a = ap.parse_args()

    series = glob.glob(os.path.join(a.vtk_dir, "*.vtm.series"))[0]
    files = json.load(open(series))["files"]
    ntot = len(files)
    idxs = np.unique(np.linspace(0, ntot - 1, a.n).astype(int))

    cmap = plt.cm.viridis
    fig, ax = plt.subplots(figsize=(7.6, 5.8))
    key = "Ei" if a.comp == "i" else "Ec"
    lbl = "非圧縮" if a.comp == "i" else "圧縮"
    kkeep = None
    for j, idx in enumerate(idxs):
        f, xs, ys, meta = read_vtk_fields(a.vtk_dir, int(idx))
        sp = energy_spectra(f["Psire"], f["Psiim"], meta, D=a.D,
                            rho_frac=a.rho_frac)
        k, E = sp["k"], sp[key]
        m = (k > 0) & (E > 0)
        c = cmap(j / max(len(idxs) - 1, 1))
        ax.loglog(k[m], E[m], "-", color=c, lw=1.8,
                  label=f"t = {meta['time']}")
        kkeep = k[m]

    # k^-5/3 の点線（慣性帯にざっくり合わせる）
    if kkeep is not None:
        band = (kkeep > 2 * kkeep.min()) & (kkeep < 0.35 * kkeep.max())
        kref = kkeep[band]
        # 最終時刻のスペクトルに合わせて縦位置を決める
        Elast = sp[key][(sp["k"] > 0)]
        klast = sp["k"][(sp["k"] > 0)]
        amp = np.median(np.interp(kref, klast, Elast) * kref**(5 / 3))
        ax.loglog(kref, amp * kref**(-5 / 3), "k--", lw=1.6,
                  label=r"$k^{-5/3}$")

    ax.set_xlabel("波数  $k$")
    ax.set_ylabel(f"{lbl}運動エネルギースペクトル  $E^{a.comp}(k)$")
    ax.set_title(f"{lbl}運動エネルギースペクトルの時間発展")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=9, ncol=2)
    fig.tight_layout()
    out = a.out or "graph/spectrum_evolution.png"
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    fig.savefig(out, dpi=140)
    print("wrote", out)


if __name__ == "__main__":
    main()
