#!/usr/bin/env python3
"""エネルギーカスケードの向きを、非圧縮運動エネルギーのフラックス Pi(k) で判定する.

  python3 tools/qt_flux.py <VTK_dir> [out.png] [--i0 A --i1 B] [--D 0.5]

原理（Numasato–Tsubota–L'vov, energy balance in k-space）：
波数 k を通って小スケール側へ流れる非圧縮運動エネルギーのフラックスは
     Pi(k, t) = - d/dt ∫_0^k E^i(k') dk'
       = -（k 以下に溜まっている非圧縮エネルギーの時間変化）
   * Pi(k) > 0（慣性領域）→ 大→小スケールの **順カスケード（direct）**
   * Pi(k) < 0            → 小→大スケールの **逆カスケード（inverse）**
2つの時刻 A,B のスペクトルの差分から中心差分で d/dt を評価する。
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
plt.rcParams["axes.unicode_minus"] = False

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qt_analysis import read_vtk_fields, energy_spectra


def spectrum_at(vtk_dir, idx, D):
    f, xs, ys, meta = read_vtk_fields(vtk_dir, idx)
    sp = energy_spectra(f["Psire"], f["Psiim"], meta, D=D)
    return sp["k"], sp["Ei"], meta["time"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("vtk_dir")
    ap.add_argument("out", nargs="?", default=None)
    ap.add_argument("--i0", type=int, default=None, help="時刻Aの series index")
    ap.add_argument("--i1", type=int, default=None, help="時刻Bの series index")
    ap.add_argument("--D", type=float, default=0.5)
    a = ap.parse_args()

    import glob
    import json
    series = glob.glob(os.path.join(a.vtk_dir, "*.vtm.series"))
    n = len(json.load(open(series[0]))["files"]) if series else 1
    i0 = a.i0 if a.i0 is not None else n // 3
    i1 = a.i1 if a.i1 is not None else 2 * n // 3

    k0, Ei0, t0 = spectrum_at(a.vtk_dir, i0, a.D)
    k1, Ei1, t1 = spectrum_at(a.vtk_dir, i1, a.D)
    dt = (t1 - t0)
    dk = k0[1] - k0[0]

    # 累積（k 以下の）非圧縮エネルギー  C(k) = ∫_0^k E^i dk'
    C0 = np.cumsum(Ei0) * dk
    C1 = np.cumsum(Ei1) * dk
    Pi = -(C1 - C0) / dt                      # フラックス（+ で順カスケード）

    # 慣性領域の代表値（両端を除く中央帯）で符号を判定
    m = (k0 > 2 * dk) & (k0 < 0.4 * k0.max())
    sign = np.mean(Pi[m])
    verdict = ("順カスケード direct（大→小）" if sign > 0
               else "逆カスケード inverse（小→大）")
    print(f"t0={t0}, t1={t1}, dt={dt}")
    print(f"慣性帯の平均 Pi = {sign:.4g}  ->  {verdict}")

    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    ax.axhline(0, color="#888", lw=1)
    ax.plot(k0, Pi, "-", color="#1f6feb", lw=2)
    ax.fill_between(k0, 0, Pi, where=(Pi > 0), color="#1f6feb", alpha=0.15)
    ax.fill_between(k0, 0, Pi, where=(Pi < 0), color="#c0392b", alpha=0.15)
    ax.set_xscale("log")
    ax.set_xlabel("wavenumber  $k$")
    ax.set_ylabel(r"incompressible energy flux  $\Pi(k)$")
    ax.set_title(f"cascade direction   t=[{t0}, {t1}]   "
                 f"mean $\\Pi$>0: direct" if sign > 0 else
                 f"cascade direction   t=[{t0}, {t1}]   mean $\\Pi$<0: inverse")
    ax.grid(True, which="both", alpha=0.25)
    fig.tight_layout()
    out = a.out or os.path.join(os.path.dirname(a.vtk_dir.rstrip("/")),
                                "energy_flux.png")
    fig.savefig(out, dpi=130)
    print("wrote", out)


if __name__ == "__main__":
    main()
