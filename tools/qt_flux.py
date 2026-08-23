#!/usr/bin/env python3
"""エネルギーカスケードの向きを、非圧縮運動エネルギーのフラックス Pi(k) で判定する.

  python3 tools/qt_flux.py <VTK_dir> [out.png] [--rho-frac 0.0] [--D 0.5]
                           [--i0 A --i1 B] [--C 1.0] [--eps-out FILE]

原理（Numasato–Tsubota–L'vov, k 空間のエネルギー収支）：
波数 k を通って小スケール側へ流れる非圧縮運動エネルギーのフラックスは
     Pi(k, t) = - d/dt ∫_0^k E^i(k') dk'      （k 以下の非圧縮エネルギーの減り）
   * Pi(k) > 0（慣性領域）→ 大→小スケールの **順カスケード（direct）**
   * Pi(k) < 0            → 小→大スケールの **逆カスケード（inverse）**
既定では全隣接時刻ペアで Pi(k,t) を求めて時間平均し（減衰乱流でも安定）、慣性帯の
プラトーを エネルギー流束 ε とみなす。Kolmogorov–Obukhov 係数 A=C·ε^{2/3} も出力し、
これを qt_spectrum* の k^{-5/3} 線の「物理的に決まる係数」として使える。
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
plt.rcParams["axes.unicode_minus"] = False

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qt_analysis import read_vtk_fields, energy_spectra


def spectrum_at(vtk_dir, idx, D, rho_frac):
    f, xs, ys, meta = read_vtk_fields(vtk_dir, idx)
    sp = energy_spectra(f["Psire"], f["Psiim"], meta, D=D, rho_frac=rho_frac)
    return sp["k"], sp["Ei"], meta["time"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("vtk_dir")
    ap.add_argument("out", nargs="?", default=None)
    ap.add_argument("--i0", type=int, default=None,
                    help="時刻Aの index（--i0/--i1 両方指定でその1ペアだけ）")
    ap.add_argument("--i1", type=int, default=None, help="時刻Bの index")
    ap.add_argument("--D", type=float, default=0.5)
    ap.add_argument("--rho-frac", type=float, default=0.0,
                    help="トラップ内部だけで計算する密度しきい値（ピーク比。トラップ系は 0.3）")
    ap.add_argument("--C", type=float, default=1.0,
                    help="Kolmogorov 定数 C（A=C·ε^{2/3}）")
    ap.add_argument("--eps-out", default=None,
                    help="ε と A=C·ε^{2/3} を書き出すファイル（k^{-5/3}係数の受け渡し用）")
    a = ap.parse_args()

    series = glob.glob(os.path.join(a.vtk_dir, "*.vtm.series"))[0]
    n = len(json.load(open(series))["files"])

    # Pi(k,t) を求める時刻ペア列。既定は全隣接ペア（時間平均で安定化）。
    if a.i0 is not None and a.i1 is not None:
        pairs = [(a.i0, a.i1)]
    else:
        pairs = [(i, i + 1) for i in range(n - 1)]

    Pi_acc, k0, spans = None, None, []
    for (i0, i1) in pairs:
        k, Ei0, t0 = spectrum_at(a.vtk_dir, i0, a.D, a.rho_frac)
        _, Ei1, t1 = spectrum_at(a.vtk_dir, i1, a.D, a.rho_frac)
        dt = t1 - t0
        if dt == 0:
            continue
        dk = k[1] - k[0]
        C0 = np.cumsum(Ei0) * dk                 # ∫_0^k E^i dk'（時刻A）
        C1 = np.cumsum(Ei1) * dk                 # 〃（時刻B）
        Pi = -(C1 - C0) / dt                     # + で順カスケード
        Pi_acc = Pi if Pi_acc is None else Pi_acc + Pi
        k0 = k; spans.append((t0, t1))
    Pi_mean = Pi_acc / len(spans)

    # 慣性帯（両端を除く中央帯）で符号と ε（プラトー）を評価
    dk = k0[1] - k0[0]
    band = (k0 > 2 * dk) & (k0 < 0.4 * k0.max())
    sign = float(np.mean(Pi_mean[band]))
    eps = float(np.mean(np.clip(Pi_mean[band], 0, None))) if sign > 0 \
        else float(np.mean(Pi_mean[band]))
    A = a.C * abs(eps) ** (2.0 / 3.0)
    direct = sign > 0
    verdict = ("順カスケード direct（大→小スケール）" if direct
               else "逆カスケード inverse（小→大スケール）")
    print(f"time span   = [{spans[0][0]}, {spans[-1][1]}]  ({len(spans)} pairs)")
    print(f"慣性帯 mean Pi = {sign:.4g}  ->  {verdict}")
    print(f"ε (energy flux) = {eps:.4g}")
    print(f"A = C·ε^(2/3)   = {A:.4g}   (C={a.C})   ->  k^-5/3 係数")

    if a.eps_out:
        os.makedirs(os.path.dirname(a.eps_out) or ".", exist_ok=True)
        with open(a.eps_out, "w") as fo:
            json.dump(dict(epsilon=eps, C=a.C, A=A, direct=direct,
                           mean_Pi=sign, span=[spans[0][0], spans[-1][1]]), fo)
        print("wrote", a.eps_out)

    fig, ax = plt.subplots(figsize=(7.4, 5.0))
    ax.axhline(0, color="#888", lw=1)
    ax.plot(k0, Pi_mean, "-", color="#1f6feb", lw=2, label=r"$\overline{\Pi(k)}$")
    ax.fill_between(k0, 0, Pi_mean, where=(Pi_mean > 0), color="#1f6feb", alpha=0.15)
    ax.fill_between(k0, 0, Pi_mean, where=(Pi_mean < 0), color="#c0392b", alpha=0.15)
    ax.axhline(eps, color="k", ls="--", lw=1.3, label=fr"$\varepsilon$ = {eps:.3g}")
    ax.axvspan(k0[band][0], k0[band][-1], color="0.85", alpha=0.4, zorder=0,
               label="慣性帯")
    ax.set_xscale("log")
    ax.set_xlabel("波数  $k$")
    ax.set_ylabel(r"非圧縮エネルギー流束  $\Pi(k)$")
    tag = "順カスケード（direct）" if direct else "逆カスケード（inverse）"
    ax.set_title(f"エネルギーカスケードの向き：{tag}   "
                 f"(mean $\\Pi$={sign:.3g}, ε={eps:.3g})")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(loc="best")
    fig.tight_layout()
    out = a.out or os.path.join(os.path.dirname(a.vtk_dir.rstrip("/")),
                                "energy_flux.png")
    fig.savefig(out, dpi=130)
    print("wrote", out)


if __name__ == "__main__":
    main()
