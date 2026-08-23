#!/usr/bin/env python3
u"""非圧縮運動エネルギーの流束 ε^i(k) を Numasato–Tsubota–L'vov (PRE 81, 016303)
の Fig.6 と同じ見方で描き、エネルギーカスケードの向きを判定する。

  python3 tools/qt_flux.py <VTK_dir> [out.png] [--t0 T --t1 T] [--rho-frac 0.0]
                           [--D 0.5] [--C 1.0] [--eps-out FILE]

■ 計算式（論文 式(15)）
     ε^i(k,t) = - ∫_{k_m}^{k}  ∂E^i_kin(k',t)/∂t  dk'
  ・k 以下（大スケール側）の非圧縮運動エネルギーが「単位時間に減った量」。
  ・∂E^i/∂t は隣接時刻のスペクトル差 (E^i(t+Δt)-E^i(t))/Δt で評価。
  ・時間窓 [t0,t1] 内の全ペアで平均（既定は最初の 20% を初期過渡として除外）。

■ グラフの見方（Fig.6 と同じ）
  ・横軸=波数 k（対数）、縦軸=ε^i(k)（線形）。灰色帯=慣性領域。
  ・慣性領域で ε^i(k) > 0 → 大→小スケールの【順カスケード direct】
                    ε^i(k) < 0 → 小→大スケールの【逆カスケード inverse】
  ・慣性領域のプラトー値が エネルギー流束 ε。k^{-5/3} 係数は A=C·ε^{2/3}。

※ 論文は「初期に大スケールへエネルギーを溜める」設定なので direct。ソリトン格子や
  干渉から作る乱流は、どの時間窓・どのスケールで測るかで符号が変わり得る。
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("vtk_dir")
    ap.add_argument("out", nargs="?", default=None)
    ap.add_argument("--t0", type=float, default=None,
                    help="時間窓の開始（既定は全体の 20%% 以降＝初期過渡を除外）")
    ap.add_argument("--t1", type=float, default=None, help="時間窓の終了")
    ap.add_argument("--D", type=float, default=0.5)
    ap.add_argument("--rho-frac", type=float, default=0.0,
                    help="トラップ内部だけで計算する密度しきい値（ピーク比。トラップ系は 0.3）")
    ap.add_argument("--C", type=float, default=1.0, help="Kolmogorov 定数 C（A=C·ε^{2/3}）")
    ap.add_argument("--eps-out", default=None, help="ε,A を書き出す JSON")
    a = ap.parse_args()

    series = glob.glob(os.path.join(a.vtk_dir, "*.vtm.series"))[0]
    files = json.load(open(series))["files"]
    times = [f["time"] for f in files]
    n = len(files)
    t0 = a.t0 if a.t0 is not None else times[0] + 0.2 * (times[-1] - times[0])
    t1 = a.t1 if a.t1 is not None else times[-1]

    # 時間窓内の隣接ペアで ε^i(k) = -∂/∂t ∫_0^k E^i dk' を求めて平均
    eps_acc, k0, used = None, None, []
    for i in range(n - 1):
        if not (t0 <= times[i] <= t1):
            continue
        k, E0, _, _ = _spec(a.vtk_dir, i, a.D, a.rho_frac)
        _, E1, _, _ = _spec(a.vtk_dir, i + 1, a.D, a.rho_frac)
        dt = times[i + 1] - times[i]
        if dt == 0:
            continue
        dk = k[1] - k[0]
        eps = -(np.cumsum(E1) - np.cumsum(E0)) * dk / dt      # ε^i(k)
        eps_acc = eps if eps_acc is None else eps_acc + eps
        k0 = k
        used.append((times[i], times[i + 1]))
    eps_k = eps_acc / len(used)

    dk = k0[1] - k0[0]
    band = (k0 > 2 * dk) & (k0 < 0.4 * k0.max())              # 慣性領域
    eps_bar = float(np.mean(eps_k[band]))
    direct = eps_bar > 0
    A = a.C * abs(eps_bar) ** (2.0 / 3.0)
    verdict = "順カスケード direct（大→小）" if direct else "逆カスケード inverse（小→大）"
    print(f"time window   = [{used[0][0]}, {used[-1][1]}]  ({len(used)} pairs)")
    print(f"慣性帯 平均 ε^i = {eps_bar:.4g}  ->  {verdict}")
    print(f"|ε| (flux)     = {abs(eps_bar):.4g}   A=C·ε^(2/3)={A:.4g} (C={a.C})")
    if a.eps_out:
        os.makedirs(os.path.dirname(a.eps_out) or ".", exist_ok=True)
        json.dump(dict(epsilon=eps_bar, abseps=abs(eps_bar), C=a.C, A=A,
                       direct=direct, window=[used[0][0], used[-1][1]]),
                  open(a.eps_out, "w"))
        print("wrote", a.eps_out)

    fig, ax = plt.subplots(figsize=(7.6, 5.4))
    ax.axhline(0, color="#444", lw=1)
    ax.axvspan(k0[band][0], k0[band][-1], color="0.88", zorder=0, label="慣性領域")
    ax.plot(k0, eps_k, "-", color="#1f6feb", lw=2.2, label=r"$\varepsilon^i(k)$")
    ax.axhline(eps_bar, color="#c0392b", ls="--", lw=1.4,
               label=fr"慣性帯平均 $\varepsilon$={eps_bar:.3g}")
    ax.set_xscale("log")
    ax.set_xlabel("波数  $k$")
    ax.set_ylabel(r"非圧縮運動エネルギー流束  $\varepsilon^i(k)$")
    tag = "順カスケード（direct, 大→小）" if direct else "逆カスケード（inverse, 小→大）"
    ax.set_title(r"$\varepsilon^i(k)=-\int_{k_m}^{k}\partial_t E^i_{kin}(k')\,dk'$"
                 f"   窓 t=[{used[0][0]:g},{used[-1][1]:g}]\n{tag}", fontsize=12)
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(loc="best", fontsize=9)
    fig.tight_layout()
    out = a.out or os.path.join(os.path.dirname(a.vtk_dir.rstrip("/")),
                                "energy_flux.png")
    fig.savefig(out, dpi=130)
    print("wrote", out)


def _spec(vtk_dir, idx, D, rho_frac):
    f, xs, ys, meta = read_vtk_fields(vtk_dir, idx)
    sp = energy_spectra(f["Psire"], f["Psiim"], meta, D=D, rho_frac=rho_frac)
    return sp["k"], sp["Ei"], sp["Ec"], meta["time"]


if __name__ == "__main__":
    main()
