#!/usr/bin/env python3
"""量子渦の本数の時間減衰を描く。横軸=時間 t、縦軸=渦点数 N_v、両対数。

  python3 tools/qt_vortex_decay.py <VTK_dir> [out.png] [--rho-frac 0.0]
                                   [--ref -1] [--semilog]

渦の数え方は soliton_hist_ensamble.f90 と同じ（プラケット位相巻き数 ±1）。
参照べき t^{--ref} を点線で重ねる（既定 N_v ∝ t^{-1}）。--semilog で片対数。
CSV も同時に書き出す。
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
from qt_analysis import read_vtk_fields, detect_vortices


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("vtk_dir")
    ap.add_argument("out", nargs="?", default=None)
    ap.add_argument("--rho-frac", type=float, default=0.0,
                    help="渦検出の密度しきい値（ピーク密度に対する割合。トラップ系は 0.3）")
    ap.add_argument("--ref", type=float, default=1.0,
                    help="参照べき N_v ∝ t^{-ref}（両対数のときだけ描く）")
    ap.add_argument("--semilog", action="store_true",
                    help="片対数（縦のみ log）。既定は両対数")
    a = ap.parse_args()

    series = glob.glob(os.path.join(a.vtk_dir, "*.vtm.series"))[0]
    files = json.load(open(series))["files"]
    nt = len(files)

    ts, nv = [], []
    for idx in range(nt):
        f, xs, ys, meta = read_vtk_fields(a.vtk_dir, idx)
        rho = f["Psire"]**2 + f["Psiim"]**2
        rho_min = a.rho_frac * float(rho.max())
        _, (npl, nmi, ntot, net) = detect_vortices(
            f["Psire"], f["Psiim"], rho_min=rho_min)
        ts.append(float(meta["time"])); nv.append(npl + nmi)
        print(f"t={meta['time']:>8}  N_v={npl+nmi:5d}  (+{npl}/-{nmi})")

    ts = np.array(ts); nv = np.array(nv, float)

    out = a.out or "graph/vortex_decay.png"
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    csv = os.path.splitext(out)[0] + ".csv"
    np.savetxt(csv, np.column_stack([ts, nv]), fmt="%.6g",
               header="time  N_vortex", comments="")

    fig, ax = plt.subplots(figsize=(7.2, 5.4))
    pos = nv > 0
    if a.semilog:
        ax.semilogy(ts[pos], nv[pos], "o-", color="#1f6feb", lw=1.8, ms=5)
    else:
        tp = ts > 0
        ax.loglog(ts[pos & tp], nv[pos & tp], "o-",
                  color="#1f6feb", lw=1.8, ms=5, label=r"$N_v(t)$")
        # 参照べき t^{-ref}（減衰域にざっくり合わせる）
        tt = ts[pos & tp]
        if tt.size > 3:
            band = tt >= tt[len(tt) // 3]
            tr = tt[band]
            amp = np.median(nv[pos & tp][band] * tr**a.ref)
            ax.loglog(tr, amp * tr**(-a.ref), "k--", lw=1.5,
                      label=rf"$t^{{-{a.ref:g}}}$")
        ax.legend()
    ax.set_xlabel("時間  $t$")
    ax.set_ylabel("量子渦の本数  $N_v$")
    ax.set_title("量子渦の本数の時間減衰")
    ax.grid(True, which="both", alpha=0.25)
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    print("wrote", out, "and", csv)


if __name__ == "__main__":
    main()
