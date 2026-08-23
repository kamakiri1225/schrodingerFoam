#!/usr/bin/env python3
u"""非圧縮運動エネルギー流束 ε^i(k) の時間変化を GIF にする（各時刻の瞬間流束）。
論文式(15) と同じ  ε^i(k,t) = -∂/∂t ∫_0^k E^i(k') dk'  を隣接時刻差で評価。

  python3 tools/qt_flux_gif.py <VTK_dir> <out.gif> [--rho-frac 0.0] [--D 0.5]
                              [--rep-out FILE] [--rep-time T]

各フレーム：横軸=k（対数）, 縦軸=ε^i(k)（線形, 全時刻共通スケール）, 0線と慣性帯,
瞬間の慣性帯符号（順/逆）をタイトルに表示。--rep-out で代表時刻の静止画も出す。
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
    ap.add_argument("--rho-frac", type=float, default=0.0)
    ap.add_argument("--rep-out", default=None, help="代表時刻の静止画の出力パス")
    ap.add_argument("--rep-time", type=float, default=None,
                    help="代表時刻（既定は時系列の中央付近）")
    ap.add_argument("--ymax", type=float, default=None,
                    help="縦軸の上下限 ±ymax（既定は発達域から自動）")
    a = ap.parse_args()

    series = glob.glob(os.path.join(a.vtk_dir, "*.vtm.series"))[0]
    files = json.load(open(series))["files"]
    times = [f["time"] for f in files]
    n = len(files)

    # 各隣接ペアの瞬間流束 ε^i(k, t_mid)
    def spec(i):
        f, xs, ys, meta = read_vtk_fields(a.vtk_dir, i)
        sp = energy_spectra(f["Psire"], f["Psiim"], meta, D=a.D, rho_frac=a.rho_frac)
        return sp["k"], sp["Ei"]

    frames = []
    E_prev = None
    for i in range(n):
        k, Ei = spec(i)
        if E_prev is not None:
            dt = times[i] - times[i - 1]
            dk = k[1] - k[0]
            eps = -(np.cumsum(Ei) - np.cumsum(E_prev)) * dk / dt
            tmid = 0.5 * (times[i] + times[i - 1])
            frames.append((tmid, k, eps))
        E_prev = Ei

    k0 = frames[0][1]
    dk = k0[1] - k0[0]
    band = (k0 > 2 * dk) & (k0 < 0.4 * k0.max())
    # 縦軸スケールは発達域（最初の20%を除く）に合わせる。初期過渡の巨大値は振り切る。
    skip = max(1, int(0.2 * len(frames)))
    if a.ymax is not None:
        hi = a.ymax
    else:
        hi = max(np.percentile(np.abs(e), 95) for _, _, e in frames[skip:])
        hi = max(hi, 1e-6)

    def draw(ax, t, k, eps):
        ax.axhline(0, color="#444", lw=1)
        ax.axvspan(k[band][0], k[band][-1], color="0.88", zorder=0, label="慣性領域")
        ax.plot(k, eps, "-", color="#1f6feb", lw=2.2, label=r"$\varepsilon^i(k)$")
        s = float(np.mean(eps[band]))
        tag = "順 direct（大→小）" if s > 0 else "逆 inverse（小→大）"
        ax.axhline(s, color="#c0392b", ls="--", lw=1.2, label=f"慣性帯平均={s:.2g}")
        ax.set_xscale("log"); ax.set_xlim(k[k > 0].min(), k.max())
        ax.set_ylim(-hi * 1.1, hi * 1.1)
        ax.set_xlabel("波数  $k$"); ax.set_ylabel(r"$\varepsilon^i(k)$")
        ax.set_title(f"非圧縮エネルギー流束  t={t:g}   [{tag}]")
        ax.grid(True, which="both", alpha=0.25); ax.legend(loc="upper left", fontsize=9)

    frames_dir = a.out + "_frames"
    os.makedirs(frames_dir, exist_ok=True)
    imgs = []
    for j, (t, k, eps) in enumerate(frames):
        fig, ax = plt.subplots(figsize=(7.4, 5.2)); draw(ax, t, k, eps)
        fig.tight_layout(); p = os.path.join(frames_dir, f"f{j:04d}.png")
        fig.savefig(p, dpi=110); plt.close(fig)
        imgs.append(Image.open(p).convert("RGB"))
    imgs[0].save(a.out, save_all=True, append_images=imgs[1:], duration=180, loop=0)
    print(f"wrote {len(imgs)} frames + {a.out}")

    if a.rep_out:
        # 代表グラフ＝発達域の複数時刻を重ね描き（順カスケードが立つ様子）
        dev = frames[skip:]
        idxs = np.unique(np.linspace(0, len(dev) - 1, 6).astype(int))
        cmap = plt.cm.viridis
        fig, ax = plt.subplots(figsize=(7.6, 5.4))
        ax.axhline(0, color="#444", lw=1)
        ax.axvspan(k0[band][0], k0[band][-1], color="0.9", zorder=0, label="慣性領域")
        for j, ii in enumerate(idxs):
            t, k, eps = dev[ii]
            ax.plot(k, eps, "-", lw=1.8, color=cmap(j / max(len(idxs) - 1, 1)),
                    label=f"t={t:g}")
        ax.set_xscale("log"); ax.set_xlim(k0[k0 > 0].min(), k0.max())
        ax.set_ylim(-hi * 1.1, hi * 1.3)
        ax.set_xlabel("波数  $k$")
        ax.set_ylabel(r"非圧縮エネルギー流束  $\varepsilon^i(k)$")
        ax.set_title(r"$\varepsilon^i(k)=-\int_{k_m}^{k}\partial_t E^i_{kin}\,dk'$"
                     "   複数時刻（慣性帯で正＝順カスケード）")
        ax.grid(True, which="both", alpha=0.25); ax.legend(fontsize=9, ncol=2)
        fig.tight_layout(); fig.savefig(a.rep_out, dpi=140); plt.close(fig)
        print(f"wrote multi-time representative -> {a.rep_out}")


if __name__ == "__main__":
    main()
