#!/usr/bin/env python3
"""量子乱流の運動エネルギースペクトル E^i(k), E^c(k) と 量子渦点(±1) を出す.

  python3 tools/qt_spectrum.py <VTK_dir> [out.png] [--index N] [--D 0.5]

<VTK_dir> は foamToVTK 出力（Psire, Psiim を含む）。--index 省略で最終時刻。
非圧縮スペクトルに Kolmogorov–Obukhov  k^{-5/3}  の目安線を重ねる。
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
from qt_analysis import read_vtk_fields, energy_spectra, detect_vortices


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("vtk_dir")
    ap.add_argument("out", nargs="?", default=None)
    ap.add_argument("--index", type=int, default=None)
    ap.add_argument("--D", type=float, default=0.5)
    a = ap.parse_args()

    flds, xs, ys, meta = read_vtk_fields(a.vtk_dir, a.index)
    Re, Im = flds["Psire"], flds["Psiim"]

    sp = energy_spectra(Re, Im, meta, D=a.D)
    winding, (npl, nmi, ntot, net) = detect_vortices(Re, Im)

    print(f"time            = {meta['time']}")
    print(f"grid            = {meta['nx']} x {meta['ny']}  (dx={meta['dx']})")
    print(f"E^i_kin (incomp)= {sp['Ei_tot']:.6g}")
    print(f"E^c_kin (comp)  = {sp['Ec_tot']:.6g}")
    print(f"vortices        = +{npl}  -{nmi}  (total {ntot}, net {net})")

    k, Ei, Ec = sp["k"], sp["Ei"], sp["Ec"]
    m = k > 0
    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    ax.loglog(k[m], Ei[m], "-", color="#1f6feb", lw=2, label=r"$E^i_{kin}(k)$ 非圧縮")
    ax.loglog(k[m], Ec[m], "-", color="#c0392b", lw=2, label=r"$E^c_{kin}(k)$ 圧縮")
    # k^-5/3 の目安線（慣性領域に合わせて縦位置を調整）
    kk = k[m]
    band = (kk > 2 * k[m].min()) & (kk < 0.3 * k[m].max())
    if band.any():
        kref = kk[band]
        amp = np.median(Ei[m][band] * kref**(5 / 3))
        ax.loglog(kref, amp * kref**(-5 / 3), "k--", lw=1.4,
                  label=r"$k^{-5/3}$（Kolmogorov–Obukhov）")
    ax.set_xlabel("wavenumber  $k$")
    ax.set_ylabel("kinetic energy spectrum")
    ax.set_title(f"2D quantum turbulence spectra   t = {meta['time']}   "
                 f"(vortices +{npl}/-{nmi})")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    fig.tight_layout()
    out = a.out or os.path.join(os.path.dirname(a.vtk_dir.rstrip("/")),
                                "energy_spectrum.png")
    fig.savefig(out, dpi=130)
    print("wrote", out)


if __name__ == "__main__":
    main()
