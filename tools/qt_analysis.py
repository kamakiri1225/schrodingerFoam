#!/usr/bin/env python3
"""2次元圧縮性量子乱流の解析（Numasato–Tsubota–L'vov, PRE 81, 016303 (2010)）.

foamToVTK 出力（Psire, Psiim）から
  * 速度場          v = 2D grad(theta) = 2D j/rho
  * 密度重み速度    w = sqrt(rho) v = 2D j/sqrt(rho)      （核でも有界）
  * Helmholtz 分解  w = w^i(非圧縮/solenoidal) + w^c(圧縮/potential)
  * 運動エネルギースペクトル  E^i(k), E^c(k)
  * 量子渦点の検出（位相の巻き数 = ±1 の位相欠陥）
を計算する共有モジュール。CLI は qt_spectrum.py / qt_flux.py。

前提：周期境界（cyclic）。j = Psire*grad(Psiim) - Psiim*grad(Psire)（hbar/m=2D）。
"""
import os
import glob
import json
import numpy as np


# ----------------------------------------------------------------------------
# 1. foamToVTK の内部場を構造格子に読み込む
# ----------------------------------------------------------------------------
def read_vtk_fields(vtk_dir, index=None, fields=("Psire", "Psiim")):
    """foamToVTK ディレクトリから指定時刻の場を (ny, nx) 配列で返す.

    vtk_dir : 例  run/<case>/VTK
    index   : None なら最終時刻。整数なら series の何番目か。
    return  : dict(field -> 2D array), x(1D), y(1D), meta(dict)
    """
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy

    series = glob.glob(os.path.join(vtk_dir, "*.vtm.series"))
    if series:
        entries = json.load(open(series[0]))["files"]
        idx = (len(entries) - 1) if index is None else index
        name = entries[idx]["name"].replace(".vtm", "")
        t = entries[idx]["time"]
        vtu = os.path.join(vtk_dir, name, "internal.vtu")
    else:                                   # 単一 .vtu を直接渡した場合
        vtu = vtk_dir
        t = None
    r = vtk.vtkXMLUnstructuredGridReader()
    r.SetFileName(vtu)
    r.Update()
    ug = r.GetOutput()
    cc = vtk.vtkCellCenters()
    cc.SetInputData(ug)
    cc.Update()
    pts = vtk_to_numpy(cc.GetOutput().GetPoints().GetData())
    xs = np.unique(np.round(pts[:, 0], 6))
    ys = np.unique(np.round(pts[:, 1], 6))
    nx, ny = len(xs), len(ys)
    ix = np.searchsorted(xs, np.round(pts[:, 0], 6))
    iy = np.searchsorted(ys, np.round(pts[:, 1], 6))

    out = {}
    for f in fields:
        arr = ug.GetCellData().GetArray(f)
        if arr is None:
            raise KeyError(f"field {f} not found in {vtu}")
        g = np.full((ny, nx), np.nan)
        g[iy, ix] = vtk_to_numpy(arr)
        out[f] = g
    meta = dict(nx=nx, ny=ny, dx=float(xs[1] - xs[0]), dy=float(ys[1] - ys[0]),
                Lx=float(nx * (xs[1] - xs[0])), Ly=float(ny * (ys[1] - ys[0])),
                time=t)
    return out, xs, ys, meta


# ----------------------------------------------------------------------------
# 2. 場の量（周期境界の中心差分）
# ----------------------------------------------------------------------------
def _ddx(f, dx):                            # x = axis 1
    return (np.roll(f, -1, 1) - np.roll(f, 1, 1)) / (2 * dx)


def _ddy(f, dy):                            # y = axis 0
    return (np.roll(f, -1, 0) - np.roll(f, 1, 0)) / (2 * dy)


def madelung_fields(a, b, dx, dy, D=0.5, rho_floor=1e-6):
    """a=Psire, b=Psiim から rho, 速度 v, 密度重み速度 w を返す（各 (vx,vy)）."""
    rho = a * a + b * b
    # 確率カレント  j = a grad b - b grad a   （hbar/m = 2D）
    jx = 2 * D * (a * _ddx(b, dx) - b * _ddx(a, dx))
    jy = 2 * D * (a * _ddy(b, dy) - b * _ddy(a, dy))
    rho_c = np.maximum(rho, rho_floor)
    vx, vy = jx / rho_c, jy / rho_c         # 流体速度 v = j/rho
    wx = np.sqrt(rho_c) * vx                # 密度重み速度 w = sqrt(rho) v
    wy = np.sqrt(rho_c) * vy
    return rho, (vx, vy), (wx, wy)


# ----------------------------------------------------------------------------
# 3. Helmholtz 分解 と 運動エネルギースペクトル
# ----------------------------------------------------------------------------
def helmholtz(wx, wy, dx, dy):
    """w を非圧縮 w^i(∇·=0) と 圧縮 w^c(∇×=0) に分解（Fourier 空間・eq.10d,e）."""
    ny, nx = wx.shape
    kx = 2 * np.pi * np.fft.fftfreq(nx, d=dx)
    ky = 2 * np.pi * np.fft.fftfreq(ny, d=dy)
    KX, KY = np.meshgrid(kx, ky)
    K2 = KX**2 + KY**2
    K2[0, 0] = 1.0                          # 平均モードは後で 0 に
    Wx, Wy = np.fft.fft2(wx), np.fft.fft2(wy)
    div = (KX * Wx + KY * Wy) / K2          # (k·W)/k^2
    Wcx, Wcy = KX * div, KY * div           # 圧縮（k に平行）
    Wix, Wiy = Wx - Wcx, Wy - Wcy           # 非圧縮（k に垂直）
    for A in (Wcx, Wcy, Wix, Wiy):
        A[0, 0] = 0.0
    return (Wix, Wiy), (Wcx, Wcy), (KX, KY)


def kinetic_spectrum(Wx, Wy, meta):
    """|W|^2 を |k| 殻でビン分けし，1次元スペクトル E(k) を返す（∫E dk = E_kin）."""
    ny, nx = meta["ny"], meta["nx"]
    dx, dy, Lx = meta["dx"], meta["dy"], meta["Lx"]
    N = nx * ny
    P = (np.abs(Wx)**2 + np.abs(Wy)**2)     # モードごとのエネルギー×2N/(dx dy)
    kxg = 2 * np.pi * np.fft.fftfreq(nx, d=dx)
    kyg = 2 * np.pi * np.fft.fftfreq(ny, d=dy)
    KX, KY = np.meshgrid(kxg, kyg)
    kmag = np.sqrt(KX**2 + KY**2)
    dk = 2 * np.pi / Lx                      # 殻幅
    nb = int(np.max(kmag) / dk) + 1
    kbin = (kmag / dk).astype(int)
    shell = np.bincount(kbin.ravel(), weights=P.ravel(), minlength=nb)
    kc = (np.arange(nb) + 0.5) * dk
    # E(k) = (dx dy)/(2N) * (1/dk) * Σ_shell |W|^2
    Ek = shell * (dx * dy) / (2.0 * N) / dk
    Ekin = np.sum(shell) * (dx * dy) / (2.0 * N)   # = ∫E dk
    return kc, Ek, Ekin


def bulk_window(rho, frac, soft=0.5):
    """トラップ内部（凝縮体があるところ）だけを残す滑らかな窓 W(x) を返す.

    frac : ピーク密度に対する割合。rho > frac*peak を 1、外側を tanh で滑らかに 0
           にする（ハードエッジによるスペクトル漏れを避ける）。frac<=0 なら全域 1。
    トラップ外の低密度ハロは位相がノイズで偽の速度・偽のスペクトルを生むため、
    これを掛けてから Helmholtz 分解・スペクトルを取る。
    """
    if frac <= 0:
        return np.ones_like(rho)
    thr = frac * float(np.nanmax(rho))
    return 0.5 * (1.0 + np.tanh((rho - thr) / (soft * thr + 1e-30)))


def energy_spectra(a, b, meta, D=0.5, rho_floor=1e-6, rho_frac=0.0):
    """一時刻の (E^i(k), E^c(k)) と積分エネルギーをまとめて返す.

    rho_frac>0 のときはトラップ内部だけを残す窓 W を w に掛けてから分解する
    （トラップ系で外側ハロの偽スペクトルを除く）。一様系は 0 のままでよい。
    """
    dx, dy = meta["dx"], meta["dy"]
    rho, _, (wx, wy) = madelung_fields(a, b, dx, dy, D, rho_floor)
    if rho_frac > 0.0:
        W = bulk_window(rho, rho_frac)
        wx, wy = wx * W, wy * W
    (Wix, Wiy), (Wcx, Wcy), _ = helmholtz(wx, wy, dx, dy)
    k, Ei, Ei_tot = kinetic_spectrum(Wix, Wiy, meta)
    _, Ec, Ec_tot = kinetic_spectrum(Wcx, Wcy, meta)
    return dict(k=k, Ei=Ei, Ec=Ec, Ei_tot=Ei_tot, Ec_tot=Ec_tot)


# ----------------------------------------------------------------------------
# 4. 量子渦点の検出（位相の巻き数 → ±1）
# ----------------------------------------------------------------------------
def detect_vortices(a, b, rho_min=0.0):
    """各プラケットの位相巻き数を返す（+1=渦, -1=反渦）.  ω=∇×v の量子化版.

    rho_min : この密度以下のプラケットは無視（トラップ外の低密度ハロで位相が
              ノイズになり偽の巻きを数えるのを防ぐ）。一様系では 0 でよい。

    return : charge (ny, nx) の整数格子（プラケット左下セル基準）,
             (n_plus, n_minus, n_total, net)
    """
    theta = np.arctan2(b, a)

    def wrap(d):
        return (d + np.pi) % (2 * np.pi) - np.pi

    # 反時計回りに 4 辺の位相差（各辺 wrap）を足す = 2π × 巻き数
    # 角: A=(y,x) B=(y,x+1) C=(y+1,x+1) D=(y+1,x)
    t   = theta
    txp = np.roll(t, -1, 1)                  # (y,   x+1)
    typ = np.roll(t, -1, 0)                  # (y+1, x)
    tcc = np.roll(typ, -1, 1)                # (y+1, x+1)
    d1 = wrap(txp - t)                       # A -> B  (+x)
    d2 = wrap(tcc - txp)                     # B -> C  (+y)
    d3 = wrap(typ - tcc)                     # C -> D  (-x)
    d4 = wrap(t - typ)                       # D -> A  (-y)
    winding = np.round((d1 + d2 + d3 + d4) / (2 * np.pi)).astype(int)
    if rho_min > 0.0:
        # 渦芯そのものは密度0なので、プラケットの「最大」密度で判定する
        # （ディスク内の渦は周囲が bulk 密度、トラップ外ハロは全て ~0）。
        rho = a * a + b * b
        rmax = np.maximum.reduce([rho, np.roll(rho, -1, 1),
                                  np.roll(np.roll(rho, -1, 1), -1, 0),
                                  np.roll(rho, -1, 0)])
        winding = np.where(rmax > rho_min, winding, 0)
    n_plus = int(np.sum(winding == 1) + np.sum(winding > 1))
    n_minus = int(np.sum(winding == -1) + np.sum(winding < -1))
    n_total = int(np.sum(np.abs(winding)))
    net = int(np.sum(winding))
    return winding, (n_plus, n_minus, n_total, net)


def detect_vortices_f90(a, b, rho_min=0.0):
    """soliton_hist_ensamble.f90 と同一式の忠実移植（検証用）.

    各プラケットで  Im( Σ log(conj(f_A)·f_B) )  を4辺足し，f90 と同じく /6 して
    整数化（切り捨て）する： uzudo = dint( Σ Im(log(...)) / 6 )。
      +2π → dint(1.047)=+1,  -2π → -1。
    位相差は cdlog の主値（(-π,π]）＝ np.angle(conj(f_A)*f_B) と同一。
    """
    f = a + 1j * b
    fxp = np.roll(f, -1, 1)                  # (y,   x+1)
    fyp = np.roll(f, -1, 0)                  # (y+1, x)
    fcc = np.roll(fyp, -1, 1)                # (y+1, x+1)
    s = (np.angle(np.conj(f) * fxp)         # A -> B
         + np.angle(np.conj(fxp) * fcc)     # B -> C
         + np.angle(np.conj(fcc) * fyp)     # C -> D
         + np.angle(np.conj(fyp) * f))      # D -> A
    winding = np.trunc(s / 6.0).astype(int)  # f90: dint(.../6)
    if rho_min > 0.0:
        rho = a * a + b * b
        rmax = np.maximum.reduce([rho, np.roll(rho, -1, 1),
                                  np.roll(np.roll(rho, -1, 1), -1, 0),
                                  np.roll(rho, -1, 0)])
        winding = np.where(rmax > rho_min, winding, 0)
    n_plus = int(np.sum(winding >= 1))
    n_minus = int(np.sum(winding <= -1))
    n_total = n_plus + n_minus
    net = int(np.sum(winding))
    return winding, (n_plus, n_minus, n_total, net)
