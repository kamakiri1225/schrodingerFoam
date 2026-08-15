# 00_2_harmonicOscillator_1D — 1次元 調和振動子ポテンシャル

線形シュレーディンガー方程式（$g=0$）を **1次元**で解き、調和トラップ $V=\tfrac12\omega^2x^2$ の中で
**コヒーレント状態**（変位させた基底状態ガウス）が形を崩さずに振動する様子を再現するデモ。

## 物理

$$
i\,\frac{\partial\psi}{\partial t} = -D\,\partial_x^2\psi + \tfrac12\omega^2x^2\,\psi,
\qquad D=\tfrac12,\ g=0,\ \omega=1
$$

- 基底状態：$\psi_0(x)\propto e^{-\omega x^2/2}$、エネルギー $E_0=\omega/2$、$|\psi_0|^2$ の幅 $\sigma=1/\sqrt{2\omega}$。
- 初期条件：基底状態を $x_0=4$ だけ**ずらした**コヒーレント状態 $\psi(x,0)=A\,e^{-\omega(x-x_0)^2/2}$（$A=(\omega/\pi)^{1/4}$）、$\psi_\mathrm{im}=0$（初速ゼロ）。
- コヒーレント状態は**波束が広がらず**、中心が古典運動 $x(t)=x_0\cos\omega t$ に従って振動する。周期 $T=2\pi/\omega\approx6.28$。

## 数値設定

- 領域 $x\in[-12,12]$、$n_x=240$（$\Delta x=0.1$）。
- y は**1セルのまま幅 3 に拡幅**（ParaView で帯状に見やすくするため。物理は 1 次元のまま）、z は 1 セル厚。y・z とも `empty`。
- 境界 x 両端：`fixedValue 0`（古典的禁止領域の遠方なので $\psi\approx0$）。
- realTime（Crank–Nicolson＋`nCorrectors 4`）、$\Delta t=0.005$、$t\le12.566$（2 周期）。
- ノルム保存を確認（`norm = 0.3` 一定＝断面積 $3\times0.1$ × 規格化 $\int|\psi|^2dx=1$）。
- ピーク密度 $\approx(\omega/\pi)^{1/2}\approx0.564$ が解析値と一致。

## 実行

```bash
openfoam2512 -c 'cd run/00_2_harmonicOscillator_1D && blockMesh && schrodingerFoam'
openfoam2512 -c 'cd run/00_2_harmonicOscillator_1D && foamToVTK -ascii'
python3 tools/render1d.py run/00_2_harmonicOscillator_1D/VTK figures/00_2_harmonicOscillator_1D magSqrPsi
```

## 結果

![1D harmonic oscillator](../../figures/00_2_harmonicOscillator_1D/magSqrPsi.gif)

放物線ポテンシャル（赤破線）の中で、ガウス波束が**形を保ったまま左右に往復**する。
波束が広がらないのがコヒーレント状態の特徴で、中心は古典粒子とまったく同じ単振動をする。

> 補足：`mode imaginaryTime`＋`normalize` にすれば、同じセットアップで**基底状態**（振動しない
> 静止ガウス）を作れる。虚時間発展の 1 次元検証にも使える。
