# 00_1_tunneling_1D — 1次元 量子トンネル効果

線形シュレーディンガー方程式（$g=0$）を **1次元**で解き、運動量を持ったガウス波束が
矩形ポテンシャル障壁に入射して**トンネルする**様子を再現するデモ。

## 物理

$$
i\,\frac{\partial\psi}{\partial t} = -D\,\partial_x^2\psi + V_\mathrm{ext}(x)\,\psi,
\qquad D=\tfrac12,\ g=0
$$

- 初期波束：$\psi(x,0)=A\,e^{-(x-x_0)^2/4\sigma^2}\,e^{ik_0x}$（$x_0=-15,\ \sigma=2,\ k_0=1.5$）
- 群速度 $v=k_0=1.5$、平均エネルギー $E=k_0^2/2=1.125$（波束のエネルギー幅の上端でも $E_\max\approx1.5$）
- 障壁：$V_\mathrm{ext}=V_0=3.0$（$|x|<0.6$、幅 $L=1.2$）。$V_0$ を波束の最大エネルギーより十分高くとり、全成分が古典的に透過不可な**純粋トンネル**にしている。

**透過率（実測 vs 理論）**：矩形障壁の平面波厳密解
$T=[1+V_0^2\sinh^2(\kappa L)/(4E(V_0-E))]^{-1}$、$\kappa=\sqrt{2(V_0-E)}\approx1.94$ ⇒ $T_\text{exact}\approx0.035$。
シミュレーション最終時刻 $t=22$ の密度積分では **$T\approx0.042$（透過4.2%）／$R\approx0.958$（反射95.8%）** で理論とよく一致。
（$e^{-2\kappa L}$ の粗い近似は薄い障壁では過小評価になるので使わない。）

> **障壁を無限大にしたいとき**：巨大な $V_0$ を入れるとポテンシャル項が硬直（stiff）化して
> $\Delta t$ を極端に小さくしないと破綻する。無限の壁はその位置で $\psi=0$ に固定する
> **ディリクレ境界（`fixedValue 0`）** で表すのが正しい（＝完全反射）。x 両端の `fixedValue 0`
> はまさに無限の壁として働いている。

## 数値設定

- 領域 $x\in[-40,40]$、$n_x=400$（$\Delta x=0.2$）。
- y は**1セルのまま幅 8 に拡幅**（ParaView で帯状に見やすくするため。物理は 1 次元のまま）、z は 1 セル厚。y・z とも `empty`。
- 境界 x 両端：`fixedValue 0`（遠方の壁。計算時間内に波束が到達しない広さ）。
- realTime（Crank–Nicolson＋Picard 反復 `nCorrectors 4`）、$\Delta t=0.01$、$t\le22$。
- ノルム保存を確認（`norm = 1.6` 一定＝断面積 $8\times0.2$ × 規格化 $\int|\psi|^2dx=1$）。

## 実行

```bash
openfoam2512 -c 'cd run/00_1_tunneling_1D && blockMesh && schrodingerFoam'
openfoam2512 -c 'cd run/00_1_tunneling_1D && foamToVTK -ascii'
python3 tools/render1d.py run/00_1_tunneling_1D/VTK figures/00_1_tunneling_1D magSqrPsi
```

## 結果

![1D tunneling](../../figures/00_1_tunneling_1D/magSqrPsi.gif)

波束が障壁（赤破線）に当たると、大部分が反射して入射波と干渉（左側の縞）し、
一部が障壁を**すり抜けて右側へ透過**する。これが量子トンネル効果。
