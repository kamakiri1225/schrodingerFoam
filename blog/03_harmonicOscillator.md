<!-- 第3回：1次元 調和振動子ポテンシャル（コヒーレント状態） -->

> **シリーズ「OpenFOAM でシュレーディンガー方程式を解く」全7回・第3回。**
> ①ソルバ自作／②量子トンネル効果／**③調和振動子ポテンシャル（本記事）**／④量子 vs 古典・波束の広がり／⑤虚時間発展で基底状態を作る／⑥走る渦対／⑦ダークソリトンの崩壊（ノイズ有無）
> リポジトリ：<https://github.com/kamakiri1225/schrodingerFoam>／ケース：`run/00_2_harmonicOscillator_1D`

第2回に続き 1 次元・線形（$g=0$）のシュレーディンガー方程式を解きます。今回は**調和振動子ポテンシャル**（放物線トラップ）$V=\tfrac12\omega^2x^2$ の中で、**コヒーレント状態**（変位させた基底状態ガウス）が**形を崩さずに単振動する**様子を再現します。

<p align="center">
  <img src="../figures/00_2_harmonicOscillator_1D/magSqrPsi.gif" width="640"><br>
  <em>放物線ポテンシャル（赤破線）の中で、ガウス波束が形を保ったまま左右に往復する</em>
</p>

## コヒーレント状態とは

量子調和振動子の**基底状態**は $\psi_0(x)\propto e^{-\omega x^2/2}$ というガウス分布で、エネルギーは $E_0=\omega/2$（ゼロ点エネルギー）です。この基底状態を丸ごと $x_0$ だけ横にずらしたものが**コヒーレント状態**です。

コヒーレント状態の面白い性質は、**時間が経っても波束が広がらず、その中心が古典粒子とまったく同じ単振動をする**こと。

$$
x(t)=x_0\cos\omega t,\qquad T=\frac{2\pi}{\omega}
$$

自由空間ならガウス波束は時間とともに必ず広がります（第4回で詳しく扱います）が、調和トラップの中では「広がろうとする効果」と「トラップが引き戻す効果」がちょうど釣り合い、**幅を保ったまま**振動します。これは「量子力学が古典力学に最も近づく状態」として知られています。

## 解く方程式と設定

$$
i\,\frac{\partial\psi}{\partial t} = -D\,\partial_x^2\psi + \tfrac12\omega^2x^2\,\psi,
\qquad D=\tfrac12,\ g=0,\ \omega=1
$$

- **基底状態**：$\psi_0(x)\propto e^{-\omega x^2/2}$、$|\psi_0|^2$ の幅 $\sigma=1/\sqrt{2\omega}$。
- **初期条件**：基底状態を $x_0=4$ だけ**ずらした**コヒーレント状態 $\psi(x,0)=A\,e^{-\omega(x-x_0)^2/2}$（$A=(\omega/\pi)^{1/4}$）、$\psi_\mathrm{im}=0$（初速ゼロ）。

初期波束（`0/Psire`）とトラップ（`0/Vext`）を `#codeStream` で与えます。

```cpp
// 0/Psire : 変位させた基底状態ガウス
const scalar omega = 1.0, x0 = 4.0;
const scalar A = Foam::pow(omega/M_PI, 0.25);
psi[i] = A*Foam::exp(-0.5*omega*Foam::sqr(x - x0));   // Psiim は 0（初速ゼロ）

// 0/Vext : 調和トラップ  V = 0.5 * omega^2 * x^2
V[i] = 0.5*Foam::sqr(omega)*Foam::sqr(x);
```

## 数値設定

- 領域 $x\in[-12,12]$、$n_x=240$（$\Delta x=0.1$）。
- y は 1 セルのまま幅 3 に拡幅（表示用。物理は 1 次元）、z は 1 セル厚。y・z とも `empty`。
- 境界 x 両端：`fixedValue 0`（古典的禁止領域の遠方なので $\psi\approx0$）。
- realTime（Crank–Nicolson＋`nCorrectors 4`）、$\Delta t=0.005$、$t\le12.566$（2 周期）。

## 結果：広がらずに単振動

放物線ポテンシャルの中で、ガウス波束が**形を保ったまま**左右に往復しました。中心は古典粒子とまったく同じ単振動 $x(t)=4\cos t$ を描き、周期 $T=2\pi\approx6.28$ で戻ってきます。

定量的にも、ノルムは計算中ずっと一定（$=0.3$＝断面積 $3\times0.1$ ×規格化 $1$）、ピーク密度は解析値 $(\omega/\pi)^{1/2}\approx0.564$ と一致しました。波束が広がらないのがコヒーレント状態の証拠です。

> **基底状態も同じセットアップで作れる**：`mode` を `imaginaryTime` にして `normalize` を有効にすると、同じトラップで**振動しない静止ガウス（基底状態）**が得られます。これは第5回で扱う虚時間発展の 1 次元検証にもなります。

## 実行手順

```bash
openfoam2512 -c 'cd run/00_2_harmonicOscillator_1D && blockMesh && schrodingerFoam'
openfoam2512 -c 'cd run/00_2_harmonicOscillator_1D && foamToVTK -ascii'
python3 tools/render1d.py run/00_2_harmonicOscillator_1D/VTK figures/00_2_harmonicOscillator_1D magSqrPsi
```

---

次回（第4回）は、トラップを**外した**とき波束がどう広がるかを 2 次元で計算し、**量子力学と古典力学の「広がり方の違い」**を比較します。今回の「広がらない」状態とは対照的に、自由空間では量子の波束が特徴的な膨張を見せます。
