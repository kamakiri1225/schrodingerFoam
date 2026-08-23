<!-- 第7回：ダークソリトンの崩壊から2次元量子乱流へ -->

> **シリーズ「OpenFOAM でシュレーディンガー方程式を解く」応用編・第7回。**
> ①ソルバ自作／②量子トンネル効果／③調和振動子ポテンシャル／④量子 vs 古典・波束の広がり／⑤虚時間発展で基底状態を作る／⑥走る渦対／**⑦2次元量子乱流（本記事）**／総集編（第0回）：全体ダイジェスト
> リポジトリ：<https://github.com/kamakiri1225/schrodingerFoam>／ケース：`run/06_solitonLattice8x8_whiteNoise_muShift`, `run/07_fourCondensate_trap_muShift`

第6回の**走る渦対**の先へ進みます。シリーズの目標だった「**ダークソリトンが波打って（snake instability）多数の量子渦に分裂し、2次元量子乱流になる**」過程を計算し、そのエネルギースペクトル・カスケードの向き・渦数の減衰までを定量化します。参考は Numasato–Tsubota–L'vov, *Direct Energy Cascade in 2D Compressible Quantum Turbulence*, PRE **81**, 016303 (2010) と、小林の学位論文（`docs/main.pdf`）です。

<p align="center">
  <img src="../figures/06_solitonLattice8x8/density_phase.gif" width="640"><br>
  <em>ソリトン格子が崩壊して量子渦の乱流になる（左＝密度、右＝位相）</em>
</p>

## 量子乱流とは — 古典乱流との決定的な違い

超流動（BEC）の流れは波動関数 $\psi=\sqrt{\rho}\,e^{i\theta}$ で表され、速度は位相の勾配 $\mathbf v = 2D\,\nabla\theta$ で決まります。渦は好き勝手な強さを取れず、**位相が渦を1周すると $2\pi$ の整数倍**という量子化条件を満たします（第6回参照）。渦の中心では密度が 0 に落ちます。

古典的な2次元乱流では、エネルギーは**大きなスケールへ**流れます（逆カスケード）。これは「エンストロフィー（渦度の2乗）」が保存するためです。ところが**量子乱流ではエンストロフィーが保存しません**——渦は量子化された点で、渦・反渦が出会えば対消滅して消えるからです。この違いのおかげで、量子乱流では古典と**逆に、エネルギーが小さなスケールへ**流れます（**順カスケード direct**）。本記事の主役はこの現象です。

## 手法A：ソリトン格子から乱流を作る（ケース06）

一様なバルクに、1次元ダークソリトン解の積 $\Psi=\prod \tanh(x-x_i)\tanh(y-y_i)$ で**格子状のダークソリトン**を置きます（$8\times8$、間隔6）。まず虚時間発展（第5回）で背景を整え、微小な白色ノイズを乗せて実時間発展させると、まっすぐな溝（ソリトン）が波打って**量子渦へ分裂**し、乱流になります。

### エネルギースペクトルに $k^{-5/3}$ が現れる

速度場をフーリエ変換して運動エネルギーを波数 $k$ ごとに分けたのが**エネルギースペクトル** $E(k)$ です。乱流が発達すると、慣性領域に Kolmogorov–Obukhov の **$k^{-5/3}$** 則が現れます。

<p align="center">
  <img src="../figures/06_solitonLattice8x8/spectrum_evolution.gif" width="600"><br>
  <em>スペクトルの時間変化。青＝非圧縮（渦）、赤＝圧縮（音波）、黒点線＝ k⁻⁵ᐟ³（固定係数）</em>
</p>

$k^{-5/3}$ 線の**係数は物理で決まります**。Kolmogorov 則 $E^i(k)=C\,\varepsilon^{2/3}k^{-5/3}$ の $\varepsilon$ はエネルギー流束で、後述のカスケード解析から測ります（点線は毎フレーム同じ位置に固定）。

### 速度を「渦」と「音」に分ける — Helmholtz 分解

密度重み速度 $\mathbf w=\sqrt{\rho}\,\mathbf v$ を、**非圧縮成分 $\mathbf w^i$（湧き出しなし＝回転＝量子渦）** と **圧縮成分 $\mathbf w^c$（回転なし＝音波）** に分解します（Helmholtz 分解）。

<p align="center">
  <img src="../figures/06_solitonLattice8x8/decomposition.gif" width="720"><br>
  <em>左＝|wⁱ| 非圧縮（渦芯に局在、±点が量子渦）／右＝|wᶜ| 圧縮（音波が全域に放射）</em>
</p>

非圧縮エネルギーは渦芯に局在し、圧縮エネルギーは音波として全域に広がります。渦が対消滅するたびに、その運動エネルギーが**音波として放出**されている——これがまさに「小スケールで非圧縮→圧縮に変換」という量子乱流の描像です。

### エネルギーカスケードの向き

波数 $k$ を通って小スケール側へ流れる非圧縮エネルギーの流束を、論文と同じ式で測ります：

$$
\varepsilon^i(k,t) = -\int_{k_m}^{k}\frac{\partial E^i_{kin}(k',t)}{\partial t}\,dk'
$$

慣性帯で $\varepsilon^i>0$ なら**順カスケード（大→小）**、$<0$ なら逆です。

<p align="center">
  <img src="../figures/06_solitonLattice8x8/cascade_direction.png" width="560"><br>
  <em>発達乱流域（t=20–100）の平均流束。慣性帯で ε≈+3.2 &gt; 0 ＝ 順カスケード direct</em>
</p>

**注意点**：全区間を平均すると逆に見えます。これは初期の「ソリトン→渦形成」の過渡（$t<16$）が支配するためで、**発達乱流域（$t\gtrsim20$）だけを見ると明確に順カスケード**になります（論文と一致）。時刻ごとの流束の変化は下のGIFのとおりです。

<p align="center">
  <img src="../figures/06_solitonLattice8x8/cascade_evolution.gif" width="560"><br>
  <em>ε ⁱ(k) の時間変化（縦軸 ±15 固定）。発達後は慣性帯で正＝順カスケード</em>
</p>

### 量子渦の数え方と「偽渦」の除去

渦は各格子プラケットの**位相の巻き数**（$\pm1$）で数えます（参考 Fortran コード `soliton_hist_ensamble.f90` と同一式）。ところが素朴に数えると、初期に白色ノイズが**ソリトンの節線（密度≈0）上**に格子スケールの偽の $2\pi$ 巻きを大量に作り、$t=0$ で 1452 本と出てしまいます（本当は渦0本）。

そこで**擬渦度** $\omega_{ps}=\partial_x\mathrm{Re}\,\partial_y\mathrm{Im}-\partial_y\mathrm{Re}\,\partial_x\mathrm{Im}$ を使います。本物の渦は $\psi=0$ の孤立点で $\nabla\mathrm{Re}\perp\nabla\mathrm{Im}$ となり $|\omega_{ps}|$ が大きい。一方ソリトン節線は $\mathrm{Im}\approx0$ なので $|\omega_{ps}|$ が小さい——これで綺麗に分離できます。

<p align="center">
  <img src="../figures/06_solitonLattice8x8/vortex_compare.gif" width="720"><br>
  <em>左＝生の巻き数（節線上に偽渦だらけ）／右＝擬渦度で偽渦除去。t=0 は 1452→0 本に</em>
</p>

偽渦を除いた渦数の時間変化を両対数で見ると、**「$t=0$ は渦0本 →（snake instability で）渦生成のピーク～250 → 対消滅で $N_v\propto t^{-1}$ 減衰」** という、生成から減衰までの物理が正しく見えます。

<p align="center">
  <img src="../figures/06_solitonLattice8x8/vortex_decay.png" width="560"><br>
  <em>灰＝生の巻き数／青＝偽渦除去後。後期は t⁻¹ 則（渦・反渦の対消滅）</em>
</p>

## 手法B：4つの凝縮体の干渉から乱流を作る（ケース07）

より実験に近い作り方も試します（学位論文 付録A）。調和トラップ $\frac12\omega^2(x^2+y^2)$ に**斥力の十字ポテンシャル** $V_0(e^{-x^2/d^2}+e^{-y^2/d^2})$ を足して、凝縮体を**4つに分離**した基底状態を作ります。

<p align="center">
  <img src="../figures/07_fourCondensate/ground_state.png" width="420"><br>
  <em>虚時間で作った4分離の基底状態（論文 図A.2 に対応）</em>
</p>

次に**斥力を切る**と、4つの凝縮体が中心へ落ちて**干渉**します。干渉縞（密度の低い線）が正方格子状のダークソリトンになり（位相は0/π交互）、それが崩壊して量子渦になります。

<p align="center">
  <img src="../figures/07_fourCondensate/density_phase.gif" width="640"><br>
  <em>4凝縮体が中心衝突→干渉→十字ダークソリトン→量子渦（左＝密度、右＝位相）</em>
</p>

<p align="center">
  <img src="../figures/07_fourCondensate/vortex_distribution.gif" width="480"><br>
  <em>干渉格子の暗い芯に乗る量子渦（t≈18 でピーク16本）</em>
</p>

こちらは**弱い・非定常な乱流**でした。渦数はトラップ振動（周期 $2\pi/\omega\approx78$）に伴って**バースト的に増減**し、$t\approx105$ でほぼ0へ。カスケードも $\varepsilon^i\approx+0.15$ と、06 の $+3.2$ より桁小さい弱い direct です。またノイズをソリトン線に乗せないので、**偽渦の問題は元々ほとんど起きません**（生の数と擬渦度除去後が一致）。

これは論文自身の指摘（`main.pdf` p.22）とも一致します——「トラップ系は渦コアに対し凝縮体が充分大きくないと統計が難しい。本研究では**一様系のソリトン格子で乱流統計をとった**」。つまり **06（一様）＝統計向き、07（干渉）＝物理的な実証**、という役割分担です。

## まとめ：06 と 07 の対比

| | 06 ソリトン格子（一様） | 07 干渉（トラップ） |
|---|---|---|
| 渦数 | 数百 → 綺麗な $t^{-1}$ 減衰 | 16ピーク → バースト → t≈105で0 |
| カスケード | 明確な direct（$\varepsilon\approx+3.2$） | 弱い direct（$\varepsilon\approx+0.15$） |
| 乱流 | 豊富・準定常（統計向き） | 弱い・非定常（小系＋振動） |
| 偽渦 | ノイズ由来多数 → 擬渦度で除去 | 元々ほぼ無し |
| 位置づけ | **統計** | **物理的実証** |

古典2次元乱流とは逆に、**量子乱流では非圧縮エネルギーが小スケールへ順カスケードし、渦の対消滅を通じて音波に変換される**——この描像を、スペクトル・流束・Helmholtz 分解・渦数減衰の4点セットで確認できました。

## 基本スケールと崩壊波数（手計算）

- 癒し長 $\xi=\sqrt{D/(g n_0)}$（$n_0=1$ で $\approx0.707$）、ダークソリトン幅 $\sqrt2\,\xi\approx1$。
- ダークソリトンは2次元で**横方向のゆらぎに不安定**（snake instability）。不安定帯は $0<q<q_c$、$q_c\xi=O(1)$、最速成長は $q_*\xi\approx0.5$。→ 崩壊してできる渦の間隔の目安 $\lambda_*=2\pi/q_*\approx4\pi\xi\approx9$。
- **エネルギーの注入スケール**：本計算は強制なしの減衰乱流で、初期状態が唯一のエネルギー源。06 ではソリトンスケール $k\sim1$ に注入され、そこから小スケール（音波）へ順カスケードします。

## 実行手順

```bash
# 06: ソリトン格子（虚時間で整える → ノイズ+実時間）
openfoam2512 -c 'cd run/06_solitonLattice8x8_whiteNoise_muShift/imaginaryTime && blockMesh && schrodingerFoam'
openfoam2512 -c 'cd run/06_solitonLattice8x8_whiteNoise_muShift/realTime && schrodingerFoam && foamToVTK -fields "(Psire Psiim)"'

# 解析（tools/）
python3 tools/qt_spectrum_gif.py  <VTK> figures/06.../spectrum_evolution.gif
python3 tools/qt_flux.py          <VTK> figures/06.../cascade_direction.png
python3 tools/qt_velocity_gif.py  <VTK> figures/06.../decomposition.gif --pseudo-frac 0.03
python3 tools/qt_vortex_decay.py  <VTK> figures/06.../vortex_decay.png   --pseudo-frac 0.03
python3 tools/qt_vortex_gif.py    <VTK> figures/06.../vortex_compare.gif --compare --pseudo-frac 0.03
# 07（トラップ系）は上記に --rho-frac 0.3（渦は --rho-frac 0.05）を付ける
```

解析ツールの詳細（各スクリプトの役割・共有モジュール `qt_analysis.py`・Helmholtz 分解や流束の式）は、付録としてこの下にまとめてあります。

---
