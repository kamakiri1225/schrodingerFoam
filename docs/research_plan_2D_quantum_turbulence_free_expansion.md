# 研究計画：2次元量子乱流BECの自由膨張

## 1. 研究テーマ

**有限体積法による2次元量子乱流BECの自由膨張解析  
— 渦分布とアスペクト比維持条件の検討 —**

英題候補：

> Free Expansion of Two-Dimensional Turbulent Bose–Einstein Condensates
> Using a Finite-Volume Gross–Pitaevskii Solver

### 中心課題

`schrodingerFoam` を用いて2次元量子乱流を動的に生成し、解放前の乱流構造と、
異方的トラップ解放後の膨張異方性・自己相似性との関係を明らかにする。

特に、次の問いを主対象とする。

1. 2次元量子乱流でもアスペクト比を保つ自己相似膨張は起こるか。
2. アスペクト比維持を決めるのは、総渦数、渦配置、偏極度、渦双極子率のどれか。
3. 非圧縮性渦運動と圧縮性の音波・密度揺らぎは、膨張へどう寄与するか。
4. 3次元乱流BECで観測された自己相似膨張と同じ機構か。

## 2. 現在地点

### 完了していること

- `schrodingerFoam` は
  
  $$
  i\partial_t\psi=
  \left[-D\nabla^2+V(\boldsymbol{x},t)+g|\psi|^2\right]\psi
  $$
  
  の実時間発展と、初期状態を作る虚時間発展に対応している。
- `releaseTime` により、指定時刻以降に外部トラップを除去できる。
- `run/00_3_release2D` では、異方的なガウス基底状態を解放し、縦横比反転を確認した。
- 同ケースの最終時刻は $t=10$。初期幅 $(4,0.5)$ に対し、OpenFOAM結果は
  $(\sigma_x,\sigma_y)=(4.123,7.786)$ となった。
- 古典比較用として「同じ初期位置密度＋等方速度分散」の解析モデルも実装した。
- 単一渦、渦対、ダークソリトン崩壊による渦生成の既存ケースがある。

### 重要な未完了点

`00_3_release2D` は $g=0$ の**線形シュレーディンガー波束**であり、相互作用するBECの
自由膨張検証ではない。研究の基準ケースには、$g>0$ の異方的トラップ基底状態を
虚時間発展で作り、その状態を実時間計算へ引き渡す必要がある。

また、$t=10$ では現在のy境界へ波束の裾が達して弱い反射縞が出る。乱流状態は通常さらに
高速・広範囲に膨張するため、領域拡大、動的再配置、または吸収境界の検討が必要である。

## 3. 物理的背景

異方的トラップ内の通常のBECを解放すると、強く閉じ込められていた短軸方向へ速く膨張する。
半径を

$$
R_j(t)=\sqrt{\langle x_j^2\rangle-\langle x_j\rangle^2}
$$

とすれば、アスペクト比は

$$
A(t)=\frac{R_x(t)}{R_y(t)}
$$

である。渦なしBECでは一般に $A(t)$ が初期値を横切り、縦横比が反転する。

主な駆動要因は次のとおり。

- 狭い方向ほど運動量幅が大きいという量子的な運動エネルギー
- 相互作用BECで異方的密度分布に蓄えられた相互作用エネルギー

一方、単一温度の古典熱平衡集団では速度分散が等方的なので、十分長い自由飛行後には
円形へ近づく。ただし、現在の比較GIF中央は完全な正準熱平衡ではなく、初期位置形状を
量子側と揃えて等方速度だけを与えた対照モデルである。

## 4. 主要な先行研究

### 4.1 3次元乱流BECの異常膨張

Hennらは振動磁場で $^{87}\mathrm{Rb}$ BECを励起し、乱雑な渦糸を伴う状態を生成した。
通常のBECで見られるアスペクト比反転が乱流状態では抑制され、形状比をほぼ保つ膨張を
報告している。実験論文であり、異常膨張の完全な機構を直接数値的に確定したものではない。

- E. A. L. Henn et al.,
  [Emergence of Turbulence in an Oscillating Bose-Einstein Condensate](https://arxiv.org/abs/0904.2564)（arXiv・無料閲覧）,
  *Phys. Rev. Lett.* **103**, 045301 (2009).

Caracanhasらは、凝縮体内に分布する渦度が通常の異方的膨張と競合するという流体モデルで、
自己相似的な異常膨張を定性的に説明した。これは個々の渦芯を解像するGP直接計算ではない。

- M. Caracanhas et al.,
  [Self-similar Expansion of the Density Profile in a Turbulent Bose-Einstein Condensate](https://arxiv.org/abs/1103.2039)（arXiv・無料閲覧）,
  *J. Low Temp. Phys.* **166**, 49–58 (2012).

近年の実験・3次元GP計算では、振動駆動されたBEC乱流が一様な渦糸タングルではなく、
強い密度揺らぎ、小さな渦リング、非一様な渦分布、圧縮性波動の混合状態であることが
示されている。

- H. A. J. Middleton-Spencer et al.,
  [Strong Quantum Turbulence in Bose-Einstein Condensates](https://arxiv.org/abs/2204.08544)（arXiv・無料閲覧）,
  *Phys. Rev. Research* **5**, 043081 (2023).
- H. A. J. Middleton-Spencer,
  [On the Expansion of Turbulent Bose-Einstein Condensates](https://theses.ncl.ac.uk/jspui/handle/10443/6318),
  PhD thesis, Newcastle University (2024).

### 4.2 2次元GP方程式による渦入り自由膨張

Tsuchitani・Tsubotaは2次元GP方程式を用い、渦なし、単一渦、渦対1組、渦対2組を含む
凝縮体の自由膨張を比較した。渦対の異方的速度場が通常の異方的膨張と競合し、条件により
アスペクト比反転を抑制し得ることを示した。

- R. Tsuchitani and M. Tsubota,
  [Expansion of a Bose-Einstein Condensate with Vortices](https://arxiv.org/abs/1312.0093)（arXiv・無料閲覧）,
  *J. Low Temp. Phys.* **174**, 223–231 (2014).

ただし、この研究は最大でも渦対2組の単純モデルであり、乱流を動的生成した計算ではない。
エネルギーフラックス、渦クラスタリング、乱流スペクトルによる分類も主題ではない。

### 4.3 2次元量子乱流の生成・減衰

Staggらは障害物後流から乱雑な渦・反渦分布を作り、障害物除去後の渦数減衰を2次元GP計算で
調べた。渦・反渦消滅、境界への流出、散逸の影響を扱うが、トラップ全体を解放した後の
アスペクト比は主題ではない。

- G. W. Stagg et al.,
  [Generation and Decay of Two-Dimensional Quantum Turbulence in a Trapped Bose-Einstein Condensate](https://arxiv.org/abs/1408.3268)（arXiv・無料閲覧）,
  *Phys. Rev. A* **91**, 013612 (2015).

## 5. 暫定的な研究ギャップ

文献調査時点での整理は次のとおり。

| 課題 | 状況 |
|---|---|
| 2次元・渦なしBECの自由膨張 | 既知 |
| 2次元・単一渦、少数渦対の自由膨張 | 既報 |
| 2次元量子乱流の生成・減衰 | 多数の研究あり |
| 動的生成した2次元量子乱流を異方的トラップから解放 | 明確な研究は限定的 |
| 乱流構造とアスペクト比維持条件の系統的対応 | 研究余地が大きい |
| 渦・音波・密度揺らぎの膨張寄与の分離 | 特に研究価値が高い |

したがって、現段階で採用する慎重な表現は次のとおり。

> 2次元量子乱流の研究は多数存在するが、動的に生成した乱流状態を異方的トラップから
> 解放し、乱流構造と膨張異方性の関係を系統的に調べた研究は限定的である。

「世界初」は現時点では使用しない。投稿前にWeb of Science、Scopus、Google Scholarで、
主要論文の被引用・引用文献を双方向に追跡する。

## 6. 研究仮説

単なる渦数よりも、解放時の速度場の異方性が膨張を直接決める可能性が高い。検証対象として
次の仮説を置く。

1. **H1：** アスペクト比維持は総渦数 $N_v$ より、非圧縮性運動エネルギーの異方性と強く相関する。
2. **H2：** 渦双極子の配向分布が膨張方向を決め、等方的配向では反転抑制が弱くなる。
3. **H3：** 同符号渦クラスタと非偏極な渦・反渦気体では、同じ渦数でも膨張率が異なる。
4. **H4：** 圧縮性エネルギーは平均アスペクト比より、密度分布の非自己相似性や波紋を強く増加させる。
5. **H5：** アスペクト比一定は自己相似膨張の必要条件だが十分条件ではない。

## 7. 「量子乱流」と呼ぶための解放前診断

### 7.1 渦統計

- 正渦数 $N_+$、負渦数 $N_-$、総渦数 $N_v=N_++N_-$
- 偏極度
  
  $$
  P=\frac{N_+-N_-}{N_++N_-}
  $$
- 最近接渦間距離と符号相関
- 渦・反渦双極子率
- 同符号渦クラスタ率とクラスタサイズ分布
- 渦の空間分布、中心からの半径分布

渦検出は、セル面またはセル中心位相から閉ループ位相差を計算し、$\pm2\pi$ 巻きを数える。
位相差は必ず $(-\pi,\pi]$ にラップする。

### 7.2 エネルギー分解

$$
\psi=\sqrt{\rho}\,e^{i\theta},
\qquad
\boldsymbol{v}=\frac{\hbar}{m}\nabla\theta
$$

とし、渦芯での速度特異性を避けるため密度重み付き速度

$$
\boldsymbol{u}=\sqrt{\rho}\,\boldsymbol{v}
$$

をHelmholtz分解する。

$$
\boldsymbol{u}=\boldsymbol{u}_{\mathrm{inc}}+
\boldsymbol{u}_{\mathrm{comp}},
\qquad
\nabla\cdot\boldsymbol{u}_{\mathrm{inc}}=0,
\qquad
\nabla\times\boldsymbol{u}_{\mathrm{comp}}=0
$$

評価量：

- 非圧縮性運動エネルギー $E_{\mathrm{inc}}$
- 圧縮性運動エネルギー $E_{\mathrm{comp}}$
- 量子圧力、相互作用エネルギー、トラップエネルギー
- $E_{\mathrm{inc}}(k)$、$E_{\mathrm{comp}}(k)$
- エネルギーフラックス
- 圧縮性比 $E_{\mathrm{comp}}/(E_{\mathrm{inc}}+E_{\mathrm{comp}})$
- 非圧縮性エネルギーテンソルと異方度
  
  $$
  K_{ij}^{\mathrm{inc}}
  =\frac12\int u_{\mathrm{inc},i}u_{\mathrm{inc},j}\,dA,
  \qquad
  Q_K=\frac{K_{xx}^{\mathrm{inc}}-K_{yy}^{\mathrm{inc}}}
  {K_{xx}^{\mathrm{inc}}+K_{yy}^{\mathrm{inc}}}
  $$

有限トラップ上のFFT・Helmholtz分解では、非周期境界によるスペクトル漏れを検証する。
窓関数、ゼロパディング、周期的大領域のどれを採用するかを明記する。

## 8. 解放後の評価量

### 8.1 半径、アスペクト比、膨張率

$$
R_x=\sqrt{\langle x^2\rangle-\langle x\rangle^2},
\qquad
R_y=\sqrt{\langle y^2\rangle-\langle y\rangle^2}
$$

$$
A(t)=\frac{R_x(t)}{R_y(t)},
\qquad
b_x(t)=\frac{R_x(t)}{R_x(t_r)},
\qquad
b_y(t)=\frac{R_y(t)}{R_y(t_r)}
$$

ここで $t_r$ はトラップ解放時刻である。

### 8.2 自己相似性

アスペクト比一定だけでは自己相似性を保証しない。規格化座標

$$
X=\frac{x-\langle x\rangle}{R_x(t)},
\qquad
Y=\frac{y-\langle y\rangle}{R_y(t)}
$$

上で密度を再配置し、解放時の規格化密度との差を例えば

$$
\epsilon_{\mathrm{ss}}(t)=
\frac{\|\widetilde\rho(X,Y,t)-\widetilde\rho(X,Y,t_r)\|_2}
{\|\widetilde\rho(X,Y,t_r)\|_2}
$$

で定量化する。内部の渦芯・音波構造を含む指標と、平滑化した包絡密度の指標を分けて出す。

### 8.3 膨張中の乱流変化

- $N_+(t)$、$N_-(t)$、渦消滅率
- クラスタ率・双極子率の時間変化
- $E_{\mathrm{inc}}/E_{\mathrm{comp}}$ の時間変化
- 渦芯間隔とヒーリング長の比
- 境界へ流出した粒子数・渦数

## 9. 計算ケース計画

### Stage 1：相互作用BEC基準計算

1. $g>0$ の異方的トラップ基底状態を虚時間発展で生成する。
2. 生成した複素波動関数を実時間ケースへ引き渡す。
3. 渦なし状態を解放し、GPスケーリング解または高精度参照解と比較する。
4. 格子、時間刻み、領域幅、境界条件への依存性を調べる。
5. 解放後の粒子数・エネルギー保存を確認する。

合格条件の案：

- 粒子数相対誤差 $<10^{-4}$
- $R_x,R_y$ の参照解との差 $<2\%$
- 境界到達前の結果のみを定量比較に使用
- 空間・時間解像度を上げたときのアスペクト比差 $<2\%$

### Stage 2：既存研究の再現

1. 渦なし
2. 単一渦
3. 配向を変えた渦対1組
4. 配向を変えた渦対2組

Tsuchitani・Tsubotaの定性的傾向と $A(t)$ を比較し、有限体積法の再現性を確認する。

### Stage 3：動的な2次元量子乱流生成

候補手法：

- 移動障害物による後流
- 楕円障害物の往復運動
- 時間変動する異方的トラップ
- 多重荷電渦の崩壊
- ダークソリトンのsnake instability

「乱流」の採用判定は、渦数だけでなく、渦統計、エネルギー分解、スペクトルを組み合わせる。

### Stage 4：パラメータマップ

変化させる量：

- 初期トラップ比 $\omega_x/\omega_y$
- 相互作用強度 $gN$
- 総渦数 $N_v$
- 偏極度 $P$
- 渦クラスタ率
- 渦双極子率と配向
- $E_{\mathrm{comp}}/E_{\mathrm{inc}}$
- 解放前の乱流発達時間

主要な応答量：

- $A(t)$ と長時間平均
- $b_x(t),b_y(t)$
- 自己相似誤差 $\epsilon_{\mathrm{ss}}$
- 境界到達時刻
- 渦数・エネルギーの減衰率

## 10. 推奨ケース構成

実装時の仮ディレクトリ名：

```text
run/
  05_0_release2D_interacting_noVortex/
  05_1_release2D_singleVortex/
  05_2_release2D_oneDipole/
  05_3_release2D_twoDipoles/
  05_4_turbulence2D_generation/
  05_5_turbulence2D_release/
tools/
  detect_vortices.py
  decompose_kinetic_energy.py
  measure_expansion.py
  measure_self_similarity.py
```

乱流生成と自由膨張を別ケースに分け、解放直前の場を明示的な初期条件として保存する。
これにより同じ乱流スナップショットを異なるトラップ解放条件で再利用できる。

## 11. 数値上の重点課題

1. **領域サイズ**：自由膨張中の境界反射を排除する。
2. **渦芯解像度**：ヒーリング長あたり最低限必要なセル数を格子収束で決める。
3. **初期状態引き渡し**：虚時間と実時間で規格化・化学ポテンシャルを一貫させる。
4. **エネルギー収支**：トラップ瞬時除去によるエネルギー変化と、その後の保存を分離する。
5. **スペクトル評価**：有限領域、非一様密度、非周期境界の影響を定量化する。
6. **散逸の扱い**：まず無散逸GPを基準とし、必要な場合のみ現象論的散逸を追加する。
7. **再現性**：乱流初期条件ごとに乱数seedを保存し、複数seedの統計を取る。

## 12. 最小実行計画

直近の順序は次とする。

1. `00_3_release2D` の領域を拡大し、$t=10$ の境界反射を除いた基準結果を作る。
2. $g>0$ の異方的基底状態を虚時間発展で作る新ケースを追加する。
3. 渦なし相互作用BECの自由膨張とアスペクト比反転を検証する。
4. 単一渦、渦対1組、渦対2組を順に再現する。
5. 渦検出と $R_x,R_y,A,b_x,b_y$ の自動計測を先に実装する。
6. その後、動的乱流生成とエネルギー分解へ進む。

最初の研究成果としては、**有限体積GPソルバーの検証＋少数渦研究の再現**を確実に行い、
次に**動的生成した2次元乱流の自由膨張マップ**を本命成果とする。

## 13. 新規性の位置づけ

| 内容 | 位置づけ |
|---|---|
| 渦なしBECの膨張 | 検証用、既知 |
| 単一渦・少数渦対の膨張 | 既報再現 |
| 多数渦を人為配置した膨張 | 中間段階 |
| 動的生成した2次元量子乱流の膨張 | 新規性が期待できる |
| 乱流タイプ別の膨張比較 | 研究価値が高い |
| アスペクト比維持条件のマップ化 | 研究価値が高い |
| 渦・音波・密度揺らぎの寄与分離 | 特に研究価値が高い |
| OpenFOAM有限体積法による実装・公開 | 数値手法上の独自性 |

## 14. 投稿前の文献調査

以下の検索語を組み合わせ、引用追跡結果を別表として残す。

```text
2D quantum turbulence free expansion
two-dimensional turbulent BEC time of flight
aspect ratio turbulent condensate expansion
Gross-Pitaevskii vortex dipole expansion
self-similar expansion quantum turbulence 2D
vortex clustering expansion Bose-Einstein condensate
```

データベース：Google Scholar、Web of Science、Scopus、INSPIRE、arXiv。

調査終了までは「世界初」ではなく、**systematic studies appear limited** または
**to the best of our knowledge** と表現する。
