<!-- 第7回：ダークソリトンの崩壊 — 種の有無で運命が変わる -->

> **シリーズ「OpenFOAM でシュレーディンガー方程式を解く」全7回・第7回（最終回）。**
> ①ソルバ自作／②量子トンネル効果／③調和振動子ポテンシャル／④量子 vs 古典・波束の広がり／⑤虚時間発展で基底状態を作る／⑥走る渦対／**⑦ダークソリトンの崩壊（本記事）**
> リポジトリ：<https://github.com/kamakiri1225/schrodingerFoam>／ケース：`run/01_darkSoliton_realTime`, `run/04_darkSoliton_noSeed`, `run/04_darkSoliton_whiteNoise`

シリーズの目標だった現象、**2 次元ダークソリトンの崩壊**を扱います。まっすぐな密度の溝（ダークソリトン）が横方向に波打ち（snake instability）、多数の**量子渦**へ分裂する過程です。最終回では「崩壊には何が必要か」を、**種（摂動）の有無**を変えた対照実験で明らかにします。

<p align="center">
  <img src="../figures/04_darkSoliton_whiteNoise/comparison_noSeed_vs_whiteNoise.gif" width="720"><br>
  <em>左：摂動なし（崩壊しない） ／ 右：白色ノイズ種あり（渦へ崩壊）</em>
</p>

## ダークソリトンとは

ダークソリトンは、一様な凝縮体に走る**密度の溝**（$|\psi|^2$ が局所的に落ち込む面）で、GP 方程式の**厳密な定常解**です。溝の中心を横切ると波動関数の符号が反転（位相が $\pi$ 跳ぶ）します。1 次元では安定ですが、**2 次元・3 次元では横方向に不安定**で、少し波打つとその揺れが成長して溝が切れ、量子渦に分裂します。これが **snake instability（蛇行不安定性）** です。

## Case 01：決定論的な種で渦へ分裂

まず、揺れの「種」を**決定論的な cos 摂動**として与え、崩壊の基本形を見ます。

- 設定：$[-16,16]^2$・$128^2$・周期境界、$\Delta t=0.005$、$t\le60$。
- 初期は $\tanh$ の溝 2 本＋横方向 cos 摂動（第1回で示した `#codeStream`）、$\psi_\mathrm{im}=0$。

| 密度 $|\psi|^2$ | 位相 $\arg\psi$ |
|:---:|:---:|
| <img src="../figures/01_darkSoliton_realTime/density.gif" width="360"> | <img src="../figures/01_darkSoliton_realTime/phase.gif" width="360"> |

2 本のまっすぐな溝が横方向に波打ち、やがて密度ゼロの節＝**渦核**に分裂していきます。位相図では渦核のまわりで位相が $2\pi$ 巻く（渦-反渦対の芽）様子が見えます。ノルムは計算中ずっと一定で、第1回で述べた前進オイラー版の発散を克服できていることが確認できます。

## Case 04：種の有無で運命が変わる

きっかけは、参考にしたカマキリ（宇宙に入ったカマキリ）記事の Fortran 計算が**白色ノイズ**で崩壊を起こしていたこと。そこで、**まっすぐで清浄なダークソリトン**が種（seed）の有無で本当に運命を変えるのかを、同一条件（$[-16,16]^2$・$128^2$・$\Delta t=0.005$・$t\le150$）で比較しました。

- **左（`04_darkSoliton_noSeed`）**：横方向摂動も乱数も一切入れない。$t=150$ まで 2 本の縞は**完全に直線のまま**。$y$ 方向の偏差は数値上ちょうど 0、ノルムは 224 を保持。
- **右（`04_darkSoliton_whiteNoise`）**：振幅 $10^{-2}$ の**白色ノイズ**を種に入れる。縞が波打ち、**離散的な点渦の集団＝量子乱流**へ崩壊。

白色ノイズは `#codeStream` の中で OpenFOAM の `Random` クラスを使って与えます（cos 摂動をノイズに置き換えるだけ）。

この対照実験から、次のことが確かめられます。

> **清浄なダークソリトンは「線形不安定だが定常」な厳密解**。放っておくだけでは崩壊せず、**不安定モードを励起する種（摂動）が必要**。右の崩壊は数値解法が勝手に作ったものではなく、入れた種が snake instability を通じて育った物理的な現象である。

カマキリ記事が決定論的な摂動ではなく白色ノイズを使うのは、全波長を一様に励起して**最速成長モード**が勝ち、**不規則・乱流的**な渦配置になるからです。今回の統一ソルバで、その機構をそのまま再現できました。

## 実行手順

```bash
# Case 01（決定論的 cos 摂動）
openfoam2512 -c 'cd run/01_darkSoliton_realTime && blockMesh && schrodingerFoam'
openfoam2512 -c 'cd run/01_darkSoliton_realTime && foamToVTK -fields "(magSqrPsi phase)" -ascii'
python3 tools/render.py run/01_darkSoliton_realTime/VTK figures/01_darkSoliton_realTime/density magSqrPsi

# Case 04（種なし / 白色ノイズ）を同条件で比較
openfoam2512 -c 'cd run/04_darkSoliton_noSeed     && blockMesh && schrodingerFoam'
openfoam2512 -c 'cd run/04_darkSoliton_whiteNoise && blockMesh && schrodingerFoam'
```

## シリーズのまとめ（全7回）

- `laplacianFoam` を **①場を2本 ②式をGPの連立 ③実時間はCN反復 ④虚時間は陰的＋規格化**の4点で改造し、GP 方程式の統一ソルバ `schrodingerFoam` を作った（第1回）。過去に発散した原因は前進オイラーの無条件不安定で、Crank–Nicolson でノルム保存させて克服した。
- 線形の検証として、**量子トンネル効果**（第2回）、**調和振動子のコヒーレント状態**（第3回）、**量子 vs 古典の自由膨張**（第4回）を 1〜2 次元で再現した。
- **虚時間発展**で基底状態を作れることを Thomas–Fermi 近似との一致で検証（第5回）。
- 相互作用ありでは、**走る渦対**（第6回）と**ダークソリトンの崩壊**（第7回）を再現。崩壊には**種（摂動）が必要**で、白色ノイズを入れると量子乱流まで崩壊することを対照実験で示した。

コード一式：<https://github.com/kamakiri1225/schrodingerFoam>
