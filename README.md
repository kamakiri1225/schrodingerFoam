# schrodingerFoam — Gross–Pitaevskii solver for OpenFOAM (v2512)

📊 **発表スライド（ブラウザで閲覧）**：
**<https://kamakiri1225.github.io/schrodingerFoam/docs/slides_schrodingerFoam.html>**
（reveal.js。→/← で章、↓/↑ で章内、`F` 全画面、`?print-pdf` を付けるとPDF化）

OpenFOAM をカスタマイズして、**Gross–Pitaevskii（非線形シュレーディンガー）方程式**
$i\,\partial_t\psi = -D\nabla^2\psi + (V_\mathrm{ext} + g|\psi|^2)\psi$
を解く統一ソルバ。**虚時間発展（初期状態づくり）** と **実時間発展（本計算）** を
`mode` 切替で1つのソルバにまとめ、2次元ダークソリトンの横方向不安定性（snake instability）
から**量子渦の生成・崩壊**までを再現する。

<p align="center">
  <img src="figures/04_darkSoliton_whiteNoise/comparison_noSeed_vs_whiteNoise.gif" width="720"><br>
  <em>左：摂動なし（崩壊しない） ／ 右：白色ノイズ種あり（渦へ崩壊）</em>
</p>

## なぜ OpenFOAM で量子力学なのか
OpenFOAM は CFD（流体解析）のツールとして知られるが、その正体は**偏微分方程式（PDE）を解くための
フレームワーク**である。そして量子力学の基礎方程式であるシュレーディンガー方程式もまた、
時間と空間の微分を含む立派な PDE だ。つまり「馴染みの薄い量子力学」も、方程式の形さえ書き下せれば、
流体解析と同じ土俵で数値的に解ける。本リポジトリはその実証として、GP 方程式を OpenFOAM に載せている。

## 一度は発散した — だから作り直した
とはいえ、素朴に移植すればよいわけではなかった。以前 `laplacianFoam` をそのまま改造した実装
（前進オイラー）は、1ステップ進むごとに波がわずかに増幅してしまい、時間刻みをどれだけ小さくしても
最後は計算が**発散**した。シュレーディンガー方程式は「波の存在確率（ノルム）が時間で増えも減りも
しない」のが物理の大前提なのに、それが数値的に壊れていたのだ。

本ソルバは、実時間発展を **Crank–Nicolson 法**（前後の状態の平均で更新する解き方）に変え、
波の大きさがステップを経ても変わらない＝**ノルムが保存**するようにして安定化した。初期状態づくりに
使う虚時間発展は、ふつうの拡散方程式と同じ形なので完全陰的に解けば無条件に安定になる。

要するに「発散したのは方程式ではなく“時間の刻み方”のせい」で、そこを直したのが作り直しの中身。
増幅率などの数式的な理由は [`notes.md`](notes.md)、読み物としての解説は `blog/` の記事にまとめた。

次期研究テーマ「2次元量子乱流BECの自由膨張」は
[`docs/research_plan_2D_quantum_turbulence_free_expansion.md`](docs/research_plan_2D_quantum_turbulence_free_expansion.md)
に、先行研究、研究ギャップ、診断量、段階的な計算計画を整理している。

## 構成
```
schrodingerFoam/   ソルバ本体（schrodingerFoam.C, createFields.H, Make/）
run/               計算ケース（各フォルダに README.md 仕様書）
  01_darkSoliton_realTime/    ダークソリトン→渦核（決定論的 cos 摂動で種付け）
  02_trap_imaginaryTime/      虚時間：調和トラップ基底状態（Thomas–Fermi 一致）
  03_vortexDipole_realTime/   走る渦対（渦-反渦ダイポール）
  04_darkSoliton_noSeed/      対照：摂動なし → 崩壊しない
  04_darkSoliton_whiteNoise/  カマキリ流：白色ノイズ種 → 不規則な渦崩壊
tools/render.py    VTK → PNG → GIF 可視化
figures/           結果 GIF・図
notes.md           ブログ用の詳しいメモ（理論・実装・比較）
```

## ビルドと実行（OpenFOAM v2512）
```bash
# ソルバのビルド
openfoam2512 -c 'cd schrodingerFoam && wmake'

# 例：ダークソリトン崩壊
openfoam2512 -c 'cd run/01_darkSoliton_realTime && blockMesh && schrodingerFoam'
```

## 参考
- カマキリ（宇宙に入ったカマキリ）「2次元ダークソリトンの崩壊」（GP方程式・Fortran）
- gucong.org "Numerical solving of non-linear Schrödinger equation"

> 元となる物理記事は筆者（リポジトリ所有者）自身のブログ記事に基づく。
