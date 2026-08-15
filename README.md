# schrodingerFoam — Gross–Pitaevskii solver for OpenFOAM (v2512)

OpenFOAM をカスタマイズして **Gross–Pitaevskii（非線形シュレーディンガー）方程式**
$$ i\,\partial_t\psi = -D\nabla^2\psi + \big(V_\mathrm{ext} + g|\psi|^2\big)\psi $$
を解く統一ソルバ。**虚時間発展（初期状態づくり）** と **実時間発展（本計算）** を
`mode` 切替で1つのソルバにまとめ、2次元ダークソリトンの横方向不安定性（snake instability）
から**量子渦の生成・崩壊**までを再現する。

<p align="center">
  <img src="figures/04_darkSoliton_whiteNoise/comparison_noSeed_vs_whiteNoise.gif" width="720"><br>
  <em>左：摂動なし（崩壊しない） ／ 右：白色ノイズ種あり（渦へ崩壊）</em>
</p>

## なぜ作り直したか
以前の試み（laplacianFoam 改造・前進オイラー）は、シュレーディンガー方程式に対して
**無条件不安定**（増幅率 $\sqrt{1+(E\Delta t)^2}>1$）で発散していた。本ソルバは実時間を
**Crank–Nicolson**（増幅率＝1、ノルム保存）で解き、虚時間を**完全陰的**で解くことで安定化した。
詳細な導出・コード解説・可視化手順は [`notes.md`](notes.md)。

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
