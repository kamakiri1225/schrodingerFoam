# 04_1_darkSoliton_whiteNoise_muShift — 白色ノイズ崩壊（化学ポテンシャルを引いた版）

`04_darkSoliton_whiteNoise` と**同一の物理・初期条件**で、実時間発展だけを
$i\partial_t\psi=(H-\mu)\psi$（$\mu$ を差し引く**回転系**）で解いた版。位相図の背景がもつ
大域位相 $e^{-i\mu t}$ の回転を取り除き、崩壊していく渦の**位相のチラつきを消す**のが目的。

## 04 との違いは1点だけ

`constant/gpProperties` に次を追加している（一様背景 $n_0=1,\ g=1$ なので $\mu=g n_0=1$）。

```c
dynamicMu       false;   // 定数を引く
muShift         1.0;     // mu = g*n0 = 1
```

- `dynamicMu true` にすれば毎ステップ $\mu=\langle\psi|H|\psi\rangle/\langle\psi|\psi\rangle$ を計算して引く。
- 既定（キー無し）は $\mu$ を引かない＝`04` と同じ。

## 仕様

`04_darkSoliton_whiteNoise/README.md` と同じ（$[-16,16]^2$・$128^2$・cyclic、$D=0.5,\ g=1$、
$\Delta t=0.005$、$t\le150$、まっすぐな $\tanh$ 対＋振幅 $10^{-2}$ の白色ノイズ）。差分は $\mu$ の差し引きのみ。

## 効果

- **密度 $\lvert\psi\rvert^2$**：$\mu$ の差し引きは大域位相なので密度・渦崩壊のダイナミクスは不変（`04` と同一）。
- **位相 $\arg\psi$**：`04` は背景が周期 $2\pi/\mu=2\pi$ で全色を巡回して点滅するが、
  本ケースは背景が静止し、**崩壊してできた各渦の $2\pi$ 巻きだけ**がくっきり見える。

理論の詳細は第1回ブログ「位相のチラつきを消す：化学ポテンシャルの差し引き（回転系）」を参照。

## 実行

```bash
openfoam2512 -c 'cd run/04_1_darkSoliton_whiteNoise_muShift && blockMesh && schrodingerFoam'
openfoam2512 -c 'cd run/04_1_darkSoliton_whiteNoise_muShift && foamToVTK -fields "(magSqrPsi phase)" -ascii'
python3 tools/render_density_phase.py run/04_1_darkSoliton_whiteNoise_muShift/VTK \
    figures/04_1_darkSoliton_whiteNoise_muShift/density_phase_muShift.gif "muShift = 1"
```

## 結果

- $t=150$ まで完走、最終ノルムは 224 を保持（`muShift = 1` を確認）。
- 白色ノイズから横方向不安定性が成長して2本のソリトン縞が崩壊、多数の渦核が生成される点は `04` と同一。
- 位相図の背景が静止するため、崩壊後の渦の位相構造が見やすい。
- 密度＋位相GIF：`figures/04_1_darkSoliton_whiteNoise_muShift/density_phase_muShift.gif`
