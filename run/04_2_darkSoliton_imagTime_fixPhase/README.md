# 04_2_darkSoliton_imagTime_fixPhase — 位相を固定した虚時間発展でダークソリトンの初期状態を作る

- 初期値：密度 1 で一様、2 本のソリトン位置（x = ±6）の間だけ位相 π（`0/Psire` = −1 / +1）。
- `constant/gpProperties`：`mode imaginaryTime`、`fixPhase true`（毎ステップ位相を初期値に戻し、密度だけを緩和）、
  `normalize true`、`targetNorm 224`（背景密度 1 の解析的な 2 本のダークソリトンのノルム）。
- 虚時間 τ = 0〜6、Δτ = 0.01。τ = 6 で厳密解 tanh²((x+6))·tanh²((x−6)) との差は最大 0.004、μ → 0.9995 ≈ g n₀。
- ふつうの（位相を固定しない）虚時間発展では、ダークソリトン（励起状態）は消えて一様な基底状態になる。
- 実行：`blockMesh && schrodingerFoam`。結果（τ = 6）は `../04_3_darkSoliton_fromImagTime_realTime/0` の元。
