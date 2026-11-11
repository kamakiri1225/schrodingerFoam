# 04_3_darkSoliton_fromImagTime_realTime — 虚時間で作った初期状態からのダークソリトン崩壊

- `0/`：`../04_2_darkSoliton_imagTime_fixPhase` の τ = 6 の結果に、加法的な複素白色ノイズ（振幅 10⁻²）を加えたもの
  （numpy の乱数、シード 20260814。ノイズの形は case 04 と同じ `no*f1*cos(2πf2)`、虚部は `no*f3*sin(2πf4)`）。
- 実時間発展（Crank–Nicolson）、t = 0〜150、Δt = 0.005。設定は case 04 と同じ。
- 巻き数 ±1 の量子渦：t = 60 で 12 個、t = 100 で 10 個（`slides/opencae2026/src/make_gifs.py` の count_vortices）。
- 実行：`blockMesh && schrodingerFoam`（約 12 分、1 並列）。
