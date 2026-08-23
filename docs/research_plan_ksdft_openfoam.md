# 計画：OpenFOAM で Kohn–Sham DFT（密度汎関数法）を解く

> 発端：X での @dc1394 さんとの議論（2026-08-23）。「シュレーディンガー方程式は
> 拡散方程式に虚数が入っただけ。多電子系では Kohn–Sham 方程式になり、ポテンシャルが
> 密度の関数（非線形）＋Hartree項は別途 Poisson 方程式で決まり、SCF で反復収束させる」。
> 助言：**水素原子から**、**まず局所密度近似（LDA）**で。ライブラリ群（統一ソルバ）で
> できるか試す。→ 本計画は `schrodingerFoam` の自然な拡張として設計する。

## 0. なぜ OpenFOAM と相性が良いか（狙い）

Kohn–Sham DFT は結局 **3つの PDE の連立＋反復** に落ちる：

1. **KS 方程式**（固有値問題）：`[-½∇² + V_eff] ψ_i = ε_i ψ_i`
   → **虚時間発展**で基底状態が得られる（`schrodingerFoam` の `imaginaryTime` が既にこれ）。
2. **Hartree ポテンシャル**（Poisson 方程式）：`∇²V_H = -4π n`
   → **OpenFOAM の本領**（`laplacianFoam` 相当、`fvm::laplacian`）。
3. **交換相関ポテンシャル** `V_xc[n]`（LDA では密度の局所代数関数）
   → `volScalarField` の式一発。

これらを **SCF（自己無撞着場）ループ**で回すだけ。「拡散方程式ソルバに機能を足す」
という本シリーズの思想そのまま。原子単位系（ħ=m_e=e=1, 4πε₀=1, 長さ=Bohr,
エネルギー=Hartree）で `D=½`。

## 1. 方程式（原子単位）

- 電子密度：`n(r) = Σ_i f_i |ψ_i(r)|²`（`f_i`＝占有数、閉殻なら2）。
- 有効ポテンシャル：`V_eff = V_ext + V_H + V_xc`。
  - 外部（原子核）：`V_ext = -Σ_A Z_A / |r - R_A|`。
  - Hartree：`V_H(r) = ∫ n(r')/|r-r'| dr'` ⇔ `∇²V_H = -4π n`。
  - 交換（LDA, Slater）：`V_x = -(3n/π)^{1/3}`、`E_x = -¾(3/π)^{1/3}∫n^{4/3}dr`。
  - 相関（LDA）：VWN もしくは Perdew–Zunger(1981) パラメタ化（一様電子ガス）。
- 全エネルギー：
  `E = Σ_i f_i ε_i - ½∫∫ n n'/|r-r'| dr dr' + (E_xc - ∫ V_xc n dr) + E_nn`
  （二重計上の補正。`E_nn` は核間反発、原子1個なら0）。

## 2. マイルストーン

### M1 — KS コア：水素原子（電子1個、e-e 相互作用なし）
- 3D メッシュ上で `[-½∇² - Z/r] ψ = ε ψ` を**虚時間発展**で解く（`imaginaryTime` を
  ほぼ流用、`Vext` を 3D Coulomb に）。規格化 `∫|ψ|²=1`。
- **検証**：`E → -0.5 Ha`、`ψ ∝ e^{-r}`（1s）。`Z=2`(He⁺) で `-2.0 Ha` も。
- **課題＝Coulomb 特異点** `-Z/r`：まず**ソフトCoulomb** `-Z/√(r²+a²)` で回避して機構を
  確立 → 収束と `a→0` 外挿、または核近傍メッシュ細分で all-electron を検討。

### M2 — Hartree 結合（Poisson）：ヘリウム（電子2個・1軌道二重占有）
- 各 SCF で `∇²V_H = -4π n` を解く（`fvm::laplacian(V_H) == -4π n`）。
- **Poisson の境界条件**：箱の境界で `V_H ≈ N_e/r`（単極子, Dirichlet）か十分大きな箱。
- **検証**：Hartree のみ（XC 抜き）で He の既知値と比較。密度が2電子で正しく広がるか。

### M3 — LDA 交換相関＋SCF 完成
- `V_xc = V_x(n) + V_c(n)`（Slater 交換＋VWN/PZ 相関）を各 SCF で更新。
- **SCF ループ**：`n 推定 → V_H, V_xc, V_eff → KS を虚時間で緩和 → 新 n → 混合 → 収束判定`。
  混合は線形（`n ← (1-α)n_old + α n_new`, α~0.3）→ 後で Pulay(DIIS)。
- **検証**：He 全エネルギー LDA ≈ **-2.83 Ha**（厳密 -2.90、LDA 誤差込み）。
  水素は LDA だと自己相互作用誤差で `~-0.45 Ha`（この差の議論自体が教材になる）。

### M4 — 複数軌道＋直交化：Li, Be …
- 占有軌道を複数持ち、虚時間の各ステップ後に **Gram–Schmidt 直交化**（低い軌道に対して）。
- スピン／占有数の扱い（まず閉殻・スピン非分極 LSDA なし）。
- **検証**：Be(1s²2s²) など小原子の全エネルギーを LDA 参照値と比較。

### M5 — 核の扱いの高度化
- **擬ポテンシャル**（norm-conserving; 内殻を除いて滑らかな `V_ext`）を導入して
  Coulomb 特異点と内殻の急変を回避（実空間 DFT の定石）。または AE＋メッシュ細分。

## 3. ソルバ構成（`ksdftFoam` もしくは `schrodingerFoam` の mode 追加）

再利用：`imaginaryTime`（軌道緩和）／`normalize`＋`targetNorm`（規格化）／`Vext`／
`fvc::laplacian`／`fvc::domainIntegrate`（`ε_i=⟨ψ|H|ψ⟩`, エネルギー積分）。

新規フィールド／処理：
- `PtrList<volScalarField> psi`（軌道 `ψ_1..ψ_N`。まず実数で可＝基底状態は実に取れる）。
- `volScalarField n`（密度）, `VH`, `Vxc`, `Veff`, `Vext`。
- **Poisson 解**：`solve(fvm::laplacian(VH) == -4.0*pi*n)`（BC は monopole Dirichlet）。
- **V_xc 更新**：`Vxc = -pow(3.0*n/pi, 1.0/3.0) + Vc_LDA(n);`
- **SCF 外側ループ**：密度混合＋収束判定（`‖Δn‖` と `ΔE`）。
- **直交化**（M4）：ステップ後に Gram–Schmidt（`domainIntegrate(psi_i*psi_j)`）。

`constant/dftProperties`（案）：`Z`, `nElectrons`, `xcType(LDA)`, `scfTol`, `mix`,
`softening a`, `poissonBC`。

## 4. 数値上の要注意点

- **Coulomb 特異点**：ソフトCoulomb→擬ポテンシャルへ。all-electron はメッシュ細分必須。
- **境界条件**：軌道 `ψ=0`（十分大きい箱）、`V_H` は単極子 Dirichlet。箱サイズ依存性を確認。
- **SCF 収束**：線形混合は不安定になりうる→ Pulay/DIIS。`α` と収束の関係を記録。
- **3D 化**：本シリーズ初の本格 3D。メッシュ生成・出力（ParaView）・後処理の整備。
- **エネルギーの二重計上補正**（1.の式）を正しく入れないと全エネルギーが合わない。
- **検証は軌道エネルギーではなく全エネルギー**で（LDA 参照値・NIST 原子データと比較）。

## 5. 検証ターゲット（数値）

| 系 | 段階 | 期待値（LDA/厳密） |
|---|---|---|
| H (Z=1, 1電子) | M1（Coulombのみ） | `E=-0.5 Ha`, `ψ∝e^{-r}` |
| He⁺ (Z=2,1電子) | M1 | `E=-2.0 Ha` |
| He (2電子,1s²) | M3（Hartree+LDA） | `E≈-2.83 Ha`（厳密 -2.90） |
| Be, Li | M4 | LDA 参照値 |

## 6. 参考

- Hohenberg–Kohn, Phys. Rev. **136**, B864 (1964)（密度が基本変数；会話中の論文）。
- Kohn–Sham, Phys. Rev. **140**, A1133 (1965)。
- LDA 相関：Vosko–Wilk–Nusair (1980) / Perdew–Zunger (1981)。
- 実空間 DFT：RSDFT（@dc1394 さんの提案）、Chelikowsky らの高次差分実空間法。
- 本リポジトリの再利用元：`schrodingerFoam`（`imaginaryTime` 基底状態探索）。

## 7. 最初の一歩（着手時）

1. 3D 箱メッシュ（例 `[-15,15]³ Bohr`）＋ソフトCoulomb `Vext=-Z/√(r²+a²)`。
2. `imaginaryTime` で 1 軌道を緩和 → `E=⟨ψ|H|ψ⟩` が `-0.5 Ha` に収束するか（水素, a→0）。
3. 通れば M2（Poisson で `V_H`）へ。ここまでで「KS＝虚時間＋Poisson＋代数XC＋SCF」の
   骨格が立つ。
