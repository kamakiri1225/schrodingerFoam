# dft02_helium — ヘリウム原子（Kohn–Sham DFT・M2：Hartree結合）

`ksdftFoam` の第2検証ケース。**電子2個（1s²）のHe**で、電子間の平均反発
（Hartreeポテンシャル）を **Poisson方程式で毎ステップ結合**する。交換相関は
まだ無し（M3）。

## 解いている式

    [ -1/2 ∇² + Vext + V_H ] ψ = ε ψ ,   n = 2ψ² ,   ∇²V_H = -4π n

- Vext = -2/√(r²+0.4²)（Z=2、ソフトCoulomb）
- V_H の境界条件：遠方単極子 **V_H = Nₑ/r = 2/r** を箱の壁に Dirichlet で与える
  （`0/VH` の `#codeStream`。境界条件内では `dict` がパッチのサブ辞書なので
  **`dict.topDict()`** でフィールドのトップ辞書に上がるのが要点）
- SCF は虚時間ステッピングが兼ねる（明示的な密度混合なしで安定に収束）
- 全エネルギー：E = 2ε − E_H（二重計上補正、E_H = ½∫V_H n dV）

## 結果（検証）

| | ε [Ha] | E_tot [Ha] | 参照との差 |
|---|---|---|---|
| ksdftFoam 64³ | -0.08294 | -1.14917 | 7.4e-3 |
| ksdftFoam 96³ | -0.08363 | -1.14550 | 3.7e-3 |
| Richardson 外挿 | | **-1.1426** | 8e-4 |
| 径方向Hartree-SCF参照（同モデル） | -0.0842 | **-1.1418** | — |

誤差は刻み比どおり半減（2次収束）。外挿後の残差 ~0.8 mHa は箱 ±12 Bohr の切断
（参照は R=25 の球）で、軌道が浅く（ε≈-0.08）広がるための有限サイズ効果。

## 物理メモ：自己相互作用誤差（SIE）

この Hartree 近似（V_H を全密度から作る）は電子が**自分自身の反発まで感じる**。
純Coulomb極限では E≈-1.95 Ha となり Hartree–Fock の -2.86 Ha から大きく外れる。
2電子系では交換が E_H のちょうど半分を打ち消すためで、**M3 の LDA 交換が
これを近似的に補正する**——DFTに V_xc が不可欠な理由が数値でそのまま見える。

## 実行

```bash
openfoam2512 -c 'cd run/dft02_helium && blockMesh && ksdftFoam'
```

前：`run/dft01_hydrogen`（M1）。次：M3（`xc slater` を有効化、LDA交換）。
解説記事：`blog/08_ksdftFoam.md`。
