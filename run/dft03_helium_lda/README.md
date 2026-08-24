# dft03_helium_lda — ヘリウム原子（Kohn–Sham DFT・M3：LDA交換）

`ksdftFoam` の第3検証ケース。M2（Hartreeのみ）に **LDA交換（Slater）** を足し、
KS-DFT の3点セット（KS方程式＋Poisson＋交換相関＋SCF）が完成する。
dft02 との差分は辞書1行だけ：

    xc              slater;   // was: none

## 解いている式

    [ -1/2 ∇² + Vext + V_H + V_x ] ψ = ε ψ
    n = 2ψ² ,  ∇²V_H = -4π n ,  V_x = -(3n/π)^{1/3}
    E = 2ε - E_H + (E_x - ∫V_x n dV) ,  E_x = -¾(3/π)^{1/3}∫n^{4/3}dV

## 結果（検証）

| | ε [Ha] | E_tot [Ha] |
|---|---|---|
| ksdftFoam 64³ | -0.27821 | **-1.62681** |
| 径方向LDA-SCF参照（同モデル a=0.4） | -0.28061 | **-1.61473** |

差 0.7% は離散化＋箱切断（M1/M2と同水準、2次収束を確認済み）。

## 物理：交換が自己相互作用誤差を補正する

| モデル（a=0.4） | E_tot [Ha] |
|---|---|
| Hartree のみ（dft02、SIEあり） | -1.142 |
| **+ LDA交換（本ケース）** | **-1.615** |

Hartree近似では電子が自分自身の反発まで感じて 0.47 Ha も浅くなっていた。
LDA交換がこれを局所近似で打ち消す（2電子系では厳密には交換が E_H の半分を
相殺）。「なぜ DFT に V_xc が要るのか」がエネルギーの数値でそのまま見える。

相関（VWN/PZ）は未実装。純Coulomb（a→0）での文献比較（LDA He ≈ -2.83 Ha）は
擬ポテンシャル or メッシュ細分（M5）と合わせて行う。

## 実行

```bash
openfoam2512 -c 'cd run/dft03_helium_lda && blockMesh && ksdftFoam'
```

前：`run/dft02_helium`（M2）。次：M4（複数軌道＋直交化）／相関汎関数。
解説記事：`blog/08_ksdftFoam.md`。
