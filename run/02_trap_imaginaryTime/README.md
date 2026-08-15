# 02_trap_imaginaryTime — 虚時間発展による調和トラップ基底状態

## 何の計算か
**虚時間発展（imaginaryTime）**で、調和ポテンシャルに捕獲されたBECの**基底状態**を求める。
統一ソルバの「初期状態を作る」機能の検証。解析解 Thomas–Fermi と比較する。

## 仕様
| 項目 | 値 |
|---|---|
| モード | `imaginaryTime`（完全陰的 backward Euler + 毎ステップ規格化） |
| 方程式 | $\partial_\tau\psi=D\nabla^2\psi-(V_{\rm ext}+g|\psi|^2-\mu)\psi$ |
| 係数 | $D=0.5,\ g=1$ |
| 外部ポテンシャル | $V_{\rm ext}=\tfrac12\omega^2 r^2,\ \omega=1$（`0/Vext` を codeStream 生成） |
| 領域 | $[-12,12]^2$、z方向 empty |
| 格子 | 128×128 |
| 境界 | 4辺 wall（$\psi=0$ fixedValue）、front/back empty |
| 時間刻み | $\Delta\tau=0.05$（陰的なので大きく取れる） |
| 規格化 | `normalize true`, `targetNorm 20`（粒子数固定） |
| 収束判定 | `convergenceTol 1e-5`（変化が止まったら自動停止） |
| 初期推定 | ガウシアン $\psi_0=e^{-r^2/2\sigma^2},\ \sigma=3$ |

## 実行
```bash
openfoam2512 -c 'blockMesh && schrodingerFoam'
```

## 結果
- 化学ポテンシャル $\mu\to5.174$ に収束（TF予測 $\mu\approx5.05$ と一致）。
- 密度は反転放物線（Thomas–Fermi）に一致、縁はヒーリング長でなめらか。
- 図: `figures/02_trap_imaginaryTime/groundstate.png`（密度＋径方向プロファイル vs TF）

## 意義
「なぜ虚時間で基底状態に落ちるのか」＝高エネルギー成分ほど $e^{-E\tau}$ で速く減衰し、
規格化を挟むと最低エネルギー状態だけが残るから（詳細は `notes.md` 参照）。
