# 01_darkSoliton_realTime — ダークソリトンの横方向不安定性→量子渦

## 何の計算か
2本のダークソリトン（低密度の溝）を初期状態にして、**実時間（realTime）**でGP方程式を解く。
横方向の不安定性（snake instability）でソリトンが数珠状にくびれ、量子渦核へと崩壊する過程を見る。
→ カマキリ記事 `animate-1280x720.gif` の再現がねらい。

## 仕様
| 項目 | 値 |
|---|---|
| モード | `realTime`（Crank–Nicolson + Picard 4反復） |
| 方程式 | $i\partial_t\psi=-D\nabla^2\psi+g|\psi|^2\psi$ |
| 係数 | $D=0.5,\ g=1$（→ $n_0=1,\ \mu=1,\ \xi=0.707$） |
| 領域 | $[-16,16]^2$、z方向 empty（2D） |
| 格子 | 128×128（$\Delta x=0.25$、≈5.6 cells/$\xi$） |
| 境界 | x,y とも **cyclic（周期）** |
| 時間刻み | $\Delta t=0.005$ |
| 終了時刻 | $t=60$（writeInterval 1 → 60フレーム） |
| 初期条件 | $\psi_0=\tanh\frac{x+x_0-s(y)}{w}\tanh\frac{x-x_0-s(y)}{w}$, $\psi_{\rm im}=0$ |
| 　パラメータ | $w=1,\ x_0=6,\ s(y)=A\cos(2\pi y/\lambda_y),\ A=0.25,\ \lambda_y=8$ |

初期の虚部を0（＝**静止したブラックソリトン**）にしている点が重要（後述の「渦が走らない」理由）。

## 実行
```bash
openfoam2512 -c 'blockMesh && schrodingerFoam'
```

## 結果
- ノルム（粒子数）は計算中ずっと保存（≈224 で一定）＝スキームがユニタリ。
- $t\sim35$ で snake instability が明瞭化、溝が密度ノードに分裂。
- 図: `figures/01_darkSoliton_realTime/density.gif`, `figures/01_darkSoliton_realTime/phase.gif`

## 注意
- 静止ブラックソリトン＋対称摂動なので、パターンは**その場で**くびれるだけで並進しない
  （渦対が「走る」のを見たいときは case 03 を参照）。
