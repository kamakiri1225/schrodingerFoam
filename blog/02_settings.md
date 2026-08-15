<!-- 第2回：OpenFOAM のファイル設定と実行方法 -->

> **全3回シリーズ・第2回。** 前回作った GP ソルバ `schrodingerFoam` を実際に動かすために、OpenFOAM 側の**ファイル設定（`0/`・`constant/`・`system/`）を1つずつ**解説し、**計算の実行方法**を紹介します。
>
> - 第1回：ソルバのカスタマイズと虚時間発展
> - 第2回（本記事）：OpenFOAM のファイル設定と計算の実行方法
> - 第3回：各ケースの計算結果（ダークソリトンの量子渦崩壊）
>
> リポジトリ：<https://github.com/kamakiri1225/schrodingerFoam>

各ケースは OpenFOAM の標準構成 `0/`（初期条件）・`constant/`（物性）・`system/`（数値設定）を持ちます。ここでは `run/01_darkSoliton_realTime` を例に、設定ファイルを1つずつ見ていきます。

## `constant/gpProperties` — 物理係数とモード

このソルバ固有の設定ファイル（自作）。物理係数と計算モードをここで決めます。

```c
// 分散係数 D = hbar / (2 m)
D               D  [0 2 -1 0 0 0 0]  0.5;

// 非線形結合定数 g（>0 で斥力）
g               g  [0 0 -1 0 0 0 0]  1.0;

// 時間発展モード: realTime | imaginaryTime
mode            realTime;

// Crank-Nicolson の Picard 反復回数（realTime）
nCorrectors     4;

// 規格化（主に imaginaryTime の基底状態探索用）
normalize       false;
targetNorm      1.0;

// 場が変化しなくなったら停止（imaginaryTime 用）。0 = 停止しない
convergenceTol  0;
```

- `mode` を `imaginaryTime` にすれば、同じソルバが虚時間発展（初期状態づくり）に変わります。
- 虚時間で基底状態を探すときは `normalize true`・`targetNorm 20`・`convergenceTol 1e-5` のように使います。

## `system/blockMeshDict` — 2次元・周期境界のメッシュ

snake instability を扱うには**周期系**が必須。$[-16,16]^2$ の正方領域を $128\times128$ に切り、z 方向は1層だけの疑似2次元にします。

```c
Lx 16; Ly 16; dz 0.25;   // 半サイズと厚み
nx 128; ny 128;          // 分割数（dx = 32/128 = 0.25）

blocks
(
    hex (0 1 2 3 4 5 6 7) ($nx $ny 1) simpleGrading (1 1 1)
);

boundary
(
    left  { type cyclic; neighbourPatch right; faces ( (0 4 7 3) ); }
    right { type cyclic; neighbourPatch left;  faces ( (1 2 6 5) ); }
    down  { type cyclic; neighbourPatch up;    faces ( (0 1 5 4) ); }
    up    { type cyclic; neighbourPatch down;  faces ( (3 7 6 2) ); }
    frontAndBack { type empty; faces ( (0 3 2 1) (4 5 6 7) ); }
);
```

- `cyclic`：x・y の両方向を**周期境界**に（ソリトンが箱の端で反射せず、無限系のように振る舞う）。
- `frontAndBack` は `empty`：z 方向には解かない＝2次元計算の指定。

> 実ファイルでは境界の各 patch を複数行に展開して読みやすくしています（上は紙面用に1行化）。

## `system/fvSchemes` — 離散化スキーム

```c
ddtSchemes        { default Euler; }                    // 時間：1次前進（陰的側で安定）
gradSchemes       { default Gauss linear; }
divSchemes        { default none; }                     // 移流項はないので none
laplacianSchemes  { default Gauss linear corrected; }   // ラプラシアン：2次・非直交補正
snGradSchemes     { default corrected; }
```

移流（div）項がないのが GP 方程式の特徴で、`divSchemes` は `none` で構いません。

## `system/fvSolution` — 線形ソルバ

虚時間モードで `fvm::` の行列を解くときに使います（実時間の `fvc::` 主体でも場の名前解決に必要）。

```c
solvers
{
    "(Psire|Psiim)"
    {
        solver          PCG;       // 共役勾配法
        preconditioner  DIC;
        tolerance       1e-9;
        relTol          0;
    }
}
```

## `system/controlDict` — 時間刻みと出力

```c
application     schrodingerFoam;
startTime       0;
endTime         60;
deltaT          0.005;      // CN の反復収束のため dt <~ dx^2/D 程度に
writeControl    runTime;
writeInterval   1;          // t を 1 進むごとに書き出し
writeFormat     binary;
```

- $\Delta t=0.005$：実時間の分離陽的ラプラシアンは Picard 収束のため $\Delta t\lesssim\Delta x^2/D$ が目安（$\Delta x=0.25,\ D=0.5$ なら $\sim0.125$、余裕を見て小さめ）。

## `0/Psire`・`0/Psiim` — 初期条件（`#codeStream` で解析生成）

初期のダークソリトン対を、OpenFOAM 内でその場でコンパイルされる `#codeStream` で作ります。2本の $\tanh$ の溝＋横方向の cos 摂動（snake の種）です。

```cpp
const scalar w = 1.0;    // ソリトン幅（= sqrt(2)*xi）
const scalar x0 = 6.0;   // 2本の溝の半間隔
const scalar A = 0.25;   // 摂動振幅
const scalar lamY = 8.0; // 横方向の摂動波長

forAll(psi, i)
{
    const scalar x = mesh.C()[i].x();
    const scalar y = mesh.C()[i].y();
    const scalar sh = A*Foam::cos(twoPi*y/lamY);   // 溝の横ずれ
    psi[i] = Foam::tanh((x + x0 - sh)/w) * Foam::tanh((x - x0 - sh)/w);
}
```

`Psiim`（虚部）は 0＝速度ゼロの「ブラックソリトン」です。カマキリ流の case 04 では、この摂動を `Random` クラスの**白色ノイズ**に置き換えています。

## 計算の実行方法

```bash
# ① メッシュ生成 → ② 計算実行 を一気に
openfoam2512 -c 'cd run/01_darkSoliton_realTime && blockMesh && schrodingerFoam'
```

- `blockMesh`：`blockMeshDict` を読んでメッシュ（`constant/polyMesh/`）を生成。
- `schrodingerFoam`：`mode` に従って計算。実行中は各ステップの `norm`（や虚時間なら `mu`）が表示され、ノルムが一定なら安定に解けている証拠です。

計算後、結果を VTK に出して Python で GIF 化します。

```bash
# ③ 場を ascii の VTK (.vtu) に書き出し
openfoam2512 -c 'cd run/01_darkSoliton_realTime && foamToVTK -fields "(magSqrPsi phase)" -ascii'

# ④ VTK → PNG連番 → GIF（tools/render.py）
python3 tools/render.py run/01_darkSoliton_realTime/VTK figures/01_darkSoliton_realTime/density magSqrPsi
python3 tools/render.py run/01_darkSoliton_realTime/VTK figures/01_darkSoliton_realTime/phase   phase
```

`render.py` は `internal.vtu` を読み、セル中心座標から一様格子に並べ替えて `matplotlib` でヒートマップ化し、`Pillow` で GIF に結合します。

---

次回（第3回）は、こうして用意したケースを実際に走らせた **各フォルダの計算結果**（ダークソリトンの snake instability →量子渦崩壊、虚時間発展での基底状態、走る渦対、種の有無の比較）を GIF つきで紹介します。
