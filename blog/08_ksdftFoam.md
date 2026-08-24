<!-- 第8回：OpenFOAMで密度汎関数法（Kohn–Sham DFT）— ksdftFoam の実装解説 -->

> **シリーズ「OpenFOAM でシュレーディンガー方程式を解く」拡張編・第8回。**
> ①ソルバ自作／②トンネル効果／③調和振動子／④量子vs古典／⑤虚時間発展／⑥走る渦対／⑦2次元量子乱流／**⑧密度汎関数法 DFT（本記事）**
> リポジトリ：<https://github.com/kamakiri1225/schrodingerFoam>／ソルバ：`ksdftFoam/`／ケース：`run/dft01_hydrogen`
> 発端：X での @dc1394 さんとの議論。「多電子系はKohn–Sham方程式になる。ポテンシャルが密度の関数（非線形）で、Hartree項は別途Poisson方程式、SCFで反復収束。水素から、まずLDAで」

GP方程式（第7回まで）に続き、**多電子系の第一原理計算＝Kohn–Sham DFT** をOpenFOAMで実装します。本記事は「どんな方程式を追加しようとしているか／実際に何を追加したか／OpenFOAMのどのファイルをどう書いたか／それが理論式とどう対応するか／どんなスキームか」を、理論とソースコードの両面から丁寧に解説します。

---

## 1. 追加しようとしている方程式（理論）

原子単位系（$\hbar=m_e=e=1$、長さ=Bohr、エネルギー=Hartree）を使います。

### 1.1 Kohn–Sham 方程式（本体）

多電子のシュレーディンガー方程式（He で6次元、タンパク質で数万次元）は直接解けないので、**電子密度 $n(\mathbf r)$ を基本変数**にする（Hohenberg–Kohn 1964）。すると「相互作用する多電子系」が「有効ポテンシャル中の**独立な1電子方程式の組**」に置き換わる（Kohn–Sham 1965）：

$$
\Big[-\tfrac{1}{2}\nabla^2 + V_{\rm eff}(\mathbf r)\Big]\psi_i(\mathbf r) = \varepsilon_i\,\psi_i(\mathbf r),
\qquad
n(\mathbf r)=\sum_i f_i\,|\psi_i(\mathbf r)|^2
$$

$f_i$ は占有数（閉殻なら2）。形は第1回で解いたシュレーディンガー方程式と同じ＝**拡散方程式＋ポテンシャル項**です。

### 1.2 有効ポテンシャル（3つの部品）

$$
V_{\rm eff} = \underbrace{V_{\rm ext}}_{\text{核}} + \underbrace{V_H}_{\text{Hartree}} + \underbrace{V_{xc}}_{\text{交換相関}}
$$

**(a) 核のCoulomb**：$V_{\rm ext}=-Z/r$。格子上で特異点は表現できないため**ソフトCoulomb**で丸める：

$$
V_{\rm ext}(\mathbf r) = -\frac{Z}{\sqrt{r^2+a^2}}
$$

**(b) Hartree ポテンシャル**（電子どうしの平均静電反発）。積分で書くと大変そうだが、実は**Poisson方程式**：

$$
V_H(\mathbf r)=\int\frac{n(\mathbf r')}{|\mathbf r-\mathbf r'|}d\mathbf r'
\quad\Longleftrightarrow\quad
\boxed{\ \nabla^2 V_H = -4\pi n\ }
$$

これは `laplacianFoam` がそのまま解ける形＝**OpenFOAMの本領**。

**(c) 交換相関（LDA）**。量子多体効果をまとめた項。局所密度近似では**密度の局所代数関数**：

$$
V_x(\mathbf r) = -\left(\frac{3\,n(\mathbf r)}{\pi}\right)^{1/3}
\qquad(\text{Slater交換。相関は後で追加})
$$

### 1.3 自己無撞着場（SCF）

$V_{\rm eff}$ が解 $\psi$ 自身の密度に依存する**非線形連立**なので、一発では解けない：

$$
n \to (V_H, V_{xc}) \to V_{\rm eff} \to \text{KS方程式を解く} \to n_{\rm new} \to \cdots \text{（収束まで反復）}
$$

### 1.4 全エネルギー（二重計上の補正）

軌道エネルギーの和は反発を二重に数えるので補正する：

$$
E = \sum_i f_i\,\varepsilon_i \;-\; E_H \;+\;\Big(E_{xc}-\!\int V_{xc}\,n\,d\mathbf r\Big),
\qquad E_H=\tfrac12\!\int V_H\,n\,d\mathbf r
$$

### 1.5 解き方：虚時間発展（第5回の再利用）

KS方程式は固有値問題。第5回と同じく**虚時間の勾配流**で基底状態へ落とす：

$$
\frac{\partial\psi}{\partial\tau} = -H\psi = \tfrac12\nabla^2\psi - V_{\rm eff}\,\psi
\qquad+\ \text{毎ステップ規格化}\ \int|\psi|^2 dV=1
$$

**ポイント**：虚時間ステップごとに $V_{\rm eff}$ を現在の密度から更新すれば、**時間発展がそのままSCF反復を兼ねる**（時間刻みが減衰付き混合の役割）。

---

## 2. 実際に追加したか（現状）

| 段階 | 内容 | 状態 |
|---|---|---|
| **M1** | 水素（核Coulombのみ、$\varepsilon\to-0.5$ Ha） | ✅ **実装・検証完了**（§6） |
| **M2** | He＋Hartree（Poisson結合） | ✅ **実装・検証完了**（§7） |
| **M3** | LDA交換＋SCF完成 | ✅ **実装・検証完了**（§7） |
| M4 | 複数軌道＋Gram–Schmidt直交化（Li, Be） | 未着手 |
| M5 | 擬ポテンシャル | 未着手 |

方針：**schrodingerFoamは拡張せず、別ソルバ `ksdftFoam` を新設**（軌道が複数・SCFが主ループ・単位系も別、と構造が違うため。GP側の見通しを壊さない）。

---

## 3. OpenFOAMのどのファイルをどう書いたか

新規作成は次の3ファイル＋ケース一式。すべて `laplacianFoam` 系の流儀です。

```
ksdftFoam/
├── Make/files        ← ビルド対象とバイナリ名
├── Make/options     ← finiteVolume をリンク
├── createFields.H   ← 場（psi, Vext, n, VH, Vxc）と定数の定義
└── ksdftFoam.C      ← メインループ（虚時間＋SCF）
run/dft01_hydrogen/  ← 水素ケース（メッシュ・初期値・辞書）
```

### 3.1 `createFields.H` — 場の定義

**KS軌道 $\psi$**（基底状態は実関数に取れるので、GPと違い実部1本でよい）：

```cpp
volScalarField psi
(
    IOobject("psi", runTime.timeName(), mesh,
             IOobject::MUST_READ, IOobject::AUTO_WRITE),
    mesh
);
```

**電子密度 $n=f\,\psi^2$**（式1.1の密度。`occupation` が $f$）：

```cpp
volScalarField n( IOobject("n", ...), occupation*sqr(psi) );
```

**Hartree／交換相関ポテンシャル**（M1では0のまま。`READ_IF_PRESENT` なのでケース側が `0/VH` で境界条件を与えられる）：

```cpp
volScalarField VH ( IOobject("VH",  ..., READ_IF_PRESENT, ...), mesh,
                    dimensionedScalar("VH", dimensionSet(0,0,-1,0,0,0,0), 0.0),
                    "zeroGradient" );
```

**次元の工夫**：OpenFOAMの次元チェッカを通すため、`laplacianFoam` の慣例（$D\sim$ m²/s、ポテンシャル$\sim$1/s、$\psi$無次元）を流用。数値的には 1 m＝1 Bohr、1 (1/s)＝1 Ha と読み替えます。Poisson用の $4\pi$ も次元付きで定義：

```cpp
const dimensionedScalar fourPi
(
    "fourPi", dimensionSet(0, -2, -1, 0, 0, 0, 0),   // 1/(m^2 s)
    4.0*Foam::constant::mathematical::pi
);
```

これで `laplacian(VH) [1/(s·m²)] == -fourPi [1/(s·m²)] * n [-]` の次元が合います。

### 3.2 `ksdftFoam.C` — メインループと理論式の対応

ループ1周＝「SCF 1回＋虚時間1ステップ」。番号は §1 の式と対応：

```cpp
while (runTime.loop())
{
    // --- 1. 密度   n = f |psi|^2                    （式1.1の右）
    n = occupation*sqr(psi);

    // --- 2. ポテンシャル更新（SCF）
    if (hartreeOn)
    {
        // lap(VH) = -4 pi n                          （式1.2b: Poisson）
        solve(fvm::laplacian(VH) + fourPi*n);
    }
    if (xcType == "slater")
    {
        // Vx = -(3 n / pi)^{1/3}                     （式1.2c: LDA交換）
        Vxc = -...*cbrt(3.0*max(n, VSMALL)/pi);
    }
    const volScalarField Veff(Vext + VH + Vxc);       //（式1.2）

    // --- 3. 虚時間1ステップ（後退Euler・完全陰解法）  （式1.5）
    solve
    (
        fvm::ddt(psi)
     ==
        D*fvm::laplacian(psi)      // +1/2 ∇²ψ   （D=0.5）
      - fvm::Sp(Veff, psi)         // −Veff ψ    （陰的な吸収項）
    );

    // --- 4. 規格化  ∫ψ²dV = 1                       （式1.5の規格化）
    const dimensionedScalar norm(fvc::domainIntegrate(sqr(psi)));
    psi *= sqrt(targetNorm/norm);

    // --- 5. 軌道エネルギー ε = <ψ|H|ψ>
    const volScalarField Hpsi(-D*fvc::laplacian(psi) + Veff*psi);
    const scalar eps = (fvc::domainIntegrate(psi*Hpsi)/targetNorm).value();

    // --- 6. 全エネルギー（二重計上補正）              （式1.4）
    scalar Etot = occupation*eps;
    if (hartreeOn) Etot -= 0.5*fvc::domainIntegrate(VH*n).value();  // −E_H
    ...

    // --- 7. 収束判定：|Δε| < tol で終了
}
```

**理論式↔コードの対応表**：

| 理論 | コード |
|---|---|
| $-\tfrac12\nabla^2\psi$ | `D*fvm::laplacian(psi)`（D=0.5、陰的） |
| $V_{\rm eff}\,\psi$ | `fvm::Sp(Veff, psi)`（係数場付きの陰的線形源） |
| $\partial\psi/\partial\tau$ | `fvm::ddt(psi)`（Euler＝後退差分） |
| $\nabla^2 V_H=-4\pi n$ | `solve(fvm::laplacian(VH) + fourPi*n)` |
| $\int(\cdot)\,dV$ | `fvc::domainIntegrate(...)`（体積重み付き総和） |
| $\varepsilon=\langle\psi\|H\|\psi\rangle$ | `domainIntegrate(psi*Hpsi)/targetNorm` |
| 規格化 | `psi *= sqrt(targetNorm/norm)` |

`fvm::`（陰的＝行列に組む）と `fvc::`（陽的＝場を評価）の使い分けが肝で、時間発展は`fvm`で無条件安定に、エネルギー評価は`fvc`で後処理的に行っています。

### 3.3 ケース側：`run/dft01_hydrogen`

**ソフトCoulomb**（`0/Vext`、`#codeStream` でセル中心座標から生成）：

```cpp
const scalar Z = 1.0;
const scalar a = 0.4;
V[i] = -Z/Foam::sqrt(x*x + y*y + z*z + a*a);
```

**初期推定**（`0/psi`）：ガウシアン `exp(-0.5 r²)`。虚時間が勝手に 1s 軌道へ整えるので雑でよい（第5回と同じ思想）。境界は箱の壁で `ψ=0`（fixedValue）。

**メッシュ**（`system/blockMeshDict`）：$[-12,12]^3$ Bohr、$64^3$ 一様六面体。

**辞書**（`constant/dftProperties`）：

```
D               0.5;      // 1/2（原子単位の運動項）
occupation      1.0;      // 水素は電子1個
hartree         false;    // M1: 裸の核のみ
xc              none;
convergenceTol  1e-8;     // |Δε| がこれ未満で停止
targetNorm      1.0;      // ∫ψ²dV = 1
```

---

## 4. スキーム（数値解法の選択と理由）

| 項目 | 選択 | 理由 |
|---|---|---|
| 時間積分 | **後退Euler**（`ddtSchemes Euler`＋全項`fvm`） | 虚時間は拡散型で剛性が強い。陰解法なら**無条件安定**で大きな $\Delta\tau=0.1$ が使え、少ステップで収束（水素は211step）。精度は不要（基底状態に落ちれば良い） |
| 空間離散化 | `Gauss linear orthogonal`（2次中心差分相当） | 一様直交六面体なので補正不要の最速・対称。**2次精度**は§6のメッシュ収束で実証 |
| 線形ソルバ | `PCG` + `DIC` 前処理 | 行列は `ddt−laplacian+Sp` で**対称正定値**。共役勾配が最適（毎step 7〜11回で収束） |
| Poisson（$V_H$） | 同じく `PCG/DIC` | Laplacianは対称。前ステップの $V_H$ が初期推定になるのでSCFが進むほど速い |
| SCF | **虚時間ステップ＝減衰付き反復** | 明示的な密度混合なしでも小さい $\Delta\tau$ が混合の役割。不安定になったら線形混合→Pulay(DIIS) を追加予定 |
| 固有値 | 虚時間＋規格化（逆冪乗法と等価） | 最低固有状態だけ欲しいので十分。励起状態（M4）は直交化を足す |

---

## 5. GP（schrodingerFoam）との違い

| | schrodingerFoam（第1〜7回） | ksdftFoam（本記事） |
|---|---|---|
| 波動関数 | 複素1本（Psire+i·Psiim） | **実数の軌道**（基底状態は実に取れる） |
| 非線形項 | $g\|\psi\|^2$（自分自身） | $V_H[n]+V_{xc}[n]$（**密度経由**、$V_H$はPoisson） |
| 主ループ | 実時間発展（物理） | **SCF収束**（虚時間は道具） |
| 単位 | $\xi$, $\hbar/gn_0$ | **原子単位**（Bohr, Hartree） |
| 答えの検証 | ソリトン/渦のダイナミクス | **全エネルギーの数値**（厳密解・文献値と比較） |

「拡散方程式＋ポテンシャル＋虚時間」という**同じ骨格**から、GPもKS-DFTも生える——というのがこのシリーズの主張です。

---

## 6. M1 検証結果（水素原子）✅

**二重の検証**をしました。まず**独立参照値**：同じソフトCoulomb（a=0.4）の径方向1次元固有値問題をPython（`scipy.eigh_tridiagonal`、4万点）で解くと $\varepsilon_{\rm exact}=-0.380872$ Ha。この径方向ソルバ自体も純Coulombで $-0.500000$（解析解）を再現することを確認済み。

その上で `ksdftFoam` のメッシュ収束：

| メッシュ | ε [Ha] | 誤差 |
|---|---|---|
| 64³ (dx=0.375) | −0.38388 | 3.0×10⁻³ |
| 96³ (dx=0.25) | −0.38224 | 1.4×10⁻³ |
| **Richardson外挿** | **−0.38093** | **6×10⁻⁵** |

誤差比 3.0/1.4≈2.2 ≒ (96/64)²=2.25 → **狙い通りの2次精度**で、外挿値は厳密値と5桁一致。虚時間の収束も単調（τ≈21、211ステップ、|Δε|<10⁻⁸）でした。

```
tau = 0.1   eps = -0.2366
tau = 5     eps = -0.38156
tau = 21.1  eps = -0.38388   ← Converged
```

---

## 7. M2・M3 検証結果（He：Hartree → LDA交換）✅

**M2（He＋Hartree）**：`occupation 2`・`Z=2`・`hartree true`。$V_H$ の境界条件は
遠方単極子 $V_H=N_e/r=2/r$ を箱の壁に Dirichlet（`0/VH` の `#codeStream`。
境界条件内では `dict` がパッチのサブ辞書なので **`dict.topDict()`** で上がるのが
実装の要点——素朴にキャストすると SEGV する）。

| M2 | E_tot [Ha] |
|---|---|
| ksdftFoam 64³ / 96³ | -1.14917 / -1.14550 |
| Richardson 外挿 | **-1.1426** |
| 径方向Hartree-SCF参照 | **-1.1418** |

SCF＝虚時間ステッピングは**明示的な密度混合なしで単調収束**した（小さい Δτ が
減衰付き混合の役割を果たす）。

**M3（＋LDA交換）**：dft02 から辞書1行 `xc slater;` だけ。

| M3 | E_tot [Ha] |
|---|---|
| ksdftFoam 64³ / 96³ | -1.62681 / -1.62099 |
| Richardson 外挿 | **-1.6163** |
| 径方向LDA-SCF参照 | **-1.61473** |

どちらも2次収束・外挿後の残差は箱切断（~1 mHa）のみ。

**物理のハイライト＝自己相互作用誤差（SIE）**：Hartree近似（$V_H$ を全密度から
作る）は電子が**自分自身の反発まで感じる**。

| モデル（a=0.4） | E_tot [Ha] |
|---|---|
| Hartree のみ（M2） | -1.142 |
| **+ LDA交換（M3）** | **-1.615** |

交換が SIE を **0.47 Ha** も補正した（2電子系では厳密には交換が $E_H$ の半分を
相殺）。「なぜ DFT に $V_{xc}$ が要るのか」が、実装するとエネルギーの数値で
そのまま見える。相関（VWN/PZ）と純Coulomb文献比較（LDA He ≈ -2.83 Ha）は
M5（擬ポテンシャル/メッシュ細分）と合わせて。

---

## 8. まとめ

- KS-DFT は「**KS方程式（虚時間）＋Poisson（Hartree）＋密度の局所関数（LDA）＋SCF**」に分解でき、どれもOpenFOAMの既存機能（`fvm::ddt/laplacian/Sp`、`domainIntegrate`）で書ける。
- M1（水素）は**2次精度収束と厳密値との5桁一致**まで検証完了。
- 計画全体は `docs/research_plan_ksdft_openfoam.md`、ケースは `run/dft01_hydrogen/README.md` 参照。

次回はM2/M3（He：Poissonの結合とLDA交換、SIEがどれだけ効くかの数値実験）を報告します。
