"""P.6 のイメージ図：エネルギー準位（高温と絶対零度の BEC）、調和振動子アニメの静止画。
出力：img/concept_levels.png、img/still_harmonic.png（sch()・gp() は未使用の旧案）
再現：python3 src/make_eq_images.py
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm

HERE = os.path.dirname(os.path.abspath(__file__)); IMG = os.path.join(os.path.dirname(HERE), "img")
for f in ("/home/kamakiri/.fonts/NotoSansCJKjp-Regular.otf",):
    if os.path.exists(f): fm.fontManager.addfont(f)
plt.rcParams.update({"font.family": "Noto Sans CJK JP", "font.size": 15, "mathtext.fontset": "cm"})
BLUE, RED, GRAY = "#2F6FB0", "#C62828", "#555555"


def frame(ax):
    ax.set_xlim(-5, 5); ax.set_ylim(-0.42, 1.32); ax.axis("off")


def sch():
    fig, ax = plt.subplots(figsize=(6.4, 2.6), dpi=200); frame(ax)
    x = np.linspace(-5, 5, 600); V = 0.035 * x ** 2
    ax.plot(x, V, color=GRAY, lw=2.4); ax.text(3.0, 0.62, "$V(x)$", color=GRAY, fontsize=17)
    rho = np.exp(-(x + 1.0) ** 2 / 0.9); rho = 0.95 * rho / rho.max()
    ax.fill_between(x, 0, rho, color=BLUE, alpha=.30, lw=0); ax.plot(x, rho, color=BLUE, lw=2.4)
    ax.plot([-1.0], [1.08], "o", ms=13, color=BLUE, mec="k")
    ax.text(-0.55, 1.08, "粒子 1 個", fontsize=16, va="center")
    ax.text(0, -0.26, "$|\\psi|^2$：その場所で見つかる確率", fontsize=15, color=BLUE, ha="center", va="center")
    fig.savefig(os.path.join(IMG, "concept_sch.png"), bbox_inches="tight", facecolor="white"); plt.close(fig)


def gp():
    fig, ax = plt.subplots(figsize=(6.4, 2.6), dpi=200); frame(ax)
    x = np.linspace(-5, 5, 600); V = 0.035 * x ** 2
    ax.plot(x, V, color=GRAY, lw=2.4); ax.text(3.75, 0.95, "$V_{\\rm ext}$", color=GRAY, fontsize=17)
    ax.text(3.0, 1.18, "トラップ", color=GRAY, fontsize=14)
    rho = np.clip(1 - (x / 3.2) ** 2, 0, None); rho = 0.95 * rho
    ax.fill_between(x, 0, rho, color=RED, alpha=.22, lw=0); ax.plot(x, rho, color=RED, lw=2.4)
    rng = np.random.default_rng(3); n = 0
    while n < 34:
        px, py = rng.uniform(-3.0, 3.0), rng.uniform(0.05, 0.9)
        if py < 0.95 * (1 - (px / 3.2) ** 2) - 0.06:
            ax.plot(px, py, "o", ms=6, color=RED, mec="k", mew=.6); n += 1
    ax.text(0, -0.26, "← $g|\\psi|^2$：互いに反発して広がる →", fontsize=15, ha="center", va="center")
    ax.text(-4.9, 1.2, "N 個の原子が同じ状態（BEC）", fontsize=15, va="center")
    fig.savefig(os.path.join(IMG, "concept_gp.png"), bbox_inches="tight", facecolor="white"); plt.close(fig)




def levels():
    """エネルギー準位：高温（いろいろな準位に分布）と絶対零度（全員が最低準位 ＝ BEC）。"""
    fig, axs = plt.subplots(1, 2, figsize=(4.8, 2.0), dpi=260)
    rng = np.random.default_rng(1)
    occ_hot = [5, 5, 4, 3, 3, 2, 1, 1]
    for ax, title, occ, col in ((axs[0], "高温", occ_hot, "#7F8C8D"), (axs[1], "絶対零度（$T=0$）", [24] + [0] * 7, RED)):
        ax.set_xlim(-0.9, 10); ax.set_ylim(-0.6, 8.6); ax.axis("off")
        for k in range(8):
            ax.plot([1.2, 9.6], [k, k], color="#9DB8D9", lw=1.4)
        for k, n in enumerate(occ):
            if n == 0: continue
            if col == RED:   # 最低準位に全員（2 段に詰める）
                for x in np.linspace(1.6, 9.2, 12): ax.plot(x, k + 0.28, "o", ms=5, color=col, mec="k", mew=.6)
                for x in np.linspace(1.9, 8.9, 12): ax.plot(x, k + 0.75, "o", ms=5, color=col, mec="k", mew=.6)
                continue
            for x in (np.linspace(2.0, 8.8, n) + rng.uniform(-0.4, 0.4, n) if n > 1 else [rng.uniform(3, 8)]): ax.plot(x, k + 0.28, "o", ms=5, color=col, mec="k", mew=.6)
        ax.annotate("", (0.6, 8.4), (0.6, -0.3), arrowprops=dict(arrowstyle="-|>", color="k", lw=1.4))
        ax.text(0.15, 8.0, "$E$", fontsize=15, va="center", ha="right")
        ax.set_title(title, fontsize=14, color=col if col == RED else "k", pad=2)
    axs[1].text(5.4, 2.3, "全原子が最低準位\n＝ ボース・アインシュタイン\n凝縮体（1 つの波 $\\psi$）", fontsize=12.5, color=RED, ha="center", va="bottom",
                bbox=dict(fc="white", ec="none", pad=2))
    fig.subplots_adjust(wspace=0.08, left=0.01, right=0.99, top=0.88, bottom=0.02)
    fig.savefig(os.path.join(IMG, "concept_levels.png"), bbox_inches="tight", facecolor="white"); plt.close(fig)


def still_harmonic():
    from PIL import Image
    im = Image.open(os.path.join(IMG, "anim_harmonic.gif")); im.seek(0); im.convert("RGB").save(os.path.join(IMG, "still_harmonic.png"))


if __name__ == "__main__":
    levels(); still_harmonic(); print("wrote img/concept_levels.png, img/still_harmonic.png")
