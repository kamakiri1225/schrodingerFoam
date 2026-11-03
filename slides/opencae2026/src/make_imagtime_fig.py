"""Imaginary-time evolution illustration for the opencae2026 slides.

1-D harmonic oscillator V = x^2/2 (hbar = m = omega = 1), finite-difference Hamiltonian.
An arbitrary initial guess is expanded in the eigenfunctions; imaginary-time evolution multiplies
each component by exp(-E_n tau) and the state is renormalised. Left: psi(tau) approaching the
ground state. Right: weight of each eigen-component |c_n(tau)|^2 / sum.
Output: slides/opencae2026/img/fig_imagtime.png
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "img", "fig_imagtime.png")

plt.rcParams.update({"font.family": ["Noto Sans CJK JP", "DejaVu Sans"], "font.size": 18,
                     "axes.titlesize": 20, "axes.labelsize": 19, "xtick.labelsize": 16,
                     "ytick.labelsize": 16, "legend.fontsize": 16})

N, L = 600, 10.0
x = np.linspace(-L, L, N)
dx = x[1] - x[0]
H = (np.diag(np.full(N, 1.0 / dx**2) + 0.5 * x**2)
     + np.diag(np.full(N - 1, -0.5 / dx**2), 1) + np.diag(np.full(N - 1, -0.5 / dx**2), -1))
E, phi = np.linalg.eigh(H)
phi /= np.sqrt(dx)
phi *= np.sign(phi[N // 2 + 5])          # fix sign of each eigenfunction (ground state positive)

psi0 = np.exp(-(x - 2.0) ** 2) + 0.9 * np.exp(-((x + 2.5) ** 2) / 0.4)   # arbitrary guess
psi0 /= np.sqrt((psi0**2).sum() * dx)
c0 = phi.T @ psi0 * dx


def psi_tau(tau):
    c = c0 * np.exp(-E * tau)
    p = phi @ c
    return p / np.sqrt((p**2).sum() * dx), c**2 / (c**2).sum()


fig, (a1, a2) = plt.subplots(1, 2, figsize=(16, 5.6), gridspec_kw={"width_ratios": [1.15, 1]})
cols = ["#9DB8D9", "#5B8FC9", "#2F6FB0", "#1F3A5F"]
for tau, c in zip([0.0, 0.3, 1.0, 4.0], cols):
    p, _ = psi_tau(tau)
    a1.plot(x, p, color=c, lw=3, label=f"τ = {tau:g}")
gs = phi[:, 0] * np.sign(phi[N // 2, 0])
a1.plot(x, gs, "--", color="#C62828", lw=2.4, label="基底状態（厳密）")
a1.set_xlim(-5, 5)
a1.set_xlabel("x")
a1.set_ylabel("ψ")
a1.set_title("ψ(τ)：どんな初期値も基底状態へ")
a1.legend(loc="upper left", frameon=False)

taus = np.linspace(0, 4, 300)
W = np.array([psi_tau(t)[1][:5] for t in taus])
labs = ["n = 0（基底）", "n = 1", "n = 2", "n = 3", "n = 4"]
cl = ["#C62828", "#2F6FB0", "#5B8FC9", "#9DB8D9", "#BBBBBB"]
for n in range(5):
    a2.plot(taus, W[:, n], color=cl[n], lw=3 if n == 0 else 2.2, label=labs[n])
a2.set_xlabel("虚時間 τ")
a2.set_ylabel("成分の割合 |cₙ|² / Σ")
a2.set_title("高いエネルギー成分ほど速く消える")
a2.set_ylim(0, 1.02)
a2.legend(loc="center right", frameon=False)
fig.tight_layout(w_pad=2.5)
fig.savefig(OUT, dpi=110, bbox_inches="tight", facecolor="white")
print("wrote", OUT, "E0..E4 =", np.round(E[:5], 4), "initial weights", np.round(c0[:5] ** 2, 3))


# ---------------------------------------------------------------- animation (anim_imagtime.gif)
def make_gif():
    import io
    from PIL import Image
    taus = np.linspace(0, 4, 61)
    Ehist = []
    for t in taus:
        p, w = psi_tau(t)
        Ehist.append((w * E).sum())
    Ehist = np.array(Ehist)
    frames = []
    for k, t in enumerate(taus):
        p, w = psi_tau(t)
        fig, (b1, b2) = plt.subplots(1, 2, figsize=(16, 5.6), gridspec_kw={"width_ratios": [1.15, 1]})
        b1.fill_between(x, p, color="#9DB8D9", lw=0)
        b1.plot(x, p, color="#1F3A5F", lw=3, label="ψ(τ)")
        b1.plot(x, gs, "--", color="#C62828", lw=2.2, label="基底状態（厳密）")
        bb = b1.twinx()
        bb.plot(x, 0.5 * x**2, color="#888888", lw=1.8, ls=":")
        bb.set_ylim(0, 12.5); bb.set_yticks([]); 
        b1.set_xlim(-5, 5); b1.set_ylim(-0.05, 0.8)
        b1.set_xlabel("x"); b1.set_ylabel("ψ")
        b1.set_title(f"虚時間 τ = {t:4.2f}：凸凹が消えて、谷底の形に落ち着く")
        b1.legend(loc="upper left", frameon=False)
        b2.plot(taus, Ehist, color="#BBBBBB", lw=2)
        b2.plot(taus[:k + 1], Ehist[:k + 1], color="#2F6FB0", lw=3)
        b2.plot([t], [Ehist[k]], "o", color="#2F6FB0", ms=12)
        b2.axhline(E[0], color="#C62828", ls="--", lw=2)
        b2.text(3.95, E[0] + 0.08, "基底エネルギー E₀ = 0.5", color="#C62828", ha="right", fontsize=17)
        b2.set_xlim(0, 4); b2.set_ylim(0, max(Ehist) * 1.08)
        b2.set_xlabel("虚時間 τ"); b2.set_ylabel("エネルギー ⟨H⟩")
        b2.set_title("エネルギーが下がり続け、最低値で止まる")
        fig.tight_layout(w_pad=2.5)
        buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=110, bbox_inches="tight", facecolor="white"); plt.close(fig)
        buf.seek(0); frames.append(Image.open(buf).convert("RGB"))
    pal = frames[len(frames) // 2].quantize(colors=96, method=Image.Quantize.MEDIANCUT)
    q = [f.quantize(palette=pal, dither=Image.Dither.NONE) for f in frames]
    out = os.path.join(HERE, "..", "img", "anim_imagtime.gif")
    q[0].save(out, save_all=True, append_images=q[1:], duration=[120] * (len(q) - 1) + [1800], loop=0, optimize=True)
    print("wrote", out, os.path.getsize(out) / 1e6, "MB", "E(0)=", Ehist[0])


make_gif()
