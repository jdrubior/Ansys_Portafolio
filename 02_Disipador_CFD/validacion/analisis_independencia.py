"""
Estudio de independencia de malla del disipador (Fluent).
Ing. Juan Diego Rubio Ruiz

Lee resultados/independencia_malla.csv, calcula el orden observado, la
extrapolación de Richardson y el GCI de la malla fina (Celik et al., 2008),
y genera resultados/imagenes/independencia_malla.png.
"""
import sys, csv
sys.stdout.reconfigure(encoding="utf-8")
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

AQUI = Path(__file__).resolve().parent
RES = AQUI.parent / "resultados"
filas = list(csv.DictReader(open(RES / "independencia_malla.csv", encoding="utf-8")))
filas.sort(key=lambda f: float(f["elementos"]))           # gruesa -> fina
N = np.array([float(f["elementos"]) for f in filas])


def gci(f, N, Fs=1.25, p_max=2.0):
    """f, N ordenados de grueso a fino (3 mallas). Devuelve p observado, p usado, extrapolado y GCI [%]."""
    f3, f2, f1 = f
    r21 = (N[2] / N[1]) ** (1 / 3)
    r32 = (N[1] / N[0]) ** (1 / 3)
    e21, e32 = f2 - f1, f3 - f2
    s = np.sign(e32 / e21)
    p = 2.0
    for _ in range(100):                                  # iteración de punto fijo (Celik)
        q = np.log((r21**p - s) / (r32**p - s))
        p = abs(np.log(abs(e32 / e21)) + q) / np.log(r21)
    p_uso = min(p, p_max)                                 # acotar al orden teórico del esquema
    ext = (r21**p_uso * f1 - f2) / (r21**p_uso - 1)
    gci_fino = Fs * abs((f1 - f2) / f1) / (r21**p_uso - 1) * 100
    return p, p_uso, ext, gci_fino


vars_ = [("R_th_K_W", "$R_{th}$", "K/W", "Resistencia térmica"),
         ("caida_presion_Pa", r"$\Delta p$", "Pa", "Caída de presión")]
fig, axs = plt.subplots(1, 2, figsize=(10, 3.8), dpi=200)
for ax, (col, sim, uni, tit) in zip(axs, vars_):
    f = np.array([float(x[col]) for x in filas])
    p, p_uso, ext, g = gci(f, N)
    print(f"{tit}: p observado = {p:.2f} (se usa {p_uso:.1f}), extrapolado = {ext:.4g} {uni}, GCI_fina = {g:.2f} %")
    ax.plot(N, f, "-o", color="#2a78d6", lw=2, ms=7, mec="white", mew=1.5, zorder=3, label="Fluent")
    ax.axhline(ext, ls="--", color="#52514e", lw=1, label=f"Richardson: {ext:.4g} {uni}")
    for n, v, nm in zip(N, f, ["gruesa", "intermedia", "fina"]):
        ax.annotate(f"{nm}\n{v:.4g}", (n, v), textcoords="offset points", xytext=(0, 9),
                    ha="center", fontsize=8, color="#52514e")
    ax.set_xscale("log"); ax.set_xlabel("Número de celdas"); ax.set_ylabel(f"{sim} [{uni}]")
    ax.set_title(tit, loc="left"); ax.grid(True, which="both", color="#e6e5e0", lw=0.7)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    ax.legend(frameon=False, fontsize=8)
axs[0].set_ylim(1.215, 1.245); axs[1].set_ylim(8.6, 9.9)
fig.tight_layout()
salida = RES / "imagenes" / "independencia_malla.png"
fig.savefig(salida, facecolor="white")
print(f"Gráfica guardada en {salida}")
