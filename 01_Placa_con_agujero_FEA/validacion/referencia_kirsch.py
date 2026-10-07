"""
Referencia analítica para la placa con agujero a tensión.
Ing. Juan Diego Rubio Ruiz

- Kirsch (placa infinita): distribución de esfuerzos alrededor del agujero.
- Heywood (placa de ancho finito): factor de concentración de esfuerzos.
- Si resultados/convergencia_malla.csv tiene datos de Ansys, grafica la convergencia de malla.
"""
import sys
sys.stdout.reconfigure(encoding="utf-8")  # consola de Windows (cp1252) no soporta griegas
from pathlib import Path
import csv
import numpy as np
import matplotlib.pyplot as plt

AQUI = Path(__file__).resolve().parent
RES = AQUI.parent / "resultados"
IMG = RES / "imagenes"
IMG.mkdir(parents=True, exist_ok=True)

# Datos del problema
W, d, sigma = 100.0, 20.0, 10.0          # mm, mm, MPa
a = d / 2

# Placa de ancho finito (Heywood)
r = d / W
Kt_neto = 2 + (1 - r) ** 3
sigma_max_ref = Kt_neto * sigma / (1 - r)
print(f"Kt neto (Heywood)       = {Kt_neto:.3f}")
print(f"Kt bruto                = {sigma_max_ref / sigma:.3f}")
print(f"sigma_max de referencia = {sigma_max_ref:.2f} MPa   (Kirsch, placa infinita: {3*sigma:.1f} MPa)")

# Kirsch: sigma_x a lo largo de la línea transversal (theta = 90°), distancia y desde el centro
y = np.linspace(a, W / 2, 300)
q = (a / y) ** 2
sigma_x_kirsch = sigma * (1 + 0.5 * q + 1.5 * q ** 2)

fig, ax = plt.subplots(figsize=(7, 4.5), dpi=150)
ax.plot(y - a, sigma_x_kirsch, lw=2, label="Kirsch (placa infinita)")
ax.axhline(sigma_max_ref, ls="--", color="k", lw=1, label=f"Heywood, ancho finito: {sigma_max_ref:.1f} MPa")
ax.set_xlabel("Distancia desde el borde del agujero [mm]")
ax.set_ylabel("σx [MPa]")
ax.set_title("Esfuerzo σx en la sección transversal del agujero")
ax.grid(alpha=0.3); ax.legend()
# Si existe el Path exportado de Ansys, superponerlo
ruta_path = RES / "path_sigma_x.csv"
if ruta_path.exists():
    filas = list(csv.DictReader(open(ruta_path, encoding="utf-8")))
    if filas:
        yp = np.array([float(f["y_mm"]) for f in filas])
        sp = np.array([float(f["sigma_x_MPa"]) for f in filas])
        ax.plot(yp - a, sp, "o", ms=3, label="Ansys (Path en x = 0)")
        ax.legend()
        # Verificación de equilibrio: la resultante en la sección neta debe ser σ·(W/2)·t
        orden_y = np.argsort(yp)
        t_placa = 5.0                                           # espesor [mm]
        F_neta = t_placa * np.trapezoid(sp[orden_y], yp[orden_y])
        F_apl = sigma * (W / 2) * t_placa
        print(f"Equilibrio: F_neta = {F_neta:.2f} N  vs  F_aplicada = {F_apl:.1f} N"
              f"  (error {abs(F_neta - F_apl) / F_apl * 100:.4f} %)")
fig.tight_layout(); fig.savefig(IMG / "referencia_kirsch.png")
print(f"Gráfica guardada en {IMG / 'referencia_kirsch.png'}")

# Convergencia de malla con los datos de Ansys
datos = list(csv.DictReader(open(RES / "convergencia_malla.csv", encoding="utf-8")))
if datos:
    n = np.array([float(f["nodos"]) for f in datos])
    s = np.array([float(f["sigma_max_MPa"]) for f in datos])
    h = np.array([float(f["tamano_global_mm"]) for f in datos])
    orden = np.argsort(n); n, s, h = n[orden], s[orden], h[orden]

    # Extrapolación de Richardson y GCI (Roache) con tres mallas de razón r = 2
    trio = [np.where(np.isclose(h, v))[0] for v in (2.0, 1.0, 0.5)]
    ext = None
    if all(len(k) for k in trio):
        f3, f2, f1 = (s[k[0]] for k in trio)          # gruesa, media, fina
        r = 2.0
        p_obs = np.log((f2 - f3) / (f1 - f2)) / np.log(r)
        ext = f1 + (f1 - f2) / (r**p_obs - 1)
        gci = 1.25 * abs((f1 - f2) / f1) / (r**p_obs - 1) * 100
        print(f"Orden de convergencia observado p = {p_obs:.2f}")
        print(f"σ_max extrapolado (Richardson)  = {ext:.2f} MPa   (GCI malla fina = {gci:.2f} %)")
        print(f"Diferencia extrapolado vs Heywood = {(ext - sigma_max_ref) / sigma_max_ref * 100:+.2f} %")

    cambio = np.abs(np.diff(s)) / s[1:] * 100     # cambio relativo entre mallas consecutivas
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2), dpi=150)
    a1.semilogx(n, s, "o-", lw=2, label="Ansys Mechanical (σx máx.)")
    a1.axhline(sigma_max_ref, ls="--", color="k", lw=1, label=f"Heywood: {sigma_max_ref:.1f} MPa")
    a1.axhspan(sigma_max_ref * 0.99, sigma_max_ref * 1.01, color="k", alpha=0.08, label="Heywood ±1 %")
    if ext is not None:
        a1.axhline(ext, ls=":", color="tab:red", lw=1.5, label=f"Richardson: {ext:.2f} MPa")
    a1.set_xlabel("Nodos"); a1.set_ylabel("σ_max [MPa]"); a1.set_title("Convergencia de malla")
    a1.grid(alpha=0.3, which="both"); a1.legend(fontsize=8)
    a2.loglog(n[1:], cambio, "o-", color="tab:orange", lw=2)
    a2.set_xlabel("Nodos"); a2.set_ylabel("Cambio respecto a la malla anterior [%]")
    a2.set_title("Variación entre mallas sucesivas"); a2.grid(alpha=0.3, which="both")
    fig.tight_layout(); fig.savefig(IMG / "convergencia_malla.png")
    for ni, si in zip(n, s):
        print(f"  {int(ni):>8d} nodos   σ_max = {si:7.3f} MPa   vs Heywood = {(si - sigma_max_ref) / sigma_max_ref * 100:+.2f} %")
    print(f"Gráfica guardada en {IMG / 'convergencia_malla.png'}")
else:
    print("convergencia_malla.csv aún está vacío: llénalo con tus resultados de Ansys y vuelve a correr el script.")
