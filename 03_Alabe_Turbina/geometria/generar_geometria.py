"""Geometria parametrica del alabe tipo NASA C3X (corte de 12.7 mm de envergadura).

Dimensiones publicadas (NASA CR-174827, alabe C3X):
  cuerda real C = 144.93 mm, cuerda axial = 78.16 mm, altura = 76.2 mm,
  RLE = 11.68 mm, RTE = 1.73 mm, 10 canales de enfriamiento radiales.
El perfil exacto NO se reproduce (la tabla de coordenadas no esta disponible
aqui): se construye un perfil parametrico con esas cuerdas y radios. La posicion
y el diametro de los 10 canales es representativa.
Salidas: perfil_exterior.csv, canales.csv, previa_geometria.png, informe de
verificacion por consola.
"""
import csv
import numpy as np
from shapely.geometry import Polygon, Point
from shapely.ops import unary_union
from scipy.optimize import brentq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

C, CAX = 144.93, 78.16
RLE, RTE = 11.68, 1.73
TMAX = 38.0            # espesor maximo (mm)
GAMMA = np.arccos(CAX/C)  # escalonamiento de la cuerda (cos g = Cax/C)
PHI_DEG = 52.0         # giro de la linea de camber (de diseno)
H_SLICE = 12.7         # espesor del corte (mm)
PITCH = 117.7          # paso entre alabes (mm), valor de la cascada C3X
WALL_MIN = 2.5         # pared minima canal-superficie (mm)
D_MAX = 6.3            # diametro maximo de canal (mm)
S_CH = [0.07, 0.15, 0.24, 0.33, 0.42, 0.51, 0.60, 0.69, 0.77, 0.85]

LC = C - RLE - RTE     # distancia entre centros de LE y TE

def camber(s, PHI, R_ARC):
    a = PHI*(s-0.5)
    x = R_ARC*np.sin(a) + 0.0
    y = -R_ARC*np.cos(a)
    return x, y, np.cos(a), np.sin(a)      # punto y tangente (tx,ty)

from scipy.interpolate import CubicSpline
_T = [0.0, 0.20, 0.40, 0.58, 0.76, 0.90, 1.0]
_H = [RLE, 16.0, TMAX/2, 17.0, 11.5, 5.2, RTE]
_SPL = CubicSpline(_T, _H, bc_type=((1, 0.0), (1, -20.0)))

def semi_espesor(s):
    """Semi-espesor a lo largo de la camber (spline C2, tangente al circulo del LE)."""
    return _SPL(s)

def construir(rot, PHI):
    R_ARC = LC/(2*np.sin(PHI/2))
    s = np.linspace(0, 1, 400)
    x, y, tx, ty = camber(s, PHI, R_ARC)
    nx, ny = -ty, tx
    h = semi_espesor(s)
    ss = np.c_[x + h*nx, y + h*ny]
    ps = np.c_[x - h*nx, y - h*ny]
    c0 = np.array([x[0], y[0]]); c1 = np.array([x[-1], y[-1]])
    poly = Polygon(np.vstack([ss, ps[::-1]])).buffer(0)
    poly = unary_union([poly, Point(*c0).buffer(RLE, 64), Point(*c1).buffer(RTE, 64)])
    poly = poly.buffer(2.0, 64).buffer(-2.0, 64)       # cierra esquinas concavas
    R = np.array([[np.cos(rot), -np.sin(rot)], [np.sin(rot), np.cos(rot)]])
    return poly, R, (x, y, h, nx, ny)

def extension(PHI, rot=-GAMMA):
    poly, R, _ = construir(rot, PHI)
    xy = np.array(poly.exterior.coords)@R.T
    return xy[:, 0].max() - xy[:, 0].min()

# angulo de rotacion para que la cuerda axial sea CAX
rot = -GAMMA
PHI = np.radians(PHI_DEG)
poly, R, (x, y, h, nx, ny) = construir(rot, PHI)
c_ext = np.array(poly.exterior.coords)@R.T
ext_x = c_ext[:, 0].max()-c_ext[:, 0].min()
shift = -c_ext.min(axis=0)             # cuerda axial desde x=0

# canales
canales = []
for sc in S_CH:
    xc = np.interp(sc, np.linspace(0, 1, 400), x)
    yc = np.interp(sc, np.linspace(0, 1, 400), y)
    hc = np.interp(sc, np.linspace(0, 1, 400), h)
    d = round(min(D_MAX if len(canales) < 8 else 3.1, 2*(hc-WALL_MIN)*0.96), 2)
    canales.append((xc, yc, d))
canales_xy = [(np.array([cx, cy])@R.T + shift, d) for cx, cy, d in canales]

# contorno exterior remuestreado por longitud de arco (sentido antihorario)
ext = Polygon((np.array(poly.exterior.coords)@R.T + shift))
n_pts = 240
dist = np.linspace(0, ext.exterior.length, n_pts, endpoint=False)
outer = np.array([ext.exterior.interpolate(d).coords[0] for d in dist])
ext_resampled = Polygon(outer)

# verificacion
print(f"Escalonamiento = {np.degrees(GAMMA):.2f} deg; giro de la camber = {np.degrees(PHI):.1f} deg")
print(f"Extension axial del contorno = {ext_x:.2f} mm (la cuerda axial C*cos(g) = {C*np.cos(GAMMA):.2f} mm entre las puntas LE y TE)")
xs, ys = ext.exterior.xy
print(f"Cuerda real aprox. (LE-TE) = {np.hypot(*(np.array([xs[np.argmin(xs)],ys[np.argmin(xs)]])-np.array([xs[np.argmax(xs)],ys[np.argmax(xs)]]))):.1f} mm")
print(f"Area del perfil = {ext.area:.1f} mm2; area remuestreada = {ext_resampled.area:.1f}")
ok = True
discs = [Point(*c).buffer(d/2, 64) for c, d in canales_xy]
for i, (c, d) in enumerate(canales_xy, 1):
    wall = ext.exterior.distance(Point(*c)) - d/2
    print(f"Canal {i:2d}: d = {d:5.2f} mm  pared minima = {wall:5.2f} mm")
    ok &= wall >= WALL_MIN - 1e-6
for i in range(len(discs)-1):
    lig = discs[i].distance(discs[i+1])
    print(f"  ligamento {i+1}-{i+2}: {lig:5.2f} mm")
    ok &= lig >= 1.5
print("Area de canales =", round(sum(p.area for p in discs), 1), "mm2")
print("RESULTADO VERIFICACION:", "OK" if ok else "REVISAR")

with open("perfil_exterior.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["x_mm", "y_mm"])
    for p in outer: w.writerow([f"{p[0]:.4f}", f"{p[1]:.4f}"])
with open("canales.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["canal", "x_mm", "y_mm", "d_mm"])
    for i, (c, d) in enumerate(canales_xy, 1):
        w.writerow([i, f"{c[0]:.4f}", f"{c[1]:.4f}", f"{d:.2f}"])

fig, ax = plt.subplots(figsize=(8, 8), dpi=120)
ax.fill(*ext.exterior.xy, color="#c9d3dc", ec="#14212f", lw=1.8)
for (c, d) in canales_xy:
    ax.add_patch(plt.Circle(c, d/2, fc="white", ec="#0f8a7c", lw=1.6))
shifted = ext.exterior.xy
ax.fill(shifted[0], np.array(shifted[1]) + PITCH, color="#e9eef2", ec="#8a97a5", lw=1, alpha=.6)
ax.set_aspect("equal"); ax.set_xlabel("x [mm] (axial)"); ax.set_ylabel("y [mm] (tangencial)")
ax.set_title("Alabe tipo C3X: perfil parametrico y 10 canales de enfriamiento")
ax.grid(alpha=.25); fig.tight_layout(); fig.savefig("previa_geometria.png")
