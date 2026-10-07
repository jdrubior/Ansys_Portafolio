"""Dominio periodico de un pasaje para el alabe tipo C3X (corte de 12.7 mm).
Cuerpos: ALABE (solido), GAS (fluido externo) y REFRIGERANTE_1..10 (fluido interno).
Fronteras periodicas c1 y c2 = c1 + (0, PASO). El alabe queda completo dentro del dominio.
Salida: dominio_alabe.step + dominio_alabe.png. Verifica con shapely y OpenCascade.
"""
import csv, numpy as np
from shapely.geometry import Polygon, LineString, Point
from shapely import affinity
import cadquery as cq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PASO, H = 117.7, 12.7
X_IN_OFF, X_OUT_OFF = 60.0, 90.0        # distancias aguas arriba / abajo del alabe

outline = np.array([[float(a), float(b)] for a, b in list(csv.reader(open("perfil_exterior.csv")))[1:]])
canales = [(float(r[1]), float(r[2]), float(r[3])) for r in list(csv.reader(open("canales.csv")))[1:]]
A = Polygon(outline)
u0 = np.array([np.cos(np.radians(57.36)), -np.sin(np.radians(57.36))])        # direccion de la cuerda
proy = outline@u0
i_le, i_te = int(np.argmin(outline[:, 0])), int(np.argmax(outline[:, 0]))   # extremos axiales (x monotona en cada arco)
n = len(outline)
def arco(i0, i1):
    idx = [i0]
    while idx[-1] != i1: idx.append((idx[-1] + 1) % n)
    return outline[idx]
a1, a2 = arco(i_le, i_te), arco(i_te, i_le)[::-1]          # ambos de LE a TE
lado_inf = a1 if a1[:, 1].mean() < a2[:, 1].mean() else a2   # arco convexo (inferior)
lado_sup = a2 if lado_inf is a1 else a1
print("LE =", outline[i_le].round(1), " TE =", outline[i_te].round(1))

def construir_c1(desp=0.0):
    """Linea media del pasaje inferior: promedio entre el arco inferior del alabe y el arco
    superior del alabe vecino de abajo (desplazado -PASO). desp mueve la curva (mm)."""
    xs = np.linspace(max(lado_inf[0, 0], lado_sup[0, 0]), min(lado_inf[-1, 0], lado_sup[-1, 0]), 70)
    yi = np.interp(xs, lado_inf[:, 0], lado_inf[:, 1])
    ys = np.interp(xs, lado_sup[:, 0], lado_sup[:, 1])
    return np.c_[xs, 0.5*(yi + ys - PASO) + desp]
def evaluar(desp):
    c1 = construir_c1(desp)
    x_in = outline[:, 0].min() - X_IN_OFF; x_out = outline[:, 0].max() + X_OUT_OFF
    pts = [(x_in, c1[0, 1])] + [tuple(p) for p in c1] + [(x_out, c1[-1, 1])]
    c2 = [(x, y + PASO) for x, y in pts]
    dom = Polygon(pts + c2[::-1]).buffer(0)
    vecinos = [affinity.translate(A, 0, PASO), affinity.translate(A, 0, -PASO)]
    ok = dom.is_valid and dom.contains(A) and not any(dom.intersects(v) for v in vecinos)
    holg = min(LineString(c1).distance(A), LineString(c2).distance(A))
    return ok, holg, dom, pts, c2

mejor = None
for d in np.linspace(-6, 6, 25):
    ok, hl, *_ = evaluar(d)
    if ok and (mejor is None or hl > mejor[1]): mejor = (d, hl)
assert mejor, "no se encontro una curva periodica valida"
desp, holg = mejor
ok, _, dom, pts, c2 = evaluar(desp)
x_in, x_out = pts[0][0], pts[-1][0]
# garganta: distancia minima entre el alabe y su vecino superior
garganta = A.distance(affinity.translate(A, 0, PASO))
print(f"desplazamiento de c1 = {desp:.1f} mm, holgura minima a c1/c2 = {holg:.1f} mm")
print(f"garganta (alabe-vecino) = {garganta:.1f} mm; paso = {PASO} mm; garganta/paso = {garganta/PASO:.3f}")
print(f"dominio: x[{x_in:.1f},{x_out:.1f}], area XY = {dom.area:.0f} mm2, valido = {dom.is_valid}")

# ---------------- CadQuery ----------------
def spl(wp, p, inc=True):
    return wp.spline([(float(x), float(y)) for x, y in p], includeCurrent=inc)

vane = cq.Workplane("XY").spline([(float(x), float(y)) for x, y in outline], periodic=True).close().extrude(H)
solido = vane
for (cx, cy, d) in canales:
    solido = solido.cut(cq.Workplane("XY").center(cx, cy).circle(d/2).extrude(H))
refrig = [cq.Workplane("XY").center(cx, cy).circle(d/2).extrude(H) for (cx, cy, d) in canales]

c1 = construir_c1(desp)
wp = cq.Workplane("XY").moveTo(x_in, float(c1[0, 1])).lineTo(float(c1[0, 0]), float(c1[0, 1]))
wp = spl(wp, c1[1:])
wp = wp.lineTo(x_out, float(c1[-1, 1])).lineTo(x_out, float(c1[-1, 1] + PASO)).lineTo(float(c1[-1, 0]), float(c1[-1, 1] + PASO))
c2a = c1[::-1][1:].copy(); c2a[:, 1] += PASO
wp = spl(wp, c2a).lineTo(x_in, float(c1[0, 1] + PASO)).close()
dominio = wp.extrude(H)
gas = dominio.cut(vane)

def area_base(wp, tol=0.01):
    """Area de las bases planas por teselado fino (Area()/Volume() de OCC son poco fiables con splines)."""
    tot = 0.0
    for f in wp.faces("<Z").vals():
        v, t = f.tessellate(tol, 0.1)
        V = np.array([(p.x, p.y) for p in v]); T = np.array(t)
        a, b = V[T[:, 1]] - V[T[:, 0]], V[T[:, 2]] - V[T[:, 0]]
        tot += 0.5*np.abs(a[:, 0]*b[:, 1] - a[:, 1]*b[:, 0]).sum()
    return tot
ar = dict(dom=area_base(dominio), gas=area_base(gas), solido=area_base(solido), refrig=sum(area_base(r) for r in refrig))
print({k: round(v, 1) for k, v in ar.items()}, "mm2 (base)")
bal = ar["gas"] + ar["solido"] + ar["refrig"]
print(f"cierre de area: gas + solido + refrigerante = {bal:.1f} vs dominio = {ar['dom']:.1f} (dif {abs(bal-ar['dom'])/ar['dom']*100:.4f} %)")
print(f"volumenes (area x {H} mm): dominio {ar['dom']*H:.0f}, gas {ar['gas']*H:.0f}, alabe {ar['solido']*H:.0f}, refrigerante {ar['refrig']*H:.0f} mm3")
print("validos:", gas.val().isValid(), solido.val().isValid(), dominio.val().isValid())
print("solidos por cuerpo (gas, alabe):", len(gas.solids().vals()), len(solido.solids().vals()))

asm = cq.Assembly(name="dominio_alabe")
asm.add(solido, name="ALABE", color=cq.Color(0.55, 0.58, 0.62))
asm.add(gas, name="GAS", color=cq.Color(0.85, 0.5, 0.2, 0.4))
for i, r in enumerate(refrig, 1):
    asm.add(r, name=f"REFRIGERANTE_{i}", color=cq.Color(0.2, 0.5, 0.9))
asm.save("dominio_alabe.step")

fig, ax = plt.subplots(figsize=(11, 9), dpi=100)
ax.fill(*dom.exterior.xy, color="#fde7d2", ec="#d9822b", lw=1.5, label="dominio de gas (periodico)")
ax.fill(*A.exterior.xy, color="#8d97a3", ec="#14212f", lw=1.5, label="alabe (solido)")
for (cx, cy, d) in canales: ax.add_patch(plt.Circle((cx, cy), d/2, fc="#9cc7f3", ec="#1f6fc5"))
for v, ls in [(affinity.translate(A, 0, PASO), ":"), (affinity.translate(A, 0, -PASO), ":")]:
    ax.plot(*v.exterior.xy, "k", lw=0.8, ls=ls, alpha=.5)
ax.plot(*zip(*pts), color="#c0392b", lw=2, label="periodica c1")
ax.plot(*zip(*c2), color="#c0392b", lw=2, ls="--", label="periodica c2 = c1 + paso")
ax.set_aspect("equal"); ax.legend(loc="upper right", fontsize=9); ax.grid(alpha=.25)
ax.set_title("Dominio periodico del pasaje (corte de 12.7 mm, vista XY)")
fig.tight_layout(); fig.savefig("dominio_alabe.png")

