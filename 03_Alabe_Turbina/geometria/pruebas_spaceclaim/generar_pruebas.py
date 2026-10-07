"""Genera guiones de SpaceClaim con geometrias 2D complejas (para extruir) y las
verifica con shapely. Cada guion dibuja contornos con SketchLine/SketchCircle y
extruye la cara principal."""
import numpy as np
from shapely.geometry import Polygon, Point, LineString
from shapely.ops import unary_union
from shapely import affinity
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PLANTILLA = '''# Python Script, API Version = V261
# -*- coding: utf-8 -*-
# {titulo}
# {desc}
# Unidades: mm. Plano XY, extrusion de {dist} mm en +Z.
# Generado por generar_pruebas.py (verificado con shapely, NO probado en SpaceClaim).
# Si la version de API no coincide, graba un guion cualquiera y copia su 1a linea.
from System.Collections.Generic import List

ClearAll()
DIST = {dist}
SUAVE = False        # True: intenta spline en lugar de poligonal para los contornos

CONTORNOS = [      # contornos cerrados (exterior primero), lista de (x, y)
{contornos}
]
CIRCULOS = [       # (x, y, radio)
{circulos}
]

ViewHelper.SetSketchPlane(Plane.PlaneXY)

def curva_cerrada(pts):
    p2 = List[Point2D]()
    for (x, y) in pts:
        p2.Add(Point2D.Create(MM(x), MM(y)))
    p2.Add(Point2D.Create(MM(pts[0][0]), MM(pts[0][1])))
    if SUAVE:
        try:
            SketchNurbs.CreateFrom2DPoints(False, p2)
            return
        except Exception as e:
            print("Spline fallo (" + str(e) + "), se usa poligonal")
    for i in range(p2.Count - 1):
        SketchLine.Create(p2[i], p2[i + 1])

for c in CONTORNOS:
    curva_cerrada(c)
for (x, y, r) in CIRCULOS:
    SketchCircle.Create(Point2D.Create(MM(x), MM(y)), MM(r))
print("Contornos: " + str(len(CONTORNOS)) + ", circulos: " + str(len(CIRCULOS)))

ViewHelper.SetViewMode(InteractionMode.Solid)
caras = GetRootPart().Bodies[0].Faces
print("Caras creadas por el dibujo: " + str(caras.Count))

def area(f):
    for get in (lambda: f.Shape.Area, lambda: f.Area):
        try:
            return get()
        except:
            pass
    return 0

# extruir solo la cara principal (la de mayor area): la pieza con sus huecos
mejor = 0
amax = -1
for i in range(caras.Count):
    a = area(caras[i])
    if a > amax:
        amax = a
        mejor = i
opts = ExtrudeFaceOptions()
try:
    sel = Selection.Create(GetRootPart().Bodies[0].Faces[mejor])
    ExtrudeFaces.Execute(sel, MM(DIST), opts)
    print("Extrusion OK (cara " + str(mejor) + ")")
except Exception as e:
    print("Extrusion fallo: " + str(e))
print("Cuerpos: " + str(GetRootPart().Bodies.Count))
print("Nota: si quedaron discos sueltos de los huecos, borralos (son regiones del dibujo).")
'''

def anillo(poly, paso=None):
    return [(round(x, 4), round(y, 4)) for x, y in list(poly.exterior.coords)[:-1]]

def resample(poly, n):
    L = poly.exterior.length
    return [poly.exterior.interpolate(d).coords[0] for d in np.linspace(0, L, n, endpoint=False)]

def fmt(lst, f):
    return ",\n".join("    " + f(x) for x in lst)

def escribir(nombre, titulo, desc, dist, contornos, circulos):
    txt = PLANTILLA.format(
        titulo=titulo, desc=desc, dist=dist,
        contornos=fmt(contornos, lambda c: "[" + ", ".join("(%.3f, %.3f)" % p for p in c) + "]"),
        circulos=fmt(circulos, lambda c: "(%.3f, %.3f, %.3f)" % c))
    open(nombre + ".py", "w", encoding="utf-8").write(txt)

def verificar(nombre, ext, holes):
    ok = ext.is_valid and ext.exterior.is_simple
    discs = [Point(x, y).buffer(r, 48) for x, y, r in holes]
    minpared = min([ext.exterior.distance(Point(x, y)) - r for x, y, r in holes] or [9e9])
    minlig = 9e9
    for i in range(len(discs)):
        for j in range(i + 1, len(discs)):
            minlig = min(minlig, discs[i].distance(discs[j]))
    solap = any(discs[i].intersects(discs[j]) for i in range(len(discs)) for j in range(i+1, len(discs)))
    dentro = all(ext.contains(d) for d in discs)
    ok = ok and dentro and not solap
    print(f"{nombre:24s} valido={ext.is_valid} area={ext.area:9.1f} huecos={len(holes):3d} "
          f"pared_min={minpared:5.2f} ligamento_min={minlig:5.2f} -> {'OK' if ok else 'REVISAR'}")
    return ok

# ---------------- 1. Engranaje recto con aligeramientos ----------------
def inv(a): return np.tan(a) - a
def engranaje(m=3.0, z=24, alpha=np.radians(20)):
    r = m*z/2; ra = r + m; rf = r - 1.25*m; rb = r*np.cos(alpha)
    semi_p = np.pi*m/2/(2*r)
    def h(rho):
        phi = np.arccos(min(rb/rho, 1.0))
        return semi_p + inv(alpha) - inv(phi)
    rhos = np.linspace(rb, ra, 14)
    pts = []
    for k in range(z):
        c = 2*np.pi*k/z
        flanco_d = [(rho*np.cos(c - h(rho)), rho*np.sin(c - h(rho))) for rho in rhos]
        punta = [(ra*np.cos(c + t), ra*np.sin(c + t)) for t in np.linspace(-h(ra), h(ra), 4)[1:-1]]
        flanco_i = [(rho*np.cos(c + h(rho)), rho*np.sin(c + h(rho))) for rho in rhos[::-1]]
        raiz_i = (rf*np.cos(c + h(rb)), rf*np.sin(c + h(rb)))
        raiz_d = (rf*np.cos(c - h(rb)), rf*np.sin(c - h(rb)))
        pts += [raiz_d] + flanco_d + punta + flanco_i + [raiz_i]
    return Polygon(pts).buffer(0), rf, r

def hacer_engranaje():
    g, rf, r = engranaje()
    holes = [(rf*0.62*np.cos(a), rf*0.62*np.sin(a), 4.5) for a in np.linspace(0, 2*np.pi, 6, endpoint=False)]
    holes.append((0.0, 0.0, 8.0))
    return g, holes, 15.0

# ---------------- 2. Impulsor centrifugo de 7 alabes curvos ----------------
def hacer_impulsor():
    n = 7
    hub = Point(0, 0).buffer(26, 96)
    vanes = []
    for k in range(n):
        t = np.linspace(0, 1, 60)
        rr = 22 + 63*t
        th = 1.55*np.log(rr/22)*1.0 + 2*np.pi*k/n
        ls = LineString(np.c_[rr*np.cos(th), rr*np.sin(th)])
        vanes.append(ls.buffer(1.9, 16))
    forma = unary_union([hub] + vanes)
    forma = forma.buffer(1.0, 24).buffer(-1.0, 24)       # empalmes en la raiz del alabe
    ext = Polygon(forma.exterior)
    holes = [(0.0, 0.0, 10.0)]
    return ext, holes, 20.0

# ---------------- 3. Disco de freno perforado ----------------
def hacer_disco():
    ext = Point(0, 0).buffer(140, 180)
    holes = [(0.0, 0.0, 32.0)]
    holes += [(55*np.cos(a), 55*np.sin(a), 6.5) for a in np.linspace(0, 2*np.pi, 5, endpoint=False)]
    holes += [(118*np.cos(a), 118*np.sin(a), 4.0) for a in np.linspace(0, 2*np.pi, 36, endpoint=False)]
    holes += [(102*np.cos(a), 102*np.sin(a), 4.0) for a in np.linspace(0, 2*np.pi, 36, endpoint=False) + np.pi/36]
    holes += [(86*np.cos(a), 86*np.sin(a), 4.0) for a in np.linspace(0, 2*np.pi, 36, endpoint=False)]
    return ext, holes, 10.0

# ---------------- 4. Placa panal (patron hexagonal de agujeros) ----------------
def hacer_panal():
    R = 95
    ext = Point(0, 0).buffer(R, 180)
    holes = [(0.0, 0.0, 9.0)]
    p = 13.0
    for j in range(-9, 10):
        for i in range(-9, 10):
            x = p*(i + 0.5*(j % 2)); y = p*np.sqrt(3)/2*j
            rho = np.hypot(x, y)
            if 22 < rho < 76: holes.append((x, y, 4.6))
    holes += [(87*np.cos(a), 87*np.sin(a), 3.0) for a in np.linspace(0, 2*np.pi, 12, endpoint=False)]
    return ext, holes, 6.0

# ---------------- 5. Costilla de ala (NACA 2412) con aligeramientos ----------------
def naca4(m, p, t, c, n=90):
    x = (1 - np.cos(np.linspace(0, np.pi, n)))/2
    yt = 5*t*(0.2969*np.sqrt(x) - 0.1260*x - 0.3516*x**2 + 0.2843*x**3 - 0.1036*x**4)  # TE cerrado
    yc = np.where(x < p, m/p**2*(2*p*x - x**2), m/(1-p)**2*((1-2*p) + 2*p*x - x**2))
    dyc = np.where(x < p, 2*m/p**2*(p - x), 2*m/(1-p)**2*(p - x))
    th = np.arctan(dyc)
    xu = x - yt*np.sin(th); yu = yc + yt*np.cos(th)
    xl = x + yt*np.sin(th); yl = yc - yt*np.cos(th)
    pts = list(zip(xu[::-1]*c, yu[::-1]*c)) + list(zip(xl[1:]*c, yl[1:]*c))
    return pts
def hacer_costilla():
    c = 200.0
    ext = Polygon(naca4(0.02, 0.4, 0.12, c)).buffer(0)
    holes = [(0.25*c, 0.0, 9.0), (0.60*c, 0.0, 5.5)]       # largueros
    holes += [(0.14*c, 0.0, 6.5), (0.40*c, 0.0, 12.0), (0.78*c, 0.0, 3.5)]
    # ajustar la y de cada agujero a la linea media local
    ajust = []
    for x, _, r in holes:
        yc = 0.02/0.4**2*(2*0.4*(x/c) - (x/c)**2) if x/c < 0.4 else 0.02/0.6**2*((1-0.8) + 0.8*(x/c) - (x/c)**2)
        rmax = ext.exterior.distance(Point(x, yc*c)) - 2.5     # pared minima 2.5 mm
        r = min(r, rmax)
        if r >= 2.0:
            ajust.append((x, yc*c, round(r, 2)))
    return ext, ajust, 3.0


# ---------------- helpers de forma ----------------
def anillo_pts(ring, maxn=700):
    p = list(ring.coords)[:-1]
    if len(p) > maxn:
        p = [ring.interpolate(d).coords[0] for d in np.linspace(0, ring.length, maxn, endpoint=False)]
    return [(round(x, 4), round(y, 4)) for x, y in p]

def semi_annulus(cx, cy, R, t, a0, a1):
    a = np.linspace(a0, a1, 60)
    ext = [(cx + R*np.cos(x), cy + R*np.sin(x)) for x in a]
    inte = [(cx + (R - t)*np.cos(x), cy + (R - t)*np.sin(x)) for x in a[::-1]]
    return Polygon(ext + inte).buffer(0)

# ---------------- 6. Turbina axial con carenado (vista frontal, alabes radiales) ----------------
def hacer_turbina_axial():
    n = 31
    rh, rs = 62.0, 100.0
    hub = Point(0, 0).buffer(rh, 144)
    shroud = Point(0, 0).buffer(106, 200).difference(Point(0, 0).buffer(rs, 200))
    blades = []
    for k in range(n):
        s = np.linspace(0, 1, 40)
        th = 2*np.pi*k/n + 0.45*s                       # inclinacion (estela del alabe)
        r = (rh - 4) + (rs - rh + 8)*s
        cx, cy = r*np.cos(th), r*np.sin(th)
        tx, ty = np.gradient(cx), np.gradient(cy)
        nrm = np.hypot(tx, ty); nx_, ny_ = -ty/nrm, tx/nrm
        w = 1.2 + 3.4*np.sin(np.pi*np.clip(s*1.05, 0, 1))**0.6      # perfil de espesor
        left = list(zip(cx + w*nx_, cy + w*ny_)); right = list(zip(cx - w*nx_, cy - w*ny_))[::-1]
        blades.append(Polygon(left + right).buffer(0))
    forma = unary_union([hub, shroud] + blades).buffer(0.8, 24).buffer(-0.8, 24)
    holes = [(0.0, 0.0, 18.0)]
    holes += [(40*np.cos(a), 40*np.sin(a), 4.5) for a in np.linspace(0, 2*np.pi, 8, endpoint=False) + np.pi/8]
    return forma, holes, 25.0

# ---------------- 7. Rueda Pelton (20 cucharas) ----------------
def hacer_pelton():
    n = 20
    disco = Point(0, 0).buffer(68, 144)
    cucharas = []
    for k in range(n):
        a = 2*np.pi*k/n
        for off in (-6.5, 6.5):
            e = affinity.scale(Point(0, 0).buffer(1, 48), 15.5, 10.5)
            e = affinity.translate(e, 82, off)
            cucharas.append(affinity.rotate(e, a, origin=(0, 0), use_radians=True))
    forma = unary_union([disco] + cucharas).buffer(1.2, 24).buffer(-1.2, 24)
    holes = [(0.0, 0.0, 14.0)]
    holes += [(28*np.cos(a), 28*np.sin(a), 3.5) for a in np.linspace(0, 2*np.pi, 6, endpoint=False)]
    holes += [(48*np.cos(a), 48*np.sin(a), 8.0) for a in np.linspace(0, 2*np.pi, 10, endpoint=False) + np.pi/10]
    return forma, holes, 30.0

# ---------------- 8. Rotor Savonius (S con traslape) ----------------
def hacer_savonius():
    r, t, g = 50.0, 5.0, 8.0
    a = semi_annulus(g, r, r, t, -np.pi/2, np.pi/2)
    b = semi_annulus(-g, -r, r, t, np.pi/2, 3*np.pi/2)
    eje = Point(0, 0).buffer(12, 64)
    forma = unary_union([a, b, eje]).buffer(0.5, 24).buffer(-0.5, 24)
    holes = [(0.0, 0.0, 5.0)]
    return forma, holes, 200.0

# ---------------- 9. Rotor Darrieus tipo H (3 alabes NACA 0018) ----------------
def hacer_darrieus():
    R, c = 400.0, 80.0
    pts = naca4(0.0, 0.5, 0.18, c)
    perfil = Polygon(pts).buffer(0)
    perfil = affinity.translate(perfil, -0.25*c, 0)              # cuarto de cuerda en el origen
    piezas = [Point(0, 0).buffer(30, 96)]
    for k in range(3):
        a = 2*np.pi*k/3
        blade = affinity.translate(perfil, 0, R)                   # cuerda tangente
        blade = affinity.rotate(blade, a, origin=(0, 0), use_radians=True)
        piezas.append(blade)
        for dz in (-0.5, 0.5):   # dos puntales por alabe (a +-25% y 75% de la envergadura aparente)
            pass
        strut = LineString([(0, 0), (0, R)]).buffer(4.0, cap_style=2)
        piezas.append(affinity.rotate(strut, a, origin=(0, 0), use_radians=True))
    forma = unary_union(piezas).buffer(1.0, 24).buffer(-1.0, 24)
    holes = [(0.0, 0.0, 10.0)]
    return forma, holes, 50.0

CASOS = [
    ("01_engranaje_recto", "Engranaje recto m=3, z=24, 20 deg, con aligeramientos", "Perfil de evolvente (poligonal fina), 6 agujeros y barreno central.", hacer_engranaje),
    ("02_impulsor_centrifugo", "Impulsor centrifugo de 7 alabes curvos (espiral logaritmica)", "Cubo central y alabes con empalme; contorno unico complejo (muchos segmentos).", hacer_impulsor),
    ("03_disco_freno_perforado", "Disco de freno perforado (3 coronas de agujeros + 5 barrenos)", "Prueba de muchos circulos dentro de un contorno circular.", hacer_disco),
    ("04_placa_panal", "Placa con patron hexagonal de agujeros (panal)", "Mas de cien circulos: sirve para ver el rendimiento de SpaceClaim.", hacer_panal),
    ("05_costilla_ala_NACA2412", "Costilla de ala NACA 2412 con aligeramientos y largueros", "Contorno aerodinamico con borde de salida afilado y 5 agujeros.", hacer_costilla),
    ("06_turbina_axial_carenada", "Turbina axial con carenado: cubo, 31 alabes inclinados y aro exterior", "Cara con 31 pasajes interiores (contornos huecos) y 9 circulos.", hacer_turbina_axial),
    ("07_rueda_pelton", "Rueda Pelton de 20 cucharas dobles", "Contorno festoneado con aligeramientos y barrenos.", hacer_pelton),
    ("08_rotor_savonius", "Rotor Savonius (S con traslape) y eje", "Dos semianillos curvos unidos por el eje: pieza delgada y curva.", hacer_savonius),
    ("09_rotor_darrieus_H", "Rotor Darrieus tipo H: 3 alabes NACA 0018 con puntales", "Caso clasico de CFD 2D de turbinas eolicas de eje vertical.", hacer_darrieus),
]

if __name__ == "__main__":
    fig, axs = plt.subplots(3, 3, figsize=(15, 15), dpi=90)
    for ax, (nom, tit, desc, fn) in zip(axs.flat, CASOS):
        ext, holes, dist = fn()
        if ext.geom_type != "Polygon":
            ext = max(ext.geoms, key=lambda g: g.area)
        ext = ext.simplify(0.04, preserve_topology=True)      # menos segmentos en SpaceClaim
        okc = verificar(nom, ext, holes)
        rings = [ext.exterior] + list(ext.interiors)
        cont = [anillo_pts(rg) for rg in rings]
        escribir(nom, tit, desc, dist, cont, [(x, y, r) for x, y, r in holes])
        ax.fill(*ext.exterior.xy, color="#c9d3dc", ec="#14212f", lw=1.0)
        for ir in ext.interiors:
            ax.fill(*ir.xy, color="white", ec="#14212f", lw=0.8)
        for x, y, r in holes:
            ax.add_patch(plt.Circle((x, y), r, fc="white", ec="#0f8a7c", lw=0.9))
        npts = sum(len(c) for c in cont)
        ax.set_aspect("equal"); ax.set_title(f"{nom}\n{len(cont)} contornos, {npts} segmentos, {len(holes)} circulos", fontsize=10)
        ax.axis("off")
    fig.tight_layout(); fig.savefig("previa_geometrias.png")
