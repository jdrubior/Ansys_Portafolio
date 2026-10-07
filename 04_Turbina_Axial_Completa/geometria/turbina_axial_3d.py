"""Turbina axial de una etapa (estator + rotor + carcasa + eje) en 3D, exportada a STEP.
Eje de giro: Z. Unidades: mm. Verificada con OpenCascade (validez, volumenes, holgura).
Dos cuerpos: ROTOR (eje + disco + alabes torcidos) y ESTATOR (carcasa + aro + alabes fijos).
"""
import math, sys, time
import numpy as np
import cadquery as cq

N_ROTOR, N_ESTATOR = 24, 16      # relacion 3:2 -> sector periodico de 45 grados (3 alabes de rotor + 2 de estator)
R_EJE, R_BARRENO_ESTATOR = 20.0, 24.0
R_HUB, R_TIP, HOLGURA_PUNTA = 60.0, 110.0, 2.0
CUERDAS_ROTOR = [30.0, 26.0, 22.0]        # cuerda en raiz, medio y punta (solidez razonable)
R_CARC_IN, R_CARC_OUT = R_TIP + HOLGURA_PUNTA, R_TIP + HOLGURA_PUNTA + 10.0
Z_EST0, Z_EST1 = 20.0, 56.0            # estator
Z_ROT0, Z_ROT1 = 66.0, 106.0           # disco del rotor (espesor axial)
Z_CARC0, Z_CARC1 = 10.0, 150.0

def naca_perfil(c, tmax_frac, camber_frac, n=22):
    """Puntos (x,y) del perfil (borde de ataque en x=0, cuerda c), NACA 4 digitos con p=0.4."""
    x = (1 - np.cos(np.linspace(0, np.pi, n)))/2
    p = 0.4
    yt = 5*tmax_frac*(0.2969*np.sqrt(x) - 0.1260*x - 0.3516*x**2 + 0.2843*x**3 - 0.1036*x**4)
    yc = np.where(x < p, camber_frac/p**2*(2*p*x - x**2), camber_frac/(1-p)**2*((1-2*p) + 2*p*x - x**2))
    dyc = np.where(x < p, 2*camber_frac/p**2*(p - x), 2*camber_frac/(1-p)**2*(p - x))
    th = np.arctan(dyc)
    up = np.c_[x - yt*np.sin(th), yc + yt*np.cos(th)]*c
    lo = np.c_[x + yt*np.sin(th), yc - yt*np.cos(th)]*c
    lo[-1] = up[-1]                      # borde de salida cerrado
    return up, lo

def seccion(r, c, tfrac, camber, stagger_deg, zc):
    """Seccion plana en x=r. Coordenadas del plano YZ: y tangencial, z axial. Centrada en zc."""
    up, lo = naca_perfil(c, tfrac, camber)
    g = math.radians(stagger_deg)
    def tr(p):
        xs, ys = p[:, 0] - 0.5*c, p[:, 1]
        z = zc + (xs*math.cos(g) - ys*math.sin(g))
        y = (xs*math.sin(g) + ys*math.cos(g))
        return [(float(a), float(b)) for a, b in zip(y, z)]
    return tr(up), tr(lo)

def alabe(r_list, c_list, t_list, cam_list, st_list, zc, giro_deg=0.0):
    wp = None
    prev = None
    for i, r in enumerate(r_list):
        up, lo = seccion(r, c_list[i], t_list[i], cam_list[i], st_list[i], zc)
        wp = cq.Workplane("YZ", origin=(r, 0, 0)) if wp is None else wp.workplane(offset=r - prev, origin=(0, 0, 0))
        wp = wp.spline(up, includeCurrent=False).spline(lo[::-1][1:] + [up[0]], includeCurrent=True).close() if False else \
             wp.polyline(up + lo[::-1][1:-1]).close()
        prev = r
    return wp.loft(combine=True, ruled=False)

def rotar_z(shape, grados):
    return shape.rotate((0, 0, 0), (0, 0, 1), grados)

def construir():
    t0 = time.time()
    # ---------------- ROTOR ----------------
    eje = cq.Workplane("XY").workplane(offset=-40).circle(R_EJE).extrude(210)
    # disco con perfil revolucionado (alma delgada y aro grueso)
    perfil_disco = [(R_EJE, Z_ROT0 + 6), (30, Z_ROT0 + 6), (44, Z_ROT0 + 14), (R_HUB - 4, Z_ROT0), (R_HUB, Z_ROT0),
                    (R_HUB, Z_ROT1), (R_HUB - 4, Z_ROT1), (44, Z_ROT1 - 14), (30, Z_ROT1 - 6), (R_EJE, Z_ROT1 - 6)]
    disco = cq.Workplane("XZ").polyline(perfil_disco).close().revolve(360, (0, 0, 0), (0, 1, 0))
    rotor = eje.union(disco)
    zc = 0.5*(Z_ROT0 + Z_ROT1)
    b = alabe([R_HUB - 1.5, 0.5*(R_HUB + R_TIP), R_TIP],
              CUERDAS_ROTOR, [0.14, 0.12, 0.09], [0.07, 0.06, 0.05], [58.0, 44.0, 30.0], zc)
    alabes = [rotar_z(b, 360.0*k/N_ROTOR) for k in range(N_ROTOR)]
    for a in alabes:
        rotor = rotor.union(a, clean=False)
    rotor = rotor.clean()
    print(f"rotor listo ({time.time()-t0:.1f} s)")

    # ---------------- ESTATOR ----------------
    carcasa = (cq.Workplane("XY").workplane(offset=Z_CARC0).circle(R_CARC_OUT).circle(R_CARC_IN).extrude(Z_CARC1 - Z_CARC0))
    aro = (cq.Workplane("XY").workplane(offset=Z_EST0).circle(R_HUB).circle(R_BARRENO_ESTATOR).extrude(Z_EST1 - Z_EST0))
    estator = carcasa.union(aro)
    zce = 0.5*(Z_EST0 + Z_EST1)
    v = alabe([R_HUB - 1.5, 0.5*(R_HUB + R_CARC_IN), R_CARC_IN + 1.5],
              [34.0, 34.0, 34.0], [0.16, 0.16, 0.16], [0.10, 0.10, 0.10], [-35.0, -35.0, -35.0], zce)
    for k in range(N_ESTATOR):
        estator = estator.union(rotar_z(v, 360.0*k/N_ESTATOR), clean=False)
    estator = estator.clean()
    print(f"estator listo ({time.time()-t0:.1f} s)")
    return rotor, estator

if __name__ == "__main__":
    rotor, estator = construir()
    r_s, e_s = rotor.val(), estator.val()
    print("rotor valido:", r_s.isValid(), "solidos:", len(rotor.solids().vals()), "volumen mm3:", round(r_s.Volume()))
    print("estator valido:", e_s.isValid(), "solidos:", len(estator.solids().vals()), "volumen mm3:", round(e_s.Volume()))
    inter = r_s.intersect(e_s)
    print("volumen de interferencia rotor-estator mm3:", round(inter.Volume(), 3))
    bb = r_s.BoundingBox()
    print(f"bbox rotor: x[{bb.xmin:.1f},{bb.xmax:.1f}] z[{bb.zmin:.1f},{bb.zmax:.1f}]")
    asm = cq.Assembly(name="turbina_axial")
    asm.add(rotor, name="ROTOR", color=cq.Color(0.15, 0.4, 0.75))
    asm.add(estator, name="ESTATOR", color=cq.Color(0.7, 0.72, 0.75))
    asm.save("turbina_axial.step")
    cq.exporters.export(rotor, "rotor.step"); cq.exporters.export(estator, "estator.step")
    print("STEP exportados")
