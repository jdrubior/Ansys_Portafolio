"""Sector periodico de 45 grados (1/8) de la turbina axial: 2 alabes de estator + 3 de rotor.
Genera:
  dominio_gas_sector.step : gas (anillo entre aro y carcasa) menos alabes -> para Fluent (MRF)
  sector_rotor_solido.step: disco + eje + 3 alabes en 30 grados -> para Mechanical (simetria ciclica)
Verifica validez, interferencias y cierre de volumen por teselado (el Volume() de OCC no es fiable con splines).
"""
import math, time
import numpy as np
import cadquery as cq
from turbina_axial_3d import (alabe, rotar_z, R_EJE, R_HUB, R_TIP, R_CARC_IN, Z_EST0, Z_EST1, Z_ROT0, Z_ROT1,
                              Z_CARC0, Z_CARC1, N_ROTOR, N_ESTATOR, CUERDAS_ROTOR)
SECTOR = 360.0/8             # 45 grados: 16/8 = 2 alabes de estator, 24/8 = 3 de rotor
Z_IN, Z_OUT = Z_CARC0, Z_CARC1

def vol_teselado(wp, tol=0.05):
    sh = wp.val() if hasattr(wp, "val") else wp
    v, t = sh.tessellate(tol, 0.2)
    if len(t) == 0:
        return 0.0
    V = np.array([(p.x, p.y, p.z) for p in v]); T = np.array(t)
    a, b, c = V[T[:, 0]], V[T[:, 1]], V[T[:, 2]]
    return abs(np.einsum("ij,ij->i", a, np.cross(b, c)).sum()/6.0)

def cuña(r0, r1, z0, z1, ang):
    """Anillo r0..r1, z0..z1 barrido 'ang' grados desde el plano XZ (verificado abajo)."""
    return (cq.Workplane("XZ").polyline([(r0, z0), (r1, z0), (r1, z1), (r0, z1)]).close()
            .revolve(ang, (0, 0, 0), (0, 1, 0)))

t0 = time.time()
gas0 = cuña(R_HUB, R_CARC_IN, Z_IN, Z_OUT, SECTOR)
bb = gas0.val().BoundingBox()
signo = 1.0 if abs(bb.ymax) > abs(bb.ymin) else -1.0
print(f"sector barrido hacia {'+Y' if signo > 0 else '-Y'}: y[{bb.ymin:.1f},{bb.ymax:.1f}] x[{bb.xmin:.1f},{bb.xmax:.1f}]")
if signo < 0:
    gas0 = gas0.mirror("XZ")                  # asegura el sector entre 0 y +30 grados

zc_r, zc_e = 0.5*(Z_ROT0 + Z_ROT1), 0.5*(Z_EST0 + Z_EST1)
b_rot = alabe([R_HUB - 1.5, 0.5*(R_HUB + R_TIP), R_TIP], CUERDAS_ROTOR, [0.14, 0.12, 0.09],
              [0.07, 0.06, 0.05], [58.0, 44.0, 30.0], zc_r)
b_est = alabe([R_HUB - 1.5, 0.5*(R_HUB + R_CARC_IN), R_CARC_IN + 1.5], [34.0]*3, [0.16]*3, [0.10]*3, [-35.0]*3, zc_e)
# Se incluyen alabes vecinos (k fuera del sector) porque los alabes muy escalonados invaden los planos periodicos:
# las porciones que sobresalen de un plano deben reaparecer en el opuesto para que el patron sea periodico.
ang_r = [(k + 0.5)*360.0/N_ROTOR for k in range(-2, 6)]
ang_e = [(k + 0.5)*360.0/N_ESTATOR for k in range(-1, 4)]
rot_blades = [rotar_z(b_rot, a) for a in ang_r]
est_vanes = [rotar_z(b_est, a) for a in ang_e]

gas = gas0
for s in rot_blades + est_vanes:
    gas = gas.cut(s)
print(f"gas listo ({time.time()-t0:.1f} s) | validez: {gas.val().isValid()} | solidos: {len(gas.solids().vals())}")

# cierre de volumen: gas + (alabes dentro del anillo) = anillo
dentro = []
for s in rot_blades + est_vanes:
    try:
        dentro.append(vol_teselado(gas0.intersect(s)))
    except Exception:
        dentro.append(0.0)       # alabe vecino que no toca el sector
v_g0, v_g = vol_teselado(gas0), vol_teselado(gas)
print(f"anillo {v_g0:.0f} mm3 | gas {v_g:.0f} | alabes dentro {sum(dentro):.0f} | cierre: {abs(v_g0 - v_g - sum(dentro))/v_g0*100:.3f} %")
print(f"bloqueo por alabes = {sum(dentro)/v_g0*100:.1f} % del volumen del anillo")

# sector solido del rotor para Mechanical
from turbina_axial_3d import R_EJE
perfil_disco = [(R_EJE, Z_ROT0 + 6), (30, Z_ROT0 + 6), (44, Z_ROT0 + 14), (R_HUB - 4, Z_ROT0), (R_HUB, Z_ROT0),
                (R_HUB, Z_ROT1), (R_HUB - 4, Z_ROT1), (44, Z_ROT1 - 14), (30, Z_ROT1 - 6), (R_EJE, Z_ROT1 - 6)]
disco = cq.Workplane("XZ").polyline(perfil_disco).close().revolve(SECTOR, (0, 0, 0), (0, 1, 0))
eje = cuña(0.0001, R_EJE, -40.0, 170.0, SECTOR)
if signo < 0:
    disco, eje = disco.mirror("XZ"), eje.mirror("XZ")
sect = disco.union(eje)
for s in rot_blades:
    sect = sect.union(s, clean=False)
sect = sect.clean()
cuna_sec = cuña(0.0001, R_TIP + 6.0, -45.0, 175.0, SECTOR)
if signo < 0: cuna_sec = cuna_sec.mirror("XZ")
sect = sect.intersect(cuna_sec)          # corte limpio en los planos 0 y 45 grados
print(f"sector rotor listo ({time.time()-t0:.1f} s) | valido: {sect.val().isValid()} | solidos: {len(sect.solids().vals())}")
bbs = sect.val().BoundingBox()
print(f"sector rotor: y[{bbs.ymin:.1f},{bbs.ymax:.1f}], x[{bbs.xmin:.1f},{bbs.xmax:.1f}]")

cq.exporters.export(gas, "dominio_gas_sector.step")
cq.exporters.export(sect, "sector_rotor_solido.step")
print("STEP exportados")
