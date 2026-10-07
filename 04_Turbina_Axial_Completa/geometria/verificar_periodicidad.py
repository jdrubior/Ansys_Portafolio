"""Comprueba que las caras periodicas (planos 0 y 45 grados) del dominio de gas coinciden al rotar 45 grados."""
import math, numpy as np, cadquery as cq
from shapely.geometry import Polygon
from shapely.ops import unary_union
g = cq.importers.importStep("dominio_gas_sector.step").val()
def plano(ang):
    n = np.array([-math.sin(math.radians(ang)), math.cos(math.radians(ang)), 0.0])
    tot = 0.0; polys = []
    for f in g.Faces():
        c = f.Center(); nrm = f.normalAt(c)
        d = abs(np.dot([c.x, c.y, c.z], n))
        if d < 1e-3 and abs(abs(np.dot([nrm.x, nrm.y, nrm.z], n)) - 1) < 1e-4:
            v, t = f.tessellate(0.05, 0.1)
            tot += sum(0.5*np.linalg.norm(np.cross(np.array([v[a].x, v[a].y, v[a].z]) - np.array([v[b].x, v[b].y, v[b].z]),
                                                  np.array([v[c2].x, v[c2].y, v[c2].z]) - np.array([v[b].x, v[b].y, v[b].z])))
                       for a, b, c2 in t)
    return tot
a0, a45 = plano(0), plano(45)
print(f"area cara periodica 0 deg: {a0:.1f} mm2 | 45 deg: {a45:.1f} mm2 | diferencia {abs(a0-a45)/a0*100:.3f} %")
