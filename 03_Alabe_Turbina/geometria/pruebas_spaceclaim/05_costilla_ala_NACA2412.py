# Python Script, API Version = V261
# -*- coding: utf-8 -*-
# Costilla de ala NACA 2412 con aligeramientos y largueros
# Contorno aerodinamico con borde de salida afilado y 5 agujeros.
# Unidades: mm. Plano XY, extrusion de 3.0 mm en +Z.
# Generado por generar_pruebas.py (verificado con shapely, NO probado en SpaceClaim).
# Si la version de API no coincide, graba un guion cualquiera y copia su 1a linea.
from System.Collections.Generic import List

ClearAll()
DIST = 3.0
SUAVE = False        # True: intenta spline en lugar de poligonal para los contornos

CONTORNOS = [      # contornos cerrados (exterior primero), lista de (x, y)
    [(200.000, 0.000), (175.926, -1.785), (148.734, -3.666), (122.533, -5.367), (101.638, -6.589), (84.156, -7.434), (67.284, -8.048), (54.516, -8.374), (42.630, -8.472), (37.087, -8.415), (31.844, -8.276), (24.599, -7.905), (20.215, -7.543), (16.213, -7.085), (12.614, -6.531), (10.971, -6.218), (8.009, -5.519), (5.492, -4.724), (3.435, -3.834), (2.582, -3.354), (1.848, -2.852), (1.235, -2.326), (0.742, -1.778), (0.372, -1.208), (0.124, -0.615), (0.000, 0.000), (0.000, 0.627), (0.126, 1.257), (0.378, 1.890), (0.756, 2.524), (1.259, 3.159), (1.887, 3.795), (3.516, 5.063), (5.637, 6.319), (8.241, 7.550), (9.720, 8.153), (13.024, 9.325), (16.777, 10.439), (20.959, 11.481), (25.551, 12.439), (30.528, 13.298), (35.864, 14.047), (41.533, 14.675), (47.503, 15.174), (53.743, 15.537), (60.220, 15.760), (66.898, 15.840), (73.743, 15.777), (80.712, 15.572), (87.725, 15.251), (94.797, 14.833), (101.891, 14.326), (112.498, 13.415), (119.492, 12.722), (129.783, 11.575), (139.730, 10.331), (149.221, 9.023), (158.151, 7.686), (166.422, 6.356), (173.944, 5.068), (182.666, 3.479), (191.239, 1.814)]
]
CIRCULOS = [       # (x, y, radio)
    (50.000, 3.437, 9.000),
    (120.000, 3.556, 5.500),
    (28.000, 2.310, 6.500),
    (80.000, 4.000, 9.090),
    (156.000, 2.396, 3.060)
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
