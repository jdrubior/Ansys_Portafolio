# Disipador de aletas - medio modelo (5 aletas) para Ansys SpaceClaim
# Ing. Juan Diego Rubio Ruiz
# Ejes: X = direccion del flujo (largo de las aletas), Y = ancho, Z = altura
# Uso: SpaceClaim -> Disenar -> Script -> pegar este codigo -> Ejecutar (boton play)

L    = 50.0                      # largo de base y aletas [mm]
W2   = 25.0                      # medio ancho de la base (medio modelo) [mm]
tb   = 5.0                       # espesor de la base [mm]
H    = 25.0                      # altura de las aletas [mm]
tf   = 1.5                       # espesor de las aletas [mm]
N    = 10                        # aletas del disipador completo
s    = (50.0 - N*tf)/(N - 1)     # separacion entre aletas = 3.889 mm
paso = tf + s                    # paso entre aletas = 5.389 mm

# Base
BlockBody.Create(Point.Create(MM(0), MM(0), MM(0)),
                 Point.Create(MM(L), MM(W2), MM(tb)),
                 ExtrudeType.ForceIndependent)

# 5 aletas; la quinta termina en y = 23.056 mm y deja media separacion (1.944 mm) hasta el plano de simetria y = 25
for i in range(5):
    y0 = i*paso
    BlockBody.Create(Point.Create(MM(0), MM(y0), MM(tb)),
                     Point.Create(MM(L), MM(y0 + tf), MM(tb + H)),
                     ExtrudeType.ForceIndependent)

# Unir base y aletas en un solo solido.
# Nota: en 2026 R1 la llamada con todos los cuerpos no los une; se hizo a mano con
# Design -> Combine (clic en la base, Ctrl+clic en cada aleta). El grabador registro:
#   targets = BodySelection.Create(Body1, Body2); Combine.Merge(targets, Info1)
