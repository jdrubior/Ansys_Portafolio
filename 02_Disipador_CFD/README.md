# Proyecto 2: Disipador aletado con aire forzado (Ansys Fluent, CHT)

Resistencia térmica de un disipador de aluminio de 10 aletas rectas con aire forzado, mediante transferencia de calor conjugada (sólido y fluido acoplados), con independencia de malla y validación analítica.

## Modelo

| Parámetro | Valor |
|-----------|-------|
| Base | 50 × 50 × 5 mm, aluminio |
| Aletas | 10 aletas rectas, 1.5 mm de espesor, 25 mm de alto |
| Calor | 10 W en la cara inferior de la base |
| Aire | 25 °C, 2 m/s, flujo paralelo a las aletas |
| Física | Laminar estacionario con energía; medio modelo por simetría |

## Independencia de malla

| Celdas | T_max base (°C) | Δp (Pa) | R_th (K/W) |
|-------:|----------------:|--------:|-----------:|
| 121 992 | 37.262 | 9.615 | 1.226 |
| 273 504 | 37.268 | 8.958 | 1.227 |
| 608 180 | 37.315 | 8.890 | 1.232 |

GCI de la malla fina: 0.68 % en R_th y 1.37 % en Δp (orden de convergencia limitado a 2). Balances de masa y energía con error menor a 0.2 %.

## Validación

R_th = 1.23 K/W queda entre las dos estimaciones analíticas (placa plana laminar con eficiencia de aleta, dos velocidades de canal): 1.31 y 1.20 K/W, con diferencias de −6 % y +3 %.

## Contenido

- `informe/informe_disipador.pdf`: informe técnico
- `informe/instructivo_disipador.pdf`: guía ilustrada paso a paso
- `geometria/disipador_mitad_spaceclaim.py`: geometría paramétrica del medio modelo
- `validacion/referencia_disipador.py`: estimación analítica de R_th
- `validacion/analisis_independencia.py`: Richardson y GCI
- `resultados/`: CSV de mallas, imágenes y video de trayectorias
