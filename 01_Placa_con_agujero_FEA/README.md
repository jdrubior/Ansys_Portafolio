# Proyecto 1: Placa con agujero a tensión (Ansys Mechanical)

Esfuerzo máximo en una placa de acero con agujero circular central bajo tensión uniaxial, con estudio de convergencia de malla y validación analítica.

## Modelo

| Parámetro | Valor |
|-----------|-------|
| Largo × ancho × espesor | 200 × 100 × 5 mm |
| Agujero | Ø20 mm (d/W = 0.2) |
| Material | Acero estructural (E = 200 GPa, ν = 0.3) |
| Carga | σ = 10 MPa en los extremos |
| Modelo | 1/4 de placa por doble simetría, esfuerzo plano 2D, cuadriláteros cuadráticos |

## Convergencia de malla

| Nodos | Elementos | σx,max (MPa) |
|------:|----------:|-------------:|
| 1 035 | 320 | 31.007 |
| 3 829 | 1 228 | 31.310 |
| 15 144 | 4 955 | 31.453 |
| 30 838 | 10 147 | 31.486 |
| 59 977 | 19 806 | 31.523 |

Extrapolación de Richardson: 31.59 MPa, GCI de la malla fina: 0.27 %.

## Validación

- Heywood (placa de ancho finito): σ_max ≈ 31.4 MPa, diferencia de 0.6 % con la malla fina.
- El perfil σx sobre la sección neta reproduce la forma de la solución de Kirsch.
- Equilibrio de la sección neta: 2500.03 N contra 2500 N aplicados (0.001 %).

## Contenido

- `informe/informe_placa_agujero.pdf`: informe técnico
- `informe/instructivo_placa_agujero.pdf`: guía ilustrada paso a paso
- `validacion/referencia_kirsch.py`: solución analítica, gráficas, Richardson/GCI y equilibrio
- `resultados/`: CSV de convergencia y del perfil de esfuerzos, imágenes de Ansys
