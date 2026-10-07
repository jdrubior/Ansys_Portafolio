# Proyecto 3 · Álabe tipo C3X: plan de simulación (Ansys Student)

Geometría: `geometria/dominio_alabe.step` (corte de 12.7 mm). Referencia: `validacion/referencia_alabe.py`.

## Decisión de modelado (por qué el refrigerante va como convección)
El flujo de los canales es radial (a lo largo de la envergadura). Un corte de 12.7 mm no reproduce un tramo desarrollado
ni el calentamiento del aire, y la condición periódica en Z no admite aporte de calor en estado estacionario. Por eso:

- **Modelo base (validable):** Fluent conjugado, gas + sólido. En las paredes de los 10 canales se impone convección:
  `h_c` de Dittus-Boelter y `T_c = 332 K` (media envergadura), valores en `resultados/referencia_analitica.csv`.
  En Workbench se **suprimen** los 10 cuerpos REFRIGERANTE_n.
- **Ampliación opcional:** canal completo con CFD en un tubo aparte (2D axisimétrico) para contrastar `h_c` con la correlación.

## Fluent
- Solver acoplado por densidad (flujo transónico, Ma salida 0.9), gas ideal, viscosidad de Sutherland.
- Turbulencia k-omega SST con y+ ~ 1: primera celda 2.3 um, 32 capas con razón 1.2.
- Entrada: presión total 203.0 kPa, Tt 689 K, intensidad de turbulencia 6.5 %. Salida: presión estática 120.0 kPa.
- Periodicidad traslacional c1/c2 con desplazamiento (0, 117.7 mm). Caras en Z: simetría (gas) y adiabática/simetría (sólido).
- Sólido: acero inoxidable tipo 310, k = 19 W/mK.
- Convergencia: residuales < 1e-5 y balance de energía < 0.5 %.

## Malla (límite Student: 1 M de celdas)
| malla | h_bulk | paso en pared | capas en Z | celdas estimadas |
|---|---|---|---|---|
| gruesa | 2.4 mm | 1.0 mm | 6 | ~108 mil |
| media | 1.6 mm | 0.7 mm | 9 | ~275 mil |
| fina | 1.1 mm | 0.5 mm | 13 | ~680 mil |

Razón de refinamiento ~1.35 entre mallas, adecuada para GCI. Quedan unas 320 mil celdas de margen bajo el límite.
Estas cifras son estimaciones: se confirman con el conteo real de Fluent Meshing/Mechanical.

## Métricas a comparar y reportar
Temperatura media y máxima del metal, calor total al refrigerante, h_g local (comparar con 700 a 1500 W/m2K),
caída de presión total del pasaje, GCI de T_metal y de q. Contraste con el balance analítico (T_metal isotermo 609 K como cota).

## Etapa multifísica (Mechanical)
Importar la temperatura del sólido a Static Structural (Workbench, "Imported Temperature"), esfuerzos térmicos y von Mises;
límite Student: 128 mil nodos (malla del sólido, ~3,700 mm2 x 12.7 mm, sobra).
