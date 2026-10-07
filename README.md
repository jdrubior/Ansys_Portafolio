# Portafolio de simulación en Ansys

**Autor:** Ing. Juan Diego Rubio Ruiz
Proyectos propios de simulación (FEA y CFD) en Ansys 2026 R1 Student, con estudio de convergencia de malla, validación analítica en Python e informe técnico en PDF.

## Proyectos

| # | Proyecto | Herramienta | Resultado principal |
|---|----------|-------------|---------------------|
| 1 | [Placa con agujero a tensión](01_Placa_con_agujero_FEA/) | Ansys Mechanical | σx,max = 31.52 MPa (5 mallas), extrapolación de Richardson 31.59 MPa, GCI 0.27 %; a 0.6 % de la estimación de Heywood |
| 2 | [Disipador aletado con aire forzado](02_Disipador_CFD/) | Ansys Fluent (CHT) | R_th = 1.23 K/W (3 mallas), GCI 0.7 % en R_th y 1.4 % en Δp; entre las dos estimaciones analíticas (1.31 y 1.20 K/W) |

## Estructura de cada proyecto

- `geometria/`: scripts de geometría (cuando aplica)
- `resultados/`: tablas en CSV e imágenes exportadas de Ansys
- `validacion/`: scripts en Python con la solución analítica y el análisis de convergencia (NumPy, Matplotlib)
- `informe/`: informe técnico y instructivo paso a paso (PDF y código fuente LaTeX)

Los archivos de proyecto de Workbench y las soluciones de Fluent (`.wbpj`, `.h5`) no se incluyen por su tamaño.

## Reproducir la validación

```
pip install numpy matplotlib
python 01_Placa_con_agujero_FEA/validacion/referencia_kirsch.py
python 02_Disipador_CFD/validacion/referencia_disipador.py
python 02_Disipador_CFD/validacion/analisis_independencia.py
```

## Software

Ansys 2026 R1 Student (límites: 128 mil nodos/elementos en estructural y 1 millón de celdas en fluidos).

## Nota sobre el uso de IA

Los informes y los instructivos se redactaron con apoyo de un asistente de IA (Claude). Los modelos, el mallado, las simulaciones en Ansys, la interpretación de los resultados y la validación son trabajo propio.
