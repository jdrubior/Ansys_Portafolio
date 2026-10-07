"""
Estimación analítica de la resistencia térmica del disipador de aletas.
Ing. Juan Diego Rubio Ruiz

Modelo: convección forzada laminar sobre placa plana (Nu = 0.664 Re^1/2 Pr^1/3)
en cada aleta + eficiencia de aleta recta con longitud corregida.
Se evalúan dos velocidades características en el canal:
  (A) solo bloqueo de las aletas en el ancho:      U_c = U W / (W - N t_f)
  (B) área real del ducto modelado (aletas + base): U_c = U A_entrada / A_libre
Ambas definen una banda de referencia para validar el resultado de Fluent.
"""
import sys
sys.stdout.reconfigure(encoding="utf-8")  # consola de Windows (cp1252)
import numpy as np

# Geometría [m]
Lb, Wb, tb = 0.050, 0.050, 0.005          # base
N, tf, H, Lf = 10, 0.0015, 0.025, 0.050    # aletas
k_al = 205.0                               # W/(m K)
Q, T_in, U = 10.0, 25.0, 2.0               # W, °C, m/s

# Aire a ~30 °C
rho, mu, cp, k_air = 1.164, 1.872e-5, 1007.0, 0.02588
Pr = mu * cp / k_air


def resistencia_termica(U_c):
    """Devuelve un diccionario con Re, h, eficiencias y R_th para la velocidad de canal U_c."""
    Re = rho * U_c * Lf / mu
    Nu = 0.664 * Re**0.5 * Pr**(1 / 3)
    h = Nu * k_air / Lf
    m = np.sqrt(2 * h / (k_al * tf))
    Hc = H + tf / 2                        # longitud corregida (punta convectiva)
    eta_f = np.tanh(m * Hc) / (m * Hc)
    A_f = N * 2 * Hc * Lf                  # área de aletas
    A_b = Lb * Wb - N * tf * Lf            # base expuesta
    eta_o = 1 - A_f / (A_f + A_b) * (1 - eta_f)
    R_conv = 1 / (eta_o * h * (A_f + A_b))
    R_base = tb / (k_al * Lb * Wb)
    return dict(U_c=U_c, Re=Re, Nu=Nu, h=h, m=m, eta_f=eta_f, eta_o=eta_o,
                R_conv=R_conv, R_base=R_base, R_th=R_conv + R_base)


s = (Wb - N * tf) / (N - 1)
casos = {
    "A (bloqueo de aletas)": U * Wb / (Wb - N * tf),
    "B (ducto modelado)":    U * (Wb / 2 * (tb + H)) / ((Wb / 2 - N / 2 * tf) * H),
}

print(f"Pr = {Pr:.3f}   separación entre aletas s = {s*1e3:.2f} mm")
delta = None
for nombre, Uc in casos.items():
    r = resistencia_termica(Uc)
    delta = 5 * Lf / np.sqrt(r["Re"])
    print(f"\nCaso {nombre}")
    print(f"  U_canal = {r['U_c']:.2f} m/s   Re_L = {r['Re']:.0f}   h = {r['h']:.1f} W/m²K")
    print(f"  η_f = {r['eta_f']:.3f}   η_o = {r['eta_o']:.3f}")
    print(f"  R_th = {r['R_conv']:.3f} + {r['R_base']:.4f} = {r['R_th']:.3f} K/W"
          f"   ->  T_base ≈ {T_in + Q * r['R_th']:.1f} °C")
    print(f"  capa límite al final de la aleta δ ≈ 5L/√Re = {delta*1e3:.1f} mm  (s/2 = {s/2*1e3:.2f} mm)")

print("\nEn ambos casos δ > s/2: las capas límite de aletas vecinas se juntan y la correlación"
      " de placa plana es solo aproximada; los dos casos forman una banda de referencia.")
