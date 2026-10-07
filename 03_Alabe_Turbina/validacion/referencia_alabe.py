"""Referencia analitica y presupuesto de malla del proyecto 3 (alabe tipo C3X, corte periodico).

Condiciones tomadas de la literatura del C3X (NASA CR-174827): Tt,in = 689 K, Ma salida = 0.90,
Re salida (cuerda real) = 1.5e6. El refrigerante es aire a 300 K (supuesto de este proyecto).
Las estimaciones de h_g y de temperatura de metal son de orden de magnitud (+-30 %): sirven para
verificar que la simulacion cae en el rango esperado, no para validar punto a punto.
"""
import csv
import numpy as np

gam, Rg, cp, Pr = 1.4, 287.0, 1004.5, 0.71
def mu_suth(T): return 1.458e-6*T**1.5/(T + 110.4)
def k_aire(T): return mu_suth(T)*cp/Pr

# ---------------- condiciones del gas ----------------
Tt, M2, Re2, C = 689.0, 0.90, 1.5e6, 0.14493
T2 = Tt/(1 + 0.5*(gam - 1)*M2**2)
a2 = np.sqrt(gam*Rg*T2); U2 = M2*a2; mu2 = mu_suth(T2)
rho2 = Re2*mu2/(U2*C); p2 = rho2*Rg*T2
pt = p2*(1 + 0.5*(gam - 1)*M2**2)**(gam/(gam - 1))
r = Pr**(1/3); Taw = T2*(1 + r*0.5*(gam - 1)*M2**2)
print(f"T2 = {T2:.1f} K | U2 = {U2:.0f} m/s | rho2 = {rho2:.3f} kg/m3 | p2 = {p2/1e3:.1f} kPa | pt = {pt/1e3:.1f} kPa | T_aw = {Taw:.1f} K")

# y+ = 1 y espesor de capa limite turbulenta a la salida
Cf = 0.0592*Re2**-0.2; ut = U2*np.sqrt(Cf/2); y1 = mu2/(rho2*ut)
delta = 0.37*C*Re2**-0.2
n_capas = int(np.ceil(np.log(1 + 1.2*delta*(1.2 - 1)/y1)/np.log(1.2)))
print(f"Cf = {Cf:.5f} | u_tau = {ut:.1f} m/s | primera celda (y+=1) = {y1*1e6:.2f} um | delta = {delta*1e3:.2f} mm | capas con razon 1.2 = {n_capas}")

# ---------------- coeficiente h del lado gas ----------------
L_humedo = 0.19                       # longitud mojada aprox. por lado (m)
Rex = rho2*U2*L_humedo/mu2
Nu_pl = 0.0296*Rex**0.8*Pr**(1/3); h_plano = Nu_pl*k_aire(T2)/L_humedo       # placa plana turbulenta (promedio)
# punto de estancamiento (cilindro, Frossling) con condiciones de entrada: Ma1 ~ 0.17
M1 = 0.17; T1 = Tt/(1 + 0.2*M1**2); U1 = M1*np.sqrt(gam*Rg*T1); p1 = pt/(1 + 0.2*M1**2)**3.5; rho1 = p1/(Rg*T1)
D_le = 2*0.01168; ReD = rho1*U1*D_le/mu_suth(T1)
h_est = 1.14*ReD**0.5*Pr**0.4*k_aire(T1)/D_le
print(f"h_g medio (placa plana turbulenta) = {h_plano:.0f} W/m2K | h_g en el borde de ataque (cilindro) = {h_est:.0f} W/m2K")
FACTOR_TU = 1.25                      # aumento por turbulencia de corriente libre (Tu ~ 6.5 %), estimacion de ingenieria
h_g = FACTOR_TU*h_plano
print(f"h_g de diseno = {FACTOR_TU} x h_plano = {h_g:.0f} W/m2K (rango esperado 700-1500 W/m2K segun la zona)")

# ---------------- canales de refrigerante: Dittus-Boelter ----------------
Tc_in, k_s = 300.0, 19.0            # K, W/mK (acero inoxidable tipo 310, como el C3X)
canales = list(csv.reader(open("../geometria/canales.csv")))[1:]
D = [float(c[3])*1e-3 for c in canales]
ReD_c = 30000.0
filas = []
for i, d in enumerate(D, 1):
    mu_c = mu_suth(Tc_in); k_c = k_aire(Tc_in)
    Nu = 0.023*ReD_c**0.8*Pr**0.4; h_c = Nu*k_c/d
    mdot = ReD_c*mu_c*np.pi*d/4
    filas.append((i, d*1e3, h_c, mdot*1e3))
    print(f"Canal {i:2d}: D = {d*1e3:4.1f} mm  Nu = {Nu:5.1f}  h_c = {h_c:5.0f} W/m2K  mdot = {mdot*1e3:5.2f} g/s")
h_c_prom = np.mean([f[2] for f in filas])

# ---------------- balance de energia del alabe (metal isotermo) ----------------
perim_gas = 2*0.19*1.0              # m de perimetro mojado (ambos lados) por m de envergadura
perim_cool = sum(np.pi*d for d in D)
T_m = (h_g*perim_gas*Taw + h_c_prom*perim_cool*Tc_in)/(h_g*perim_gas + h_c_prom*perim_cool)
q_lin = h_g*perim_gas*(Taw - T_m)    # W por metro de envergadura
q_slice = q_lin*0.0127
mdot_tot = sum(f[3] for f in filas)*1e-3
dTc = q_lin*0.0762/(mdot_tot*cp)
print(f"Metal isotermo (cota baja de gradientes): T_m = {T_m:.0f} K | Tw/Tg = {T_m/Tt:.2f}")
print(f"Flujo de calor por corte de 12.7 mm = {q_slice:.1f} W | calentamiento del refrigerante en 76.2 mm = {dTc:.1f} K")
Tc_ref = Tc_in + 0.5*dTc            # temperatura media del refrigerante a media envergadura
print(f"T del refrigerante a media envergadura (para la condicion de convección en el corte) = {Tc_ref:.0f} K")
print(f"Biot gas (t = 10 mm): {h_g*0.010/k_s:.2f}  -> el metal NO es isotermo; esperar gradientes de decenas de K")

# ---------------- presupuesto de malla (Fluent Student: 1 M de celdas) ----------------
A_gas, A_sol, perim = 23623.0, 3685.0, 380.0      # mm2 en XY (del dominio), perimetro mojado total (mm)
def celdas(h_bulk, h_wall, nz, capas=n_capas, ratio=1.0):
    h = h_bulk*ratio; hw = h_wall*ratio
    bulk = A_gas/h**2 + A_sol/(h*0.6)**2
    pared = perim/hw*capas
    return int((bulk + pared)*nz)
niveles = [("gruesa", 2.4, 1.0, 6), ("media", 1.6, 0.7, 9), ("fina", 1.1, 0.5, 13)]
print("\nPresupuesto de celdas (estimacion con celdas hexaedricas por barrido):")
filas_m = []
for nom, hb, hw, nz in niveles:
    n = celdas(hb, hw, nz)
    filas_m.append((nom, hb, hw, nz, n)); print(f"  {nom:7s}: h_bulk = {hb} mm, paso en pared = {hw} mm, {nz} capas en Z -> {n:,d} celdas {'(OK < 1 M)' if n < 1e6 else '(EXCEDE)'}")
with open("../resultados/presupuesto_malla.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["malla", "h_bulk_mm", "paso_pared_mm", "capas_z", "celdas_estimadas"]); w.writerows(filas_m)
with open("../resultados/referencia_analitica.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["variable", "valor", "unidad"])
    for k, v, u in [("Tt_in", Tt, "K"), ("Ma_salida", M2, "-"), ("Re_salida", Re2, "-"), ("T_estatica_salida", T2, "K"), ("p_estatica_salida", p2, "Pa"),
                    ("pt_entrada", pt, "Pa"), ("T_aw", Taw, "K"), ("y1_y+=1", y1, "m"), ("capas_pared", n_capas, "-"),
                    ("h_g_medio", h_g, "W/m2K"), ("h_g_estancamiento", h_est, "W/m2K"), ("h_c_prom_DB", h_c_prom, "W/m2K"),
                    ("T_metal_isotermo", T_m, "K"), ("Tc_media_envergadura", Tc_ref, "K"), ("h_c_canal_6.3mm", filas[0][2], "W/m2K"), ("h_c_canal_3.1mm", filas[-1][2], "W/m2K"), ("Re_D_canal", ReD_c, "-")]:
        w.writerow([k, f"{v:.6g}", u])
