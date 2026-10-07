"""Referencia analitica/numerica 1D para el disco de freno (proyecto 3).

Modelo de validacion independiente de Ansys:
  * Energia de frenado (100 -> 0 km/h, desaceleracion constante)
  * Modelo concentrado (masa unica) -> temperatura media
  * Conduccion 1D en el espesor de la banda de friccion (diferencias finitas)
    con flujo de calor por friccion variable en el tiempo y conveccion.
Los parametros son REPRESENTATIVOS de un disco ventilado de auto compacto,
no medidas de un vehiculo especifico. Ajustar en PARAMS cuando se tengan.
"""
import csv
import numpy as np

PARAMS = dict(
    m_veh=1450.0,        # kg, vehiculo + carga
    v0=100/3.6,          # m/s
    a=0.8*9.81,          # m/s2, desaceleracion constante
    f_front=0.70,        # fraccion de frenado en el eje delantero
    r_in=0.100, r_out=0.140,   # m, radios de la pista de friccion
    t_disc=0.025,        # m, espesor total del disco ventilado
    t_wall=0.0085,       # m, espesor de cada placa de friccion
    m_disc=8.0,          # kg, masa del disco (hierro fundido)
    rho=7150.0, k=50.0, c=500.0,   # hierro fundido gris
    h=60.0,              # W/m2K convectivo medio (se reemplazara con Fluent)
    T_inf=25.0,
    gamma=0.90,          # fraccion de calor que entra al disco (resto, balata)
)

def potencia(t, p):
    """Potencia de frenado por disco en t (W)."""
    t_stop = p["v0"]/p["a"]
    v = np.maximum(p["v0"] - p["a"]*t, 0.0)
    P_total = p["m_veh"]*p["a"]*v
    return P_total*p["f_front"]/2.0*p["gamma"]*(t <= t_stop)

def energia_por_disco(p):
    return 0.5*p["m_veh"]*p["v0"]**2*p["f_front"]/2.0*p["gamma"]

def modelo_concentrado(p):
    Q = energia_por_disco(p)
    return Q/(p["m_disc"]*p["c"])

def simular_1d(p, n_stops=1, t_total=None, dt=5e-4, nx=60, t_cool=30.0):
    """Conduccion 1D en una placa de friccion: cara caliente (x=0) con flujo q,
    cara opuesta (x=L) adiabatica por simetria con los canales (aprox.).
    Conveccion en la cara caliente, h aplicada sobre toda la superficie."""
    L = p["t_wall"]
    dx = L/(nx-1)
    alpha = p["k"]/(p["rho"]*p["c"])
    assert alpha*dt/dx**2 < 0.5, "paso inestable"
    A = np.pi*(p["r_out"]**2-p["r_in"]**2)   # una cara
    t_stop = p["v0"]/p["a"]
    cycle = t_stop + t_cool
    if t_total is None:
        t_total = n_stops*cycle
    T = np.full(nx, p["T_inf"])
    ts, Ts, Tm = [], [], []
    t = 0.0
    while t < t_total:
        tc = t % cycle
        q = potencia(np.array(tc), p)/(2*A) if tc <= t_stop else 0.0
        Tn = T.copy()
        Tn[1:-1] = T[1:-1] + alpha*dt/dx**2*(T[2:]-2*T[1:-1]+T[:-2])
        # nodo 0: balance con flujo y conveccion (semi-celda)
        Tn[0] = T[0] + 2*dt/(p["rho"]*p["c"]*dx)*(q - p["h"]*(T[0]-p["T_inf"])
                          + p["k"]*(T[1]-T[0])/dx)
        Tn[-1] = T[-1] + 2*alpha*dt/dx**2*(T[-2]-T[-1])   # adiabatica
        T = Tn
        t += dt
        if int(round(t/dt)) % 100 == 0:
            ts.append(t); Ts.append(T[0]); Tm.append(T.mean())
    return np.array(ts), np.array(Ts), np.array(Tm)

if __name__ == "__main__":
    p = PARAMS
    t_stop = p["v0"]/p["a"]
    print(f"t_frenado = {t_stop:.2f} s")
    print(f"Energia por disco = {energia_por_disco(p)/1e3:.1f} kJ")
    print(f"P_max por disco   = {potencia(np.array(0.0), p)/1e3:.1f} kW")
    dT = modelo_concentrado(p)
    print(f"Modelo concentrado: dT media tras 1 frenada = {dT:.1f} K "
          f"-> T = {p['T_inf']+dT:.1f} C")
    ts, Ts, Tm = simular_1d(p, n_stops=1, t_total=t_stop+30.0)
    print(f"1D, 1 frenada: T_sup max = {Ts.max():.1f} C a t = {ts[Ts.argmax()]:.2f} s")
    ts5, Ts5, Tm5 = simular_1d(p, n_stops=5)
    print(f"1D, 5 frenadas (30 s entre): T_sup max = {Ts5.max():.1f} C")
    with open("../resultados/referencia_1d.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["t_s", "T_sup_C", "T_media_C"])
        for r in zip(ts5, Ts5, Tm5): w.writerow([f"{x:.4f}" for x in r])
