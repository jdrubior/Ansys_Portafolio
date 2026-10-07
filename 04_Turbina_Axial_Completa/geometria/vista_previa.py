import numpy as np, cadquery as cq, pyvista as pv
pv.OFF_SCREEN = True

def malla(path, tol):
    s = cq.importers.importStep(path).val()
    v, t = s.tessellate(tol, 0.3)
    V = np.array([(p.x, p.y, p.z) for p in v]); T = np.array(t)
    faces = np.hstack([np.full((len(T), 1), 3), T]).ravel()
    return pv.PolyData(V, faces)

rotor = malla("rotor.step", 0.15); estator = malla("estator.step", 0.25)
def escena(corte, vista, salida, tam=(1300, 900)):
    p = pv.Plotter(off_screen=True, window_size=tam)
    p.set_background("white")
    est = estator.clip(normal=(0, -1, 0), origin=(0, 0, 0)) if corte else estator
    p.add_mesh(est, color="#aab4c0", smooth_shading=True, specular=0.4, opacity=1.0)
    p.add_mesh(rotor, color="#1f6fc5", smooth_shading=True, specular=0.6, specular_power=20)
    p.camera_position = vista
    p.enable_anti_aliasing()
    p.screenshot(salida)
escena(True, [(300, -380, 300), (0, 0, 70), (0, 0, 1)], "previa_turbina_3d_iso.png")
escena(False, [(0, 0, 520), (0, 0, 70), (0, 1, 0)], "previa_turbina_3d_frontal.png", (900, 900))
