import numpy as np, cadquery as cq, pyvista as pv
pv.OFF_SCREEN = True
def malla(path, tol):
    s = cq.importers.importStep(path).val()
    v, t = s.tessellate(tol, 0.3)
    V = np.array([(p.x, p.y, p.z) for p in v]); T = np.array(t)
    return pv.PolyData(V, np.hstack([np.full((len(T), 1), 3), T]).ravel())
gas = malla("dominio_gas_sector.step", 0.25); rot = malla("sector_rotor_solido.step", 0.15)
p = pv.Plotter(off_screen=True, window_size=(1300, 900)); p.set_background("white")
p.add_mesh(gas, color="#f0a35e", opacity=0.28, smooth_shading=True)
p.add_mesh(gas.extract_feature_edges(40), color="#b45f06", line_width=1)
p.add_mesh(rot, color="#1f6fc5", smooth_shading=True, specular=0.6, specular_power=20)
p.camera_position = [(330, 60, 250), (70, 35, 75), (0, 0, 1)]
p.enable_anti_aliasing(); p.screenshot("previa_sector_45.png")
