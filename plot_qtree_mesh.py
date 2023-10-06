from matplotlib.colors import ListedColormap 
import matplotlib.pyplot as plt

import pyvista
try: 
    mesh = pyvista.read('./Output/parv/QtreeMesh.vtu')
    pyvista.set_jupyter_backend('client')
    boring_cmap = plt.matplotlib.colormaps["rainbow"]
    mesh.plot(scalars='material_set_id', show_edges=True, cpos="xy", color=True, cmap=boring_cmap)
except FileNotFoundError:
    print("File is not found.")

