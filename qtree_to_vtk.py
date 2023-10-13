import sys
import vtk
sys.path.append("./src/")
from meshData import MeshData

input_filename = ""
vtk_filename = "./Output/parv/QtreeMesh.vtu"

if len(sys.argv) == 2:
    input_filename = sys.argv[1]
    vtk_filename = "./Output/parv/" + input_filename.split('/')[-1].replace(".txt",".vtu")
print(vtk_filename)

nodes_filename = './Output/mesh/scor.txt'
elements_filename = './Output/mesh/selm.txt'

# Load qtree mesh data
qtreeData = MeshData(nodes_filename, elements_filename)
nodes = qtreeData.nodes
elements = qtreeData.elements
number_of_nodes = qtreeData.number_of_nodes
number_of_elements = qtreeData.number_of_elements


writer = vtk.vtkXMLUnstructuredGridWriter()
writer.SetFileName(vtk_filename)
writer.SetInputData(qtree_vtk_dataset)
writer.Write()
