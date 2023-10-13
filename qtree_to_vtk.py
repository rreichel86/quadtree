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

# Create VTK object for the data model
qtree_vtk_dataset = vtk.vtkUnstructuredGrid() 

# Create the points by defining their coordinates
points = vtk.vtkPoints()
for id in range(number_of_nodes):
    points.InsertPoint(id, [*nodes[id].coords, 0])
qtree_vtk_dataset.SetPoints(points)

# Create the cells by specifying connectivity
qtree_vtk_dataset.Allocate(number_of_elements)
for id in range(number_of_elements):
    point_ids = list(map(lambda a: a - 1, elements[id].node_ids))
    num_point_ids = elements[id].number_of_nodes
    qtree_vtk_dataset.InsertNextCell(7,num_point_ids,point_ids)

# Create data arrays (on cells)
array = vtk.vtkIntArray()
array.SetName("material_set_id")
array.SetNumberOfValues(number_of_elements)
for id in range(number_of_elements):
    array.SetValue(id, elements[id].material_set_id)
qtree_vtk_dataset.GetCellData().AddArray(array)

writer = vtk.vtkXMLUnstructuredGridWriter()
writer.SetFileName(vtk_filename)
writer.SetInputData(qtree_vtk_dataset)
writer.Write()
