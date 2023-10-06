import sys
import vtk
sys.path.append("./src/")
from elmtData import Node_2d, Element

list_of_nodes = []
list_of_elements = []

# Read elements from a file
# nro, nnodes, mat_nro, nodes
with open('./Output/mesh/selm.txt') as elmt_data:
    lines = elmt_data.readlines()
    for line in lines:
        elmt = list(map(int, line.split()))
        id, number_of_nodes, material_set_id, *node_ids = elmt
        node_ids = node_ids[:number_of_nodes]
        list_of_elements.append( Element(id, number_of_nodes, material_set_id, node_ids) )

# Read nodes from a file
# nro, x-coord, y-coord
with open('./Output/mesh/scor.txt') as node_data:
    lines = node_data.readlines()
    for line in lines:
        first, *rest = line.split()
        id = int(first)
        coords = list(map(float, rest))
        list_of_nodes.append( Node_2d(id, *coords) )

number_of_elements = len(list_of_elements)
number_of_nodes = len(list_of_nodes)

# Create VTK object for the data model
qtree_vtk_dataset = vtk.vtkUnstructuredGrid() 

# Create the points by defining their coordinates
points = vtk.vtkPoints()
for id in range(number_of_nodes):
    points.InsertPoint(id, [*list_of_nodes[id].coords, 0])
qtree_vtk_dataset.SetPoints(points)

# Create the cells by specifying connectivity
qtree_vtk_dataset.Allocate(number_of_elements)
for id in range(number_of_elements):
    point_ids = list(map(lambda a: a - 1, list_of_elements[id].node_ids))
    num_point_ids = list_of_elements[id].number_of_nodes
    qtree_vtk_dataset.InsertNextCell(7,num_point_ids,point_ids)

# Create data arrays (on cells)
array = vtk.vtkIntArray()
array.SetName("material_set_id")
array.SetNumberOfValues(number_of_elements)
for id in range(number_of_elements):
    array.SetValue(id, list_of_elements[id].material_set_id)
qtree_vtk_dataset.GetCellData().AddArray(array)

writer = vtk.vtkXMLUnstructuredGridWriter()
writer.SetFileName("output.vtu")
writer.SetInputData(qtree_vtk_dataset)
writer.Write()

print("Done!")
