import vtk

class Element:
    def __init__(self, id, number_of_nodes, material_set_id, node_ids=None):
        self._id = id
        self._material_set_id = material_set_id
        self._number_of_nodes = number_of_nodes
        if node_ids is None:
            self._node_ids = []
        else:
            self._node_ids = node_ids
    
    @property
    def id(self):
        return self._id
        
    @property
    def material_set_id(self):
        return self._material_set_id

    @property
    def number_of_nodes(self):
        return self._number_of_nodes

    @property
    def node_ids(self):
        return self._node_ids

    @node_ids.setter
    def node_ids(self, node_ids):
        if self._node_ids is None:
            self._node_ids = node_ids
        
class Node_2d:
    def __init__(self, nro, x_coord, y_coord):
        self._nro = nro
        self._x = x_coord
        self._y = y_coord

    @property
    def nro(self):
        return self._nro
    
    @property
    def coords(self):
        return [self._x, self._y]

class MeshData:
    def __init__(self, nodes_filename, elements_filename):
        self._list_of_nodes = []
        self._list_of_elements = []
        self.load_nodes(nodes_filename)
        self.load_elements(elements_filename)

    # Read elements from a file
    # nro, nnodes, mat_nro, nodes
    def load_elements(self, elements_filename): 
        with open(elements_filename, 'r') as elmt_data: 
            lines = elmt_data.readlines()
            for line in lines:
                elmt = list(map(int, line.split()))
                id, number_of_nodes, material_set_id, *node_ids = elmt
                node_ids = node_ids[:number_of_nodes]
                self._list_of_elements.append( Element(id, number_of_nodes, material_set_id, node_ids) )

    # Read nodes from a file
    # nro, x-coord, y-coord
    def load_nodes(self, nodes_filename):
        with open(nodes_filename, 'r') as node_data:
            lines = node_data.readlines()
            for line in lines:
                first, *rest = line.split()
                id = int(first)
                coords = list(map(float, rest))
                self._list_of_nodes.append( Node_2d(id, *coords) )

    @property
    def nodes(self):
        return self._list_of_nodes

    @property
    def elements(self):
        return self._list_of_elements

    @property
    def number_of_nodes(self):
        return len(self._list_of_nodes)

    @property
    def number_of_elements(self):
        return len(self._list_of_elements)

    def create_vtk_object(self): 
        # Create VTK object for the mesh data
        vtk_dataset = vtk.vtkUnstructuredGrid() 

        # Create the points by defining their coordinates
        points = vtk.vtkPoints()
        for id in range(self.number_of_nodes):
            points.InsertPoint(id, [*self.nodes[id].coords, 0])
        vtk_dataset.SetPoints(points)

        # Create the cells by specifying connectivity
        vtk_dataset.Allocate(self.number_of_elements)
        for id in range(self.number_of_elements):
            point_ids = list(map(lambda a: a - 1, self.elements[id].node_ids))
            num_point_ids = self.elements[id].number_of_nodes
            vtk_dataset.InsertNextCell(7,num_point_ids,point_ids)

        # Create data arrays (on cells)
        array = vtk.vtkIntArray()
        array.SetName("material_set_id")
        array.SetNumberOfValues(self.number_of_elements)
        for id in range(self.number_of_elements):
            array.SetValue(id, self.elements[id].material_set_id)
        vtk_dataset.GetCellData().AddArray(array)
        return vtk_dataset

