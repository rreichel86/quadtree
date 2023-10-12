
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
    def __inti__(self):
        self._list_of_nodes = []
        self._list_of_elements = []
        load_elements()
        load_nodes()

    def load_elements(self):
        with open('./Output/mesh/selm.txt') as elmt_data:
    # Read elements from a file
    # nro, nnodes, mat_nro, nodes
            lines = elmt_data.readlines()
            for line in lines:
                elmt = list(map(int, line.split()))
                id, number_of_nodes, material_set_id, *node_ids = elmt
                node_ids = node_ids[:number_of_nodes]
                self._list_of_elements.append( Element(id, number_of_nodes, material_set_id, node_ids) )

    def load_nodes(self):
        with open('./Output/mesh/scor.txt') as node_data:
    # Read nodes from a file
    # nro, x-coord, y-coord
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
