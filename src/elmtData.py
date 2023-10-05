
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

