# 
# Copyright (c) 2009-2018 fem2ufo
# 


# Python stdlib imports
from typing import Tuple, Dict, List, ClassVar, Iterable, MutableMapping
from collections.abc import Mapping


class BasicShell:
    """
    """
    __slots__ = ['_materials', '_points', '_type', '_nodes']
                # '_elements',  'overlap', 
                # '_properties', '_offsets',  '_section',
                # '_mesh' # '_points',
    
    def __init__(self, shell_type:str, points: ClassVar,
                 materials: ClassVar) -> None:
        """
        """
        self._type: str = shell_type
        self._points: ClassVar = points
        self._materials: ClassVar = materials
        #
        # self._properties = BeamProperties()
    #
    @property
    def corner_points(self) -> Tuple:
        """
        """
        if len(self._nodes) == 3:
            return self._points[self._nodes[0]],\
                   self._points[self._nodes[1]],\
                   self._points[self._nodes[2]]
        else:
            return self._points[self._nodes[0]],\
                   self._points[self._nodes[1]],\
                   self._points[self._nodes[2]],\
                   self._points[self._nodes[4]]  

    @corner_points.setter
    def corner_points(self, points: List):
        """
        """
        _nodes_no = len(points)
        if _nodes_no < 3:
            raise IOError('   *** Joint number {:} < 3'
                          .format(_nodes_no))
        
        try:
            self._nodes
            raise Warning("points already exists")
        except AttributeError:
            self._nodes = []
            if isinstance(points, (list, tuple)):
                for point_name in points:
                    try:
                        self._points[point_name]
                        self._nodes.append(point_name)
                    except KeyError:
                        raise IOError('point {:} does not exist'.format(point_name))
            elif isinstance(points, dict):
                pass
            else:
                raise IOError('   *** Joint input format {:} not recognized'.format(points))
        #      
    #
    @property
    def materials(self) -> ClassVar:
        """
        """
        return self._materials    
    #
    @property
    def type(self):
        """
        """
        return self._type

#
#
class Shells(Mapping):
    """
    """
    __slots__ = ['_shells', '_points', '_shell_type',
                 '_materials']
    
    def __init__(self, shell_type:str, points:ClassVar, 
                 materials:MutableMapping) -> None:
        """
        """
        self._shell_type:str = shell_type
        self._points: ClassVar = points
        self._shells: Dict = {}
        self._materials:MutableMapping = materials
    #
    def __getitem__(self, shell_name: str)-> ClassVar:
        """
        """
        return self._shells[shell_name]
    #
    def __setitem__(self, shell_name: str, corner_nodes: List[float]) -> None:
        """
        """
        self._shells[shell_name] = BasicShell(shell_type=self._shell_type,
                                              points = self._points,
                                              materials= self._materials)
        if corner_nodes:
            self._shells[shell_name].corner_points = corner_nodes
    #
    #
    def __len__(self) -> int:
        return len(self._shells)

    
    def __iter__(self)-> Iterable:
        """
        """
        return iter(self._shells)
        #

    def __contains__(self, value) -> bool:
        return value in self._shells
    #    
    
#
#
class ShellsXX:
    """
    FE Shell Element class  
    
    Element 
        |_ name
        |_ number
        |_ type
        |_ material [mt[0],..]
        |_ section [st[0],...]
        |_ node [nd[0], nd[1], nd[2],..., nd[n]]
        |_ offset [ecc[0], ecc[1], ecc[2],..., ecc[n]]
        |_ releases [rel[0], rel[1], rel[2],..., rel[n]]
        |_ guidepoint [v[0], v[1], v[2],..., v[n]]
        |_ mass
        |_ stress
        |_ code_check
        |_ concept
        |_ hydrodynamics
        |_ aerodynamics
    
    **Parameters**:  
      :number:  integer internal number 
      :name:  string node external name
    """
    __slots__ = ('number', 'name', 'type', 'material', 'node', 'load',
                 'section', 'offsets', 'releases', 'guidepoint', 'mass',
                 'stress', 'concept', 'code_check', 'group', 'overlap',
                 'properties')
    #
    #_nodes = Node()

    #
    def __init__(self, Name, Number=None, Type='N/A'):
        self.number = Number
        self.name = Name
        self.type = Type
        self.material = []
        self.node = []
        self.section = []
        self.offsets = []
        self.releases = []
        self.guidepoint = []
        self.stress = []
        self.mass = []
        self.load = {}
        self.overlap = False
        self.properties = None

    #
    def get_sum_distributed_load(self, load_id, load_factor=1):
        """
        Flats udl and sum all the udl loads within the list
        """
        if not self.load[load_id].pressure:
            return False

        _udl = []
        #
        for _pressure in self.load[load_id].pressure:
            _udl.append(_pressure)
        #
        #
        try:
            _total_0 = [sum(i) for i in zip(*_udl)]
            _total = [_load * load_factor for _load in _total_0]
            return _total
        except IndexError:
            return []

    #
    def line_length(self, line):
        """
        When line load is specified, the relation between local node numbers and loaded line will be:
        4 nodes:
            LINE =1 means line load between node 1 and 2 
            LINE =2 means line load between node 2 and 3 
            LINE =3 means line load between node 3 and 4 
            LINE =4 means line load between node 4 and 1
        3 nodes:
            Line 1 means line load between nodes 2 and 3. 
            Line 2 means line load between nodes 1 and 3. 
            Line 3 means line load between nodes 1 and 2.
        """
        _number_nodes = len(self.node)
        end_1, end_2 = get_node_end(_number_nodes, line)
        return get_node2node_length(self._nodes, self.node, end_1, end_2)

        # if _number_nodes == 3:
        # print('triangle')
        # if line == 1:
        #    return get_node2node_length(self.node, end_1=1, end_2=2)
        # elif line == 2:
        #    return get_node2node_length(self.node, end_1=0, end_2=2)
        # else:
        #    return get_node2node_length(self.node, end_1=0, end_2=1)
        # else:
        # print('quad')
        # if line == 1:
        #    return get_node2node_length(self.node, end_1=0, end_2=1)
        # elif line == 2:
        #    return get_node2node_length(self.node, end_1=1, end_2=2)
        # elif line == 3:
        #    return get_node2node_length(self.node, end_1=2, end_2=3)
        # else:
        #    return get_node2node_length(self.node, end_1=3, end_2=0)
        #
        # return


#