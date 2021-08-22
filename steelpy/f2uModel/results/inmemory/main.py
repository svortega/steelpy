#
# Copyright (c) 2009-2021 fem2ufo
#

# Python stdlib imports
from array import array
from collections import defaultdict, Mapping
from typing import NamedTuple, Dict, List, Tuple, Union
#
#
# package imports
from steelpy.f2uModel.results.operations.output import (Results, print_deflections,
                                                        print_member_forces, print_node_reactions)
#
#
# --------------------------
#
#
class ResultInmemory:
    
    __slots__ = ['_displacement', '_beam_force', 
                 '_basic_loads', '_reaction']

    def __init__(self, basic_loads):
        """
        """
        self._basic_loads = basic_loads
        self._displacement = NodeDisplacement()
        self._beam_force = BeamForce()
        self._reaction = BoundaryNodes()
    #
    #
    @property
    def node_displacement(self):
        """
        [node, load_tile, system, x, y, z, rx, ry, rz]
        """
        return self._displacement._get_node_displacements(self._basic_loads)
    #
    @node_displacement.setter
    def node_displacement(self, values:List):
        """
        values = [node, load_tile, system, x, y, z, rx, ry, rz]
        """
        for value in values:
            self._displacement[value[0]] = value[1:]
        #print("-->")
    #
    def print_node_displacement(self):
        """
        """
        print_deflections(self.node_displacement)
    #
    @property
    def element_force(self):
        """
        [element, load_title, node_name, position, system, fx, fy, fz, mx, my, mz]
        """
        return self._beam_force._get_element_forces(self._basic_loads)

    @element_force.setter
    def element_force(self, values:List):
        """
        values = [element, load_title, node_name, position, system, fx, fy, fz, mx, my, mz]
        """
        for value in values:
            self._beam_force[value[0]] = value[1:]        
        #print("-->")
    #
    def print_element_forces(self):
        """ """
        print_member_forces(self.element_force)
    #
    #
    #
    @property
    def node_reaction(self):
        """
        [node, load_tile, system, x, y, z, rx, ry, rz]
        """
        return self._reaction._get_node_displacements(self._basic_loads)
    #
    @node_reaction.setter
    def node_reaction(self, values:List):
        """
        values = [node, load_tile, system, x, y, z, rx, ry, rz]
        """
        for value in values:
            self._reaction[value[0]] = value[1:]
        #print("-->")
    #
    def print_node_reactions(self):
        """
        """
        print_node_reactions(self.node_reaction)
    #
#
#
def iter_items(rows):
    """ """
    member_load = {}
    for row in rows:
        try:
            member_load[row[0]].extend(list(row[2:]))
        except KeyError:
            member_load[row[0]] = list(row[2:])
    return member_load
#
#
#
class NodeBasic:
    __slots__ = ['_system', '_labels', '_system',
                 '_load_title', '_load_number',
                 '_x', '_y', '_z', '_rx', '_ry', '_rz',]

    def __init__(self):
        """
        """
        self._x: array  = array('f',[])
        self._y: array  = array('f',[])
        self._z: array  = array('f',[])
        self._rx: array = array('f',[])
        self._ry: array = array('f',[])
        self._rz: array = array('f',[])
        #
        self._labels:  array = array("I", [])
        self._load_title: List[str] = []
        self._load_number: array = array("I", [])
        self._system: List[str] = []
    #
    def __setitem__(self, node_number:int, 
                    result:List) -> None:
        """
        results = [load_title, system, x, y, z, rx, ry, rz]
        """
        self._labels.append(node_number)
        self._load_title.append(result[0])
        #
        self._system.append(result[1])
        self._x.append(result[2])
        self._y.append(result[3])
        self._z.append(result[4])
        self._rx.append(result[5])
        self._ry.append(result[6])
        self._rz.append(result[7])
    #
    #
    #def __getitem__(self, node_number: int)-> List[Tuple]:
    #    """
    #    """
    #    print("-->")
    #
    def _get_node_displacements(self, basic_loads):
        """ """
        bloads = {item.title:key for key, item in basic_loads.basic.items()}
        basic = self._get_displacement(bloads, "basic")
        #
        step = len(basic)+1
        bloads = {item.title:key for key, item in basic_loads.combination.items()}
        comb = self._get_displacement(bloads, "combination", step=step)
        #
        displacement = {"basic": basic, "combination": comb}
        return displacement
    #
    def _get_displacement(self, bloads:Dict, ltype:str, step:int=1):
        """ """
        node_res = defaultdict(list)
        for index, lname in enumerate(self._load_title):
            if lname not in bloads:
                continue
            name = self._labels[index]
            node_res[lname].append([name, self._x[index], self._y[index], self._z[index],
                                    self._rx[index], self._ry[index], self._rz[index]])
        #
        basic = {}
        for x, (key, item) in enumerate(node_res.items()):
            lname = bloads[key]
            basic[key] = Results(name=lname, number=x+step, title=key,
                                 load_type=ltype, items=item)
        return basic
#
#
class NodeDisplacement(NodeBasic):
    #__slots__ = ['_system', '_labels', '_system',
    #             '_load_title', '_load_number',
    #             '_x', '_y', '_z', '_rx', '_ry', '_rz',]

    def __init__(self):
        """
        """
        super().__init__()
    #

#
#
class BeamForce:
    __slots__ = ['_system', '_labels', '_system',
                 '_load_title', '_load_number',
                 '_fx', '_fy', '_fz', 
                 '_mx', '_my', '_mz',
                 '_node', '_position']

    def __init__(self):
        """
        """
        self._fx: array  = array('f',[])
        self._fy: array  = array('f',[])
        self._fz: array  = array('f',[])
        self._mx: array = array('f',[])
        self._my: array = array('f',[])
        self._mz: array = array('f',[])
        #
        self._node : array = array("I", [])
        self._position: array = array('f',[]) 
        #
        self._labels:  array = array("I", [])
        self._load_title: List[str] = []
        self._load_number: array = array("I", [])
        self._system: List[str] = []
    #
    def __setitem__(self, element_number:int, 
                    result:List) -> None:
        """
        results = [load_title, system, node, fx, fy, fz, mx, my, mz]
        """
        self._labels.append(element_number)
        self._load_title.append(result[0])
        # node and node position
        self._node.append(result[1])
        self._position.append(result[2])
        self._system.append(result[3])
        # force
        self._fx.append(result[4])
        self._fy.append(result[5])
        self._fz.append(result[6])
        self._mx.append(result[7])
        self._my.append(result[8])
        self._mz.append(result[9])
    #
    #def __getitem__(self, element_number: int)-> List[Tuple]:
    #    """
    #    """
    #    print("-->")
    #
    def _get_element_forces(self, basic_loads):
        """ """
        bloads = {item.title:key for key, item in basic_loads.basic.items()}
        memf_basic = self._get_force(bloads, "basic")
        #
        step = len(memf_basic)+1
        bloads = {item.title:key for key, item in basic_loads.combination.items()}
        memf_comb = self._get_force(bloads, "combination", step=step)
        #
        forces = {"basic":memf_basic, "combination":memf_comb}
        return forces
    #
    def _get_force(self, bloads:Dict, ltype:str, step:int=1):
        """ """
        member_res = defaultdict(list)
        for index, lname in enumerate(self._load_title):
            if lname not in bloads :
                continue
            #if self._system[index] != system:
            #    continue
            name = self._labels[index]
            member_res[lname].append([name, self._node[index], 
                                      self._position[index], self._system[index],
                                      self._fx[index], self._fy[index], self._fz[index],
                                      self._mx[index], self._my[index], self._mz[index]])
        #
        basic = {}
        for x, (key, item) in enumerate(member_res.items()):
            lname = bloads[key]
            basic[key] = Results(name=lname, number=x+step, title=key,
                                 load_type=ltype, items=item)
        return basic
    #
    #def _get_forces_items(self):
    #    """ """
    #    pass
#
#
class BoundaryNodes(NodeBasic):
    
    def __init__(self):
        """
        """
        super().__init__()
    #
    #
    #def _get_node_reaction(self):
    #    """" """
    #    pass


#
class Boundary:
    
    def __init__(self) -> None:
        """
        """
        self._nodes = BoundaryNodes()
    #
    @property
    def node(self):
        """"""
        return self._nodes
    
    @node.setter
    def node(self, values):
        """"""
        for value in values:
            self._nodes[value[0]] = value[1:]
#
