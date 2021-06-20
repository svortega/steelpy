# 
# Copyright (c) 2009-2021 fem2ufo
# 
# Python stdlib imports
from itertools import chain
from collections import Counter
from collections import defaultdict
import functools
import re
from typing import NamedTuple, Dict, Union, Tuple, List
from math import  isclose, dist
#
# package imports
#from steelpy.trave3D.preprocessor.assemble import get_bandwidth

#
def node_renumbering(nodes, elements):
    """ """
    connectivities = elements.get_connectivities
    nodes_conn = get_nodes_connected(nodes, connectivities)
    node_degrees = get_node_degrees(connectivities)
    new_nodes_number = get_node_renumber(nodes, node_degrees, nodes_conn)
    return new_nodes_number

#
def get_nodes_connected(nodes, connectivities):
    """
    """
    # split nodes in two lists
    connectivity = list(map(list, zip(*connectivities)))
    # get nodes connected
    nodes_conn = {}
    for key in nodes.keys():
        node_ends = [[n for n,x in enumerate(col) if x==key]
                     for col in connectivity]
        flat_list = list(chain(*node_ends))
        #
        nodes_conn[key] = [col[index] for col in connectivity
                           for index in flat_list if col[index] != key]
    return nodes_conn
#
def get_node_degrees(connectivities):
    """ """
    # Step 1 - get nodes degrees.
    degree = get_node_degree(connectivities)
    try:
        items = degree[0]
        raise IOError(f" free nodes {items}")
    except:
        # Step 2 - pick a starting node (i.e. the node with low degree).
        first_degree = next(iter(degree.values()))
        return first_degree
#
def get_node_renumber(nodes, first_degree, nodes_conn):
    """ """
    node_labels = [key for key in nodes.keys()]
    # Pick up one node from first level
    levels = []
    start_node = first_degree[0]
    levels.append([start_node])
    # remove node from list
    node_labels = del_list_inplace(node_labels, levels[-1])
    # define level 2
    levels.append(nodes_conn[start_node])
    node_labels = del_list_inplace(node_labels, levels[-1])
    #
    # define next levels
    rem = list(chain(*levels))
    steps = len(node_labels)
    _iter = 0
    while node_labels:
        _iter += 1
        rem = get_level(levels, nodes_conn, rem)
        del_list_inplace(node_labels, levels[-1])
        if _iter > steps:
            raise RuntimeError(" Node renumbering fail")
    #
    #_node_number = 0
    #_memb_number = 0
    #_memb_rem = []
    #for level in levels:
        # renumber node
        #for _node_name in sorted(level):
            #_node_number += 1
            #_index = nodes._labels.index(_node_name)
            #nodes._number[_index] = _node_number
            #nodes[_node_name].number = _node_number
            #nodes.update_number(node_name=_node_name, number=_node_number)
            # renumber member
            #for _member_name in sorted(nodes_conn[_node_name][0]):
            #    if _member_name in _memb_rem:
            #        continue
            #    _memb_rem.append(_member_name)
            #    _memb_number += 1
            #    #elements.update_item(element_number=_member_name, item="number", 
            #    #                     value=_memb_number)
            #    elements[_member_name].number = _memb_number
    #print('end renumbering nodes')
    new_nodes_number = list(chain(*levels))
    return new_nodes_number
#
#
def get_node_degree(connectivities):
    """
    Scan all the nodes and order them according to their dregree.
    The degree of a node is the number of nodes connected to it
    """
    flat_nodes = list(chain.from_iterable(connectivities))
    shared_nodes = Counter(flat_nodes)
    degree = defaultdict(list)
    for key, value in sorted(shared_nodes.items()):
        degree.setdefault(value, []).append(key)
    return degree
#
#
def get_level(levels, nodes_conn, rem):
    """
    """
    cases = {}
    rem2 = []
    for item in levels[-1]:
        #nodes_rem = set(nodes_conn[item]) - set(rem)
        cases[item] = set(nodes_conn[item]) - set(rem)
        rem2.extend(cases[item])
    rem2 = list(set(rem2))
    levels.append(rem2)
    return rem + rem2
#
def del_list_inplace(lst, id_to_del):
    """ """
    for item in set(id_to_del):
        while item in lst:
            lst.remove(item)
    return lst
#
#
#
@functools.lru_cache(maxsize=2048) 
def find_node_data(word_in: str) -> str:
    """
    Identify beam data from user
    """
    _key: Dict = {"number": r"\b(number|mesh)\b",
                  "name": r"\b(name|label)\b",
                  "elements": r"\b(element|member|item(s)?)\b",
                  "group": r"\b(group|set)(s)?\b",
                  "boundary": r"\b(boundar(y|ies))\b",
                  "z": r"\b(z)\b",
                  "y": r"\b(y)\b",
                  "x": r"\b(x)\b",
                  "coordinates": r"coordinates"}
    
    _match = common.find_keyword(word_in, _key)
    return _match

@functools.lru_cache(maxsize=2048) 
def find_element_data(word_in: str) -> str:
    """
    Identify beam data from user
    """
    _key: Dict = {"number": r"\b(number|mesh)\b",
                  "name": r"\b(name|label)\b",
                  "connectivity": r"\b(connectivity|node(s)?|joint(s)?)\b",
                  "material": r"\b(material)\b",
                  "section": r"\b(geometry|section|shape)\b",
                  "hinges": r"\b(hinge(s)?)\b",
                  "group": r"\b(group|set)(s)?\b",
                  # "boundary": r"\b(boundar(y|ies))\b",
                  "type": r"type"}
    
    _match = common.find_keyword(word_in, _key)
    return _match

#
def get_args(args, items, item_class, item_type):
    """
    """
    if len(args) == 1:
        for _arg in args[0]:
            if type(_arg) in [str, int, float]:
                try:
                    items[_arg]
                    raise Exception('{:} {:} already exist'
                                    .format(_arg, item_type))
                except KeyError:
                    items[_arg] = item_class(*args[0])
                break
            else:
                try:
                    items[_arg[0]]
                    raise Exception('{:} {:} already exist'
                                    .format(_arg[0], item_type))
                except KeyError:
                    items[_arg[0]] = item_class(*_arg)
    else:
        try:
            items[args[0]]
            raise Exception('{:} {:} already exist'
                            .format(args[0], item_type))
        except KeyError:
            items[args[0]] = item_class(*args)
#
#
#
#
class CoordCartesian(NamedTuple):
    """ Cartesian coordinate system"""
    x: float
    y: float
    z: float
    name: int
    number: int
    index:int
    system:str ="cartesian"
    #boundaries: Iterable
    #sets: List[Tuple]
    #
    #def get_coordinates(self, units:str='metre'):
    #    """
    #    """
    #    return "{:14.0f} {: 14.5f} {: 14.5f} {: 14.5f}".format(self.number,
    #                                                           self.x, self.y, self.z)
    #
    #@property
    #def boundary(self) -> Tuple:
    #    """
    #    """
    #    return self.boundaries[self.number]
    #
    #@boundary.setter
    #def boundary(self, value: List) -> None:
    #    """
    #    """
    #    self.boundaries[self.number] = value

    def __str__(self) -> str:
        return "{:12d} {: 12.5f} {: 12.5f} {: 12.5f}\n"\
            .format(self.name, self.x, self.y, self.z)

    def __eq__(self, other) -> bool:
        """
        """
        if (isclose(self.x, other.x, abs_tol=1e-03)
            and isclose(self.y, other.y, abs_tol=1e-03)
            and isclose(self.z, other.z, abs_tol=1e-03)):
            return True
        return False
    #
    def distance(self, other):
        """ distance between two nodes"""
        #print("here")
        node1 = [self.x, self.y, self.z]
        return dist(node1[:3], other[:3])      

class CoordCylindrical(NamedTuple):
    """
    """
    r: float
    theta: float
    z: float
    name: int
    number: int
    index: int
    system:str ="cylindrical"

class CoordSpherical(NamedTuple):
    """
    """
    r: float
    theta: float
    phi: float
    name: int
    number: int
    index: int
    system:str ="spherical"
#
#
def get_coordinate_system(system):
    """
    """
    if 'cylindrical' in system.lower():
        return CoordCylindrical
    elif 'spherical' in system.lower():
        return CoordSpherical
    else:
        return CoordCartesian
#
#
def check_point_list(data, steps:int=6)->List[float]:
    """ """
    new_data = []
    for x in range(steps):
        try:
            try:
                new_data.append(data[x].value)
            except AttributeError:
                new_data.append(data[x])
        except IndexError:
            new_data.append(0.0)
    return new_data
#
def check_point_dic(data)->List[float]:
    """ """
    new_data = [0,0,0]
    for key, item in data.items():
        if re.match(r"\b(x)\b", str(key), re.IGNORECASE):
            new_data[0] = item.value
        elif re.match(r"\b(y)\b", str(key), re.IGNORECASE):
            new_data[1] = item.value
        elif re.match(r"\b(z)\b", str(key), re.IGNORECASE):
            new_data[2] = item.value
    return new_data
#
#
    


