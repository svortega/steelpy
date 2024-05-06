# 
# Copyright (c) 2009 steelpy
# 
#
# Python stdlib imports
from __future__ import annotations
from array import array
from dataclasses import dataclass
from itertools import chain, count
from collections import Counter
from collections import defaultdict
from collections.abc import Mapping
import functools
import re
from typing import NamedTuple
from math import isclose, dist
#
# package imports
from steelpy.utils.math.operations import zeros, to_matrix
from steelpy.utils.dataframe.main import DBframework


#
#
class NodeBasic(Mapping):
    """
    This is a fem2ufo model node class


    Parameters
    ----------
    boundaries: object
        f2u boundary object


    Attributes
    ----------
    labels : array
       node internal number
    x : array
       coordinate x
    y : array
       coordinate y
    z : array
       coordinate y
    sets : List[tuple]
        set with node/element
    """
    __slots__ = ['_labels', '_system',  #'_number',
                 '_db', '_plane']  # '_jbc',

    def __init__(self, system: str) -> None:
        """
        system : cartesian/cylindrical/spherical
        """
        #self._labels : list = []
        #self._number : array = array('I', [])
        self._system = system  # get_coordinate_system(system)
        #
        #self._jbc: array = array('I', [])
        self._db = DBframework()

    #
    def __len__(self) -> int:
        return len(self._labels)

    def __iter__(self):
        """
        """
        nodes = sorted(self._labels)
        return iter(nodes)

    def __contains__(self, value) -> bool:
        return value in self._labels

    def __str__(self, units: str = "si") -> str:
        """ """
        lenght = ' m'
        space = " "
        output = "\n"
        output += "{:}\n".format(80 * "_")
        output += "\n"
        output += f"{33 * space}NODES\n"
        output += "\n"
        output += (f"Node {12 * space} x  [{lenght}] {4 * space} y  [{lenght}] {4 * space} z  [{lenght}] ")
        output += "\n"
        output += "{:}\n".format(80 * ".")
        output += "\n"
        #
        for key in self._labels:
            node = self.__getitem__(key)
            output += node.__str__()
        return output
        #

    #
    # ----------------------------------
    #
    @property
    def plane(self) -> NamedTuple:
        """ """
        return self._plane

    #
    @plane.setter
    def plane(self, plane: NamedTuple) -> None:
        """ """
        self._plane = plane

    #
    # ----------------------------------
    #
    @property
    def system(self) -> tuple:
        """
        """
        return self._system

    #@system.setter
    #def system(self, value:str) -> None:
    #    """
    #    """
    #    self._system = get_coordinate_system(value)
    #
    def _get_coordinates(self, coordinates):
        """ """
        if isinstance(coordinates, (list, tuple)):
            coordinates = check_point_list(coordinates, steps=3)
        elif isinstance(coordinates, dict):
            coordinates = check_point_dic(coordinates)
        else:
            raise Exception('Node input format not valid')
        return coordinates

    #
    #
    def get_new_point(self, coordinates: list | tuple):
        """ """
        #create a new point
        while True:
            node_number = next(self.get_number())
            try:
                self._labels.index(node_number)
            except ValueError:
                break
        self.__setitem__(node_number, coordinates)
        #
        try:
            node_id = coordinates.name
        except AttributeError:
            node_id = node_number
        #
        return node_id

    #
    def get_point_name(self, coordinates,
                       tol: float = 0.01, rel_tol: float = 1e-6) -> int:
        """
        tol: absolte tolerance in metres (0.010 m default)
        """
        # get index of x coord location in existing database
        coord = self._get_coordinates(coordinates)
        #
        items = self._isclose(key='*', item='x', value=coord[0],
                              abs_tol=tol, rel_tol=rel_tol)
        # check if y and z coord match
        if items:
            for item in items:
                if isclose(coord[1], item[4], abs_tol=tol, rel_tol=rel_tol):
                    if isclose(coord[2], item[5], abs_tol=tol, rel_tol=rel_tol):
                        return item[-1]
        raise IOError('   error coordinate not found')

    #
    def get_number(self, start: int = 1):
        """
        """
        try:
            n = max(self._labels) + 1
        except ValueError:
            n = start
        #
        while True:
            yield n
            n += 1
    #
    def orphan(self, elements):
        """
        Returns a list of the names of nodes that are not attached to any elements.
        """
        beams = elements.beam()
        connectivities = beams.get_connectivities()
        nodes_conn = get_nodes_connected(nodes=list(self._labels),
                                         connectivities=connectivities)
        single_nodes = [key for key, item in nodes_conn.items()
                        if not item]
        return single_nodes
    #
    # ----------------------------------
    #

#
#
#
def node_renumbering(nodes, elements):
    """ """
    # FIXME elements general not beams
    beams = elements.beam()
    connectivities = beams.get_connectivities()
    nodes_conn = get_nodes_connected(nodes=list(nodes.keys()),
                                     connectivities=connectivities)
    node_degrees = get_node_degrees(nodes_conn)
    new_nodes_number = get_node_renumber(nodes, node_degrees, nodes_conn)
    #single_nodes = [key for key, item in nodes_conn.items()
    #                if len(item) == 1]
    return new_nodes_number  #, single_nodes


#
def get_nodes_connected(nodes:list|tuple, connectivities:tuple|list):
    """
    nodes: list[int]
    connectivities [element_id, node_id, node_end]
    """
    # transpose connectivities [[element_id], [node_id], [node_end]]
    conn = list(map(list, zip(*connectivities)))
    # get elements end nodes
    memb_conn = {key: [conn[1][x]
                       for x, row in enumerate(conn[0])
                       if row == key]
                 for key in conn[0]}
    # group nodes connected each other.
    nodes_conn = defaultdict(list)
    for memb, items in memb_conn.items():
        for key in nodes:
            if key in items:
                nodes_conn[key].extend(items)
    # remove redundant items
    node_ends = defaultdict(list)
    for key, items in nodes_conn.items():
        temp = [item for item in items if item != key]
        node_ends[key] = list(set(temp))
    #
    # get nodes connected
    #nodes_conn = {}
    #for key in nodes:
    #    #
    #    nodes_conn[key] = [row for row in conn[1]
    #                       if row == key]
        #
        #node_ends = [[n for n, x in enumerate(col) if x == key]
        #             for col in connectivity]        
        #
        #flat_list = list(chain(*node_ends))
        #
        #nodes_conn[key] = [col[index] for col in connectivity
        #                   for index in flat_list if col[index] != key]
    #1 / 0
    return node_ends


#
def get_node_degrees(nodes_conn):
    """ """
    # Step 1 - get nodes degrees.
    degree = get_node_degree(nodes_conn)
    #shared_nodes = Counter({key: len(items)
    #                        for key, items in nodes_conn.items()})
    #
    #degree = [item[0] for item in shared_nodes.most_common()]
    #degree.reverse()
    #
    try:
        items = degree[0]
        raise IOError(f" orphan nodes {items}")
    except:
        # Step 2 - pick a starting node (i.e. the node with low degree).
        #first_degree = next(iter(degree.values()))
        #first_degree = next(iter(degree))
        #return first_degree
        return degree


#
#
def get_node_renumber(nodes, node_degrees: list[int],
                      nodes_conn: dict):
    """ """
    node_labels = [key for key in nodes.keys()]
    # Pick up one node from first level
    levels = []
    start_node = node_degrees[0]
    levels.append([start_node])
    # remove node from list
    node_labels = del_list_inplace(node_labels, levels[-1])
    # define level 2
    levels.append(nodes_conn[start_node])
    node_labels = del_list_inplace(node_labels, levels[-1])
    #
    # define next levels
    rem = list(chain(*levels))
    maxsteps = len(node_labels)
    step = 0
    while node_labels:
        step += 1
        rem = get_level(levels, nodes_conn, rem)
        del_list_inplace(node_labels, levels[-1])
        if step > maxsteps:
            raise RuntimeError(" Node renumbering fail")
    #
    new_nodes_number = list(chain(*levels))
    return new_nodes_number


#
#
def get_node_degree(nodes_conn:dict):
    """
    Scan all the nodes and order them according to their degree.
    The degree of a node is the number of nodes connected to it
    """
#    flat_nodes = list(chain.from_iterable(connectivities))
#    shared_nodes = Counter(flat_nodes)
#    degree = defaultdict(list)
#    for key, value in sorted(shared_nodes.items()):
#        degree.setdefault(value, []).append(key)
#    return degree
    shared_nodes = Counter({key: len(items)
                            for key, items in nodes_conn.items()})
    degree = [item[0] for item in shared_nodes.most_common()]
    degree.reverse()
    return degree
#

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
    _key: dict = {"number": r"\b(number|mesh)\b",
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
    _key: dict = {"number": r"\b(number|mesh)\b",
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
    name: int | str
    number: int
    index: int
    boundary: tuple | None
    system: str = "cartesian"

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
    #
    # ----------------------------------
    #
    def __str__(self) -> str:
        return "{:12d} {: 12.5f} {: 12.5f} {: 12.5f}\n" \
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
    # ----------------------------------
    #
    def distance(self, other):
        """ distance between two nodes"""
        #print("here")
        node1 = [self.x, self.y, self.z]
        return dist(node1[:3], other[:3])

    #
    # ----------------------------------
    #
    @property
    def fixity(self) -> str:
        """ """
        boundary = tuple(self.boundary[:6])
        match boundary:
            case (1, 1, 1, 1, 1, 1):
                return 'fixed'
            case (0, 0, 0, 0, 0, 0):
                return 'free'
            case (1, 1, 1, 1, 0, 0):
                return 'pinned'
            case (1, 0, 1, 1, 0, 0):
                return 'guide'
            case (1, 1, 0, 1, 0, 0):
                return 'guide'
            case (0, 1, 1, 1, 0, 0):
                return 'rolled'
            case _:
                return 'user_defined'
    #
    #def dof(self):
    #    """ """
    #    return self.index


class CoordCylindrical(NamedTuple):
    """
    """
    r: float
    theta: float
    z: float
    name: int | str
    number: int
    index: int
    boundaries: tuple
    system: str = "cylindrical"


class CoordSpherical(NamedTuple):
    """
    """
    r: float
    theta: float
    phi: float
    name: int | str
    number: int
    index: int
    boundaries: tuple
    system: str = "spherical"


#
#
def get_coord_system(system):
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
@dataclass
class NodePoint:
    name: int
    component: int
    number: int
    coord_system: str
    x: float | None
    y: float | None
    z: float | None
    r: float | None
    theta: float | None
    phi: float | None
    title: str | None
    index: int | None
    boundary: tuple | None

    #
    def system(self):
        """ """
        if 'cylindrical' in self.coord_system.lower():
            return CoordCylindrical
        elif 'spherical' in self.coord_system.lower():
            return CoordSpherical
        else:
            return CoordCartesian(x=self.x, y=self.y, z=self.z,
                                  name=self.name, number=self.number,
                                  index=self.index, boundary=self.boundary)


#
#
def check_point_list(data, steps: int = 6) -> list[float]:
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
def check_point_dic(data) -> list[float]:
    """ """
    new_data = [0, 0, 0]
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
