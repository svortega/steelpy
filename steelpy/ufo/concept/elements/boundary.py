# 
# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
from array import array
from collections.abc import Mapping
import re
from typing import NamedTuple
#
#
# package imports
from steelpy.ufo.process.elements.nodes import NodePoint
from steelpy.ufo.process.elements.boundary import BoundaryNode, BoundaryItem

#
class BoundaryNodes(BoundaryNode):
    """
    FE Fixity

    Boundary
        |_ type
        |_ constrain [x, y, z, mx, my, mz]

    **Parameters**:
      :name : coupling, pile
      :type : free, fix, dependent, prescribed, supernode, master
      :constrain :
            0- free
            1- fixed
            2- prescribed displacement, temperature, different from zero
            3- linear dependent
            4- supernode (link)
    """
    __slots__ = ['_x', '_y', '_z',
                 '_rx', '_ry', '_rz',
                 '_labels', '_number', '_title']

    def __init__(self) -> None:
        """
        """
        self._labels: list[int|str] = []
        self._number : array = array('I', [])
        self._x : array = array('i', [])
        self._y : array = array('i', [])
        self._z : array = array('i', [])
        self._rx : array = array('i', [])
        self._ry : array = array('i', [])
        self._rz : array = array('i', [])
        #
        self._title : list[int|str] = []
        #
        super().__init__()
        #

    def __setitem__(self, node_id:int|str,
                    value:list|tuple|dict|str) -> None:
        """
        1 : fix
        0 : free
        """
        try:
            # TODO : update data
            self._labels.index(node_id)
            raise Warning('    *** warning node {:} already exist'.format(node_id))
        except ValueError:
            #if isinstance(value, dict):
            #    self._x.append(value['x'])
            #    self._y.append(value['y'])
            #    self._z.append(value['z'])
            #    self._rx.append(value['rx'])
            #    self._ry.append(value['ry'])
            #    self._rz.append(value['rz'])
            #    try:
            #        title = value['title']
            #    except KeyError:
            #        title = "NULL"
            #else:
            title = self.setup_data(value)
            #
            self._labels.append(node_id)
            self._number.append(self._labels.index(node_id))
            self._title.append(title)
            # _type = _bound_type[_type]
    #
    #
    def __getitem__(self, node_id: None | str | int) -> tuple | bool:
        """
        """
        try:
            _index = self._labels.index(node_id)
            return BoundaryItem(x=self._x[_index], y=self._y[_index], z=self._z[_index],
                                rx=self._rx[_index], ry=self._ry[_index], rz=self._rz[_index],
                                number=self._number[_index], name=self._title[_index],
                                node=node_id)
        except ValueError:
            return False
            # raise IndexError

    #
    #
    def __delitem__(self, node_id: int) -> None:
        """
        """
        try:
            i = self._labels.index(node_id)
            self._labels.remove(node_id)
            self._number.pop(i)
            # self._type.pop(i)
            self._x.pop(i)
            self._y.pop(i)
            self._z.pop(i)
            self._rx.pop(i)
            self._ry.pop(i)
            self._rz.pop(i)
        except IndexError:
            # logging.warning('  ** boundary {:} does not exist'.format(node_id))
            raise Warning('  ** boundary {:} does not exist'.format(node_id))

    #
    #
    #
    def setup_data(self, value):
        """ """
        title = None
        value = self._get_fixity(value)
        # update data
        self._x.append(value[0])
        self._y.append(value[1])
        self._z.append(value[2])
        self._rx.append(value[3])
        self._ry.append(value[4])
        self._rz.append(value[5])
        return title
    #
    #
    def update_data(self, value:list|tuple):
        """ """
        title = "NULL"
        if isinstance(value, list):
            node_name = value[0]
            if isinstance(value[1], str):
                value = self.get_boundary(value[1])
            else:
                1/0
        elif isinstance(value, BoundaryItem):
            node_name = value.node
            title = value.name
            value = value[:6]

        else:
            1/0
        #
        index = self._labels.index(node_name)
        number = self._number[index]
        # delete
        self.__delitem__(node_name)
        # update
        self._labels.insert(index, node_name)
        self._number.insert(index, number)
        self._title.insert(index, title)
        self._x.insert(index, value[0])
        self._y.insert(index, value[1])
        self._z.insert(index, value[2])
        self._rx.insert(index, value[3])
        self._ry.insert(index, value[4])
        self._rz.insert(index, value[5])
    #
    #
    def get_number(self, start: int = 1):
        """
        """
        try:
            n = len(self._labels) + 1
        except ValueError:
            n = start
        #
        while True:
            yield n
            n += 1
#
#
# TODO : remove redundant code
#
class BoundaryJoint:

    def __init__(self) -> None:
        """
        """
        self._nodes = BoundaryNodes()
    #
    #
    def supports(self, values:list|None=None):
        """"""
        return self.point(values)    
    #
    @property
    def point(self):
        """"""
        return self._nodes

    @point.setter
    def point(self, values):
        """"""
        for value in values:
            self._nodes[value[0]] = value[1:]


#
#
class BoundaryType:
    __slots__ = ["_cls", "_boundary_name"]
    
    def __init__(self, cls, boundary_name):
        """
        """
        self._cls = cls
        self._boundary_name = boundary_name
    #
    @property
    def support(self):
        """ """
        return self._cls._nodes.point[self._boundary_name]
    
    @support.setter
    def support(self, conditions):
        """ """
        self._cls._nodes.point[self._boundary_name] = conditions
        #print('--')
    #
    @property
    def point(self):
        """ """
        index = self._cls._labels.index(self._boundary_name)
        return self._cls._points[index]
#
#
class ConceptBoundaries:
    
    __slots__ = ["_labels", "_nodes",
                 "f2u_points", "_points"]
    
    def __init__(self):
        """
        """
        self._nodes = BoundaryJoint()
        self._labels: list[str|int] = []
        self._points: list[tuple[float]] = []
    
    def __setitem__(self, support_name: int|str,
                    coordinates: list[float]|dict[str, float]) -> None:
        """
        """
        try:
            self._labels.index(support_name)
            raise Exception('boundary name {:} already exist'.format( support_name))
        except ValueError:
            self._labels.append(support_name)
            try:
                self._points.append((coordinates[0],
                                     coordinates[1],
                                     coordinates[2]))
            except IndexError:
                self._points.append((coordinates[0],
                                     coordinates[1], 0))
    
    def __getitem__(self, support_name:int) -> tuple:
        """
        node
        """
        return BoundaryType(cls=self, boundary_name=support_name)
    
    #
    #
    def __len__(self) -> float:
        return len(self._labels)

    def __iter__(self) -> Iterator:
        """
        """
        return iter(self._labels)

    def __contains__(self, value) -> bool:
        return value in self._labels
    #
    #
    def df(self, df, columns:dict|None=None):
        """ """
        self._labels = df.name.tolist()
        self._points = df[["x", "y", "z"]].values.tolist()
        points = df[["name", "support"]].values.tolist()
        self._supports.node = points
        #print('---')
#
# 
class BoundaryConcept:
    """
    Boundary Condition :
    
    Point : free/fixed/supernode/manual
    Line
    
    """
    __slots__ = ['_labels', '_number', '_supports',
                 '_component']
    
    def __init__(self, points, component:str|int):
        """
        """
        self._component = component
        self._supports = BoundariesSupports(points,
                                            self._component)
    #
    #
    #
    #def __len__(self) -> float:
    #    return len(self._labels)
    #def __iter__(self):
    #    """
    #    """
    #    return iter(self._labels)
    #def __contains__(self, value) -> bool:
    #    return value in self._labels    
    #
    #
    def support(self):
        """ """
        return self._supports
    #
    #def support(self, dof:list|tuple|str = 'fixed',
    #            name:str|None = None):
    #    """Boundary condition along an edge"""
    #    bname = name
    #    mnumber = next(self.get_number())
    #    if not bname:
    #        bname = f"sp_{mnumber}"        
    #    try:
    #        self._labels.index(bname)
    #        raise IOError(f'    *** warning support {bname} already exist')
    #    except ValueError:
    #        self._labels.append(bname)
    #        self._number.append(mnumber)
    #        fixity = self.get_fixity(dof)
    #        self._nodes[bname] = BoundariesJoint(fixity=fixity)
    #        return self._nodes[bname]
    #
    #def point(self, dof:list|tuple|str):
    #    """
    #    Boundary condition inserted at points
    #    """
    #    self._fixity = self.get_fixity(dof)
    #    #print('-->', fixity)
    #    return self._nodes
    #
    # ----------------------------
    # Operations
    # ----------------------------
    #
    #
    #def get_number(self, start:int=1):
    #    """
    #    """
    #    try:
    #        n = max(self._number) + 1
    #    except ValueError:
    #        n = start
    #    #
    #    while True:
    #        yield n
    #        n += 1
    #
    #
    #
    def df(self, df, columns:dict|None=None):
        """ """
        # types: support/line 
        grptype = df.groupby(['type'])
        #
        # support
        #
        support = grptype.get_group('support')
        grpboundary = support.groupby(['boundary'])
        for key, items in grpboundary:
            boundary = key[0]
            # FIXME: in case input is a list
            if boundary.lower() == 'free':
                continue
            self._supports[boundary] = boundary
            for item in items.itertuples():
                self._supports[boundary].points = item.x, item.y, item.z
        #
        # Lines 
        #print('---')    
#
#
class BoundariesSupports(Mapping):
    """
    """
    __slots__ = ['_labels', '_type', '_number', '_fixity',
                 '_nodes', '_line', '_nitems', '_litems',
                 '_points', '_component']
    
    def __init__(self, points, component:str|int):
        """
        """
        # concept points
        self._points = points
        self._component = component
        #
        self._fixity:list = []
        self._labels:list[str|int] = []
        self._type:list[str|int] = []
        self._number:list[int] = []
        #
        self._nodes = BoundaryNodes()
        self._line = []
        #
        self._nitems: dict = {}
        self._litems: dict = {}        
    #
    def __setitem__(self, name: int|str,
                    dof:list|tuple|str) -> None:
        """
        """
        try:
            self._labels.index(name)
            raise IOError(f'    *** warning support {name} already exist')
        except ValueError:
            self._labels.append(name)
            self._nitems[name] = []
            self._litems[name] = []           
        #
        self._fixity.append(self.get_fixity(dof))
        #print('-->')
        
    #
    def __getitem__(self, name: int|str):
        """
        """
        try:
            index = self._labels.index(name)
            return SupportItems(cls=self, name=name)
        except ValueError:
            raise IndexError(f' ** Support {name} no valid')
        
        
    #
    def __len__(self) -> float:
        return len(self._labels)
    
    def __iter__(self):
        return iter(self._labels)
    
    def __contains__(self, value) -> bool:
        return value in self._labels    
    #    
    #
    # ----------------------------
    # Operations
    # ----------------------------
    #    
    #
    def get_fixity(self, fixity):
        """ """
        if isinstance(fixity, str):
            if re.match(r"\b(fix(ed)?)\b", fixity, re.IGNORECASE):
                return [1,1,1,1,1,1]
            elif re.match(r"\b(pinn(ed)?|roll)\b", fixity, re.IGNORECASE):
                return [1,1,1,0,0,0]
            elif re.match(r"\b(free)\b", fixity, re.IGNORECASE):
                return None
            else:
                raise IOError("boundary type {:} not implemented".format(fixity))
        elif isinstance(fixity, (list, tuple)):
            return fixity
        elif isinstance(fixity, dict):
            return [fixity['x'], fixity['y'], fixity['z'], 
                    fixity['rx'], fixity['ry'], fixity['rz']]
        else:
            raise Exception('   *** Boundary input format not recognized')
    #
#
#
class SupportItems:
    """
    """
    __slots__ = ['cls', '_name']
    
    def __init__(self, cls, name):
        """
        """
        self._name = name
        self.cls = cls
    # 
    #
    def line(self):
        """Boundary condition along an edge"""
        raise NotImplemented()
    #
    #
    @property
    def points(self):
        """ """
        name = self._name
        index = self.cls._labels.index(name)
        fixity =  self.cls._fixity[index]
        #number = self.cls._number[index]
        items = self.cls._nitems[name]
        return PointItem(*fixity, name, items)
    
    @points.setter
    def points(self, coordinates:tuple|list):
        """Boundary condition inserted at support points"""
        index = self.cls._labels.index(self._name)
        #
        #try:
        #    pname = self.cls._points.get_point_name(coordinates)
        #except OSError:
        #    pname = self.cls._points.get_new_point(coordinates)
        #else:
        #    raise IOError('Point data not valid')
        #
        try:
            nname = coordinates.name
            #npoint = coordinates
        except AttributeError:
            nodes = self.cls._nodes
            while True:
                nname = next(nodes.get_number())
                if nname in nodes._labels:
                    continue
                break
            # create a temporary node
            coord = self.cls._nodes._get_coordinates(coordinates)
            pname = self.cls._points.get_point_name(coordinates)
            
            bound = BoundaryItem(*self.cls._fixity[index],
                                 name=self._name,
                                 points=[nname])
            coord = NodePoint(name=nname,
                              component=self.cls._component,
                              number=nname,
                              coord_system="cartesian",
                              x=coord[0],y=coord[1],z=coord[2],
                              theta=None, phi= None, r=None,
                              title=None, index=None,
                              boundary=bound)
            coordinates = coord.system()
        #
        #nname = self.cls._points.get_name(pname)
        #npoint = self.cls._points[nname]
        #
        #npoint = self.cls._nodes[nname]
        #npoint = coordinates
        self.cls._nodes[nname] = self.cls._fixity[index]
        self.cls._nitems[self._name].append(coordinates)
        #print('bound')
    #
    #
    def _set_item(self, b_name, b_type):
        """ """
        try:
            self._labels.index(b_name)
            raise IOError('   error boundary {:} already exist'.format(b_name))
        except ValueError:
            boundary_type = b_type
            self._labels.append(b_name)
            self._type.append(boundary_type)
            b_number = next(self.get_number())
            self._number.append(b_number)
        #
        return b_number
    #
    def get_number(self, start:int=1):
        """
        """
        try:
            n = max(self._number) + 1
        except ValueError:
            n = start
        #
        while True:
            yield n
            n += 1
    #
    #
    def __str__(self, units:str="si") -> str:
        """ """
        lenght = ' m'
        space = " "
        output = "\n"
        output += "{:}\n".format(80*"_")
        output += "\n"
        output += f"{33*space}BOUNDARIES\n"
        output += "\n"
        output += (f"Node {14*space} x {6*space} y {6*space} z {5*space} mx {5*space} my {5*space} mz {5*space} title")
        output += "\n"
        output += "{:}\n".format(80*".")
        output += "\n"
        for key, node in self._nodes.items():
            #if sum(node[:6]) == 0:
            #    continue
            output += node.__str__()
        return output
#
#
class PointItem(NamedTuple):
    """
    """
    x: float
    y: float
    z: float
    rx: float
    ry: float
    rz: float
    #number: int
    name: str|None
    points:list
    
    #def __str__(self) -> str:
    #    if (name := self.name) == None:
    #        name = ""
    #    return "{:12d} {: 8.0f} {: 8.0f} {: 8.0f} {: 8.0f} {: 8.0f} {: 8.0f} {:>12s}\n"\
    #        .format(self.node, self.x, self.y, self.z, self.rx, self.ry, self.rz, name)
#
#

