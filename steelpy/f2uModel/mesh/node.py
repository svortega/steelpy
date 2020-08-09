# 
# Copyright (c) 2009-2019 fem2ufo
#


# Python stdlib imports
from array import array
from collections.abc import Mapping
import logging
from math import isclose
from typing import NamedTuple, Tuple, List, Iterator, Iterable, Union, Dict
import re

# package imports
from steelpy.process.units.main import Units
#from fem2ufo.f2u_model.femodel.sets.sets import MeshGroupCases


class CoordCartesian(NamedTuple):
    """ Cartesian coordinate system"""
    x: Units
    y: Units
    z: Units
    number: int
    name:Union[str,int]
    system:str ="cartesian"
    #boundaries: Iterable
    #sets: List[Tuple]
    #
    # def add_element(self, element_number):
    # """
    # """
    # self.elements.append(element_number)
    #
    def get_coordinates(self, units:str='metre'):
        """
        """
        return "{:14.0f} {: 14.5f} {: 14.5f} {: 14.5f}".format(self.number,
                                                               self.x.convert(units).value,
                                                               self.y.convert(units).value,
                                                               self.z.convert(units).value)

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
        return "{:14.0f} {: 14.5f} {: 14.5f} {: 14.5f}".format(self.number, self.x, self.y, self.z)

    def __eq__(self, other) -> bool:
        """
        """
        if (isclose(self.x.value, other.x.value, abs_tol=1e-03)
                and isclose(self.y.value, other.y.value, abs_tol=1e-03)
                and isclose(self.z.value, other.z.value, abs_tol=1e-03)):
            return True
        return False


class CoordCylindrical(NamedTuple):
    """
    """
    r: float
    theta: float
    z: float
    number: int
    name:Union[str,int]
    system:str ="cylindrical"

class CoordSpherical(NamedTuple):
    """
    """
    r: float
    theta: float
    phi: float
    number: int
    name:Union[str,int]
    system:str ="spherical"
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
class Nodes(Mapping):
    """
    This is a fem2ufo model node class


    Parameters
    ----------
    boundaries: object
        f2u boundary object


    Attributes
    ----------
    _labels : array
        node internal number
    _x : array
        coordinate x
    _y : array
        coordinate y
    _z : array
        coordinate y
    _sets : List[tuple]
        set with node/element
    """
    __slots__ = ['_labels', '_number', '_x', '_y', '_z', 
                 '_system',  '_sets', 'f2u_units']
    
    def __init__(self) -> None:
        """
        """
        global f2u_units
        f2u_units = Units()
        #self._boundaries = boundaries
        self.system: str = 'cartesian'
        # set variables
        self._labels: List[Union[str,int]] = []
        self._number : List[int] = array('I', [])
        self._x: List[float] = array('f', [])
        self._y: List[float] = array('f', [])
        self._z: List[float] = array('f', [])
    #
    #
    # ---------------------------------
    #
    def __setitem__(self, node_name: Union[int, str],
                    coordinates: Union[List[float], Dict[str, float]]) -> None:
        """
        """
        try:
            self._labels.index(node_name)
            raise Exception('    *** warning point {:} already exist'
                            .format(node_name))
        except ValueError:
            coordinates = self._get_coordinates(coordinates)
            #
            self._labels.append(node_name)
            self._number.append(self._labels.index(node_name))
            self._x.append(coordinates[0])
            self._y.append(coordinates[1])
            self._z.append(coordinates[2])            

    def __getitem__(self, node_name: int) -> Tuple:
        """
        node_name : node number
        """
        try:
            _index = self._labels.index(node_name)
            return self.system(x=self._x[_index] * f2u_units.m, 
                               y=self._y[_index] * f2u_units.m, 
                               z=self._z[_index] * f2u_units.m,
                               number=self._number[_index], name=node_name)
        except ValueError:
            raise IndexError('   *** node {:} does not exist'.format(node_name))
    #
    def __delitem__(self, node_name: int) -> None:
        """
        """
        try:
            i = self._labels.index(node_name)
            #self._number.remove(node_name)
            self._number.pop(i)
            self._labels.pop(i)
            self._x.pop(i)
            self._y.pop(i)
            self._z.pop(i)
            #self._sets.pop(i)
        except IndexError:
            logging.warning(' delete -- node {:} does not exist'.format(node_name))
            return

    def __len__(self) -> float:
        return len(self._labels)

    def __iter__(self) -> Iterator:
        """
        """
        return iter(self._labels)

    def __contains__(self, value) -> bool:
        return value in self._labels
    
    @property
    def system(self) -> Tuple:
        """
        """
        return self._system
    
    @system.setter
    def system(self, value:str) -> None:
        """
        """
        #if 'cylindrical' in value:
        #    self._system = CoordinatesCylindrical
        #
        #elif 'spherical' in value:
        #    self._system = CoordinatesSpherical
        #
        #else:
        #    self._system = CoordinatesCartesian
        self._system = get_coordinate_system(value)
    #
    #
    #def scale(self, scalar):
    #    """Return the product of self and numeric object alpha."""
    #    if not isinstance(scalar, (float, int)):
    #        raise TypeError("must be a scalar")
    #    self._x = array('f',[x * scalar for x in self._x])
    #    self._y = array('f',[x * scalar for x in self._y])
    #    self._z = array('f',[x * scalar for x in self._z])
    #    #return array('f',result)
    #
    #
    @property
    def _get_maxmin(self):
        """
        """
        max_x = max(self._x)
        min_x = min(self._x)
        #
        max_y = max(self._y)
        min_y = min(self._y)
        #
        max_z = max(self._z)
        min_z = min(self._z)
        return [max_x, max_y, max_z], [min_x, min_y, min_z]
    #
    def _get_coordinates(self, coordinates):
        """ """
        if isinstance(coordinates, (list, tuple)):
            coordinates = check_point_list(coordinates, steps=3)
        elif isinstance(coordinates, dict):
            coordinates = check_point_dic(coordinates)
        else:
            raise Exception('   *** Node input format not recognized')
        return coordinates
    #
    def get_point_name(self, coordinates, tol:float=0.01):
        """ 
        tol: absolte tolerance in metres (0.010 m default)
        """
        # check if point already exist
        if isinstance(coordinates, self.system):
            return coordinates.name
        # get index of x coord location in existing database
        coord = self._get_coordinates(coordinates)
        indeces = [index for index, item in enumerate(self._x)
                   if isclose(coord[0], item, abs_tol=tol)]
        # check if y and z coord match 
        if indeces:
            for index in indeces:
                if isclose(coord[1], self._y[index], abs_tol=tol):
                    if isclose(coord[2], self._z[index], abs_tol=tol):
                        return self._labels[index]
        #create a new point
        while True:
            node_name = "pnt_{:}".format(str(next(self.get_number())))
            try:
                self._labels.index(node_name)
            except ValueError:
                break
        self.__setitem__(node_name, coordinates)
        return node_name
    #
    def get_number(self, start:int=1)-> Iterable[int]:
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