# 
# Copyright (c) 2009-2021 fem2ufo
#


# Python stdlib imports
from array import array
from collections.abc import Mapping
import logging
from math import isclose, dist
from typing import NamedTuple, Tuple, List, Iterator, Iterable, Union, Dict
import re

# package imports
from steelpy.process.units.main import Units
from steelpy.f2uModel.mesh.operations.nodes import (check_point_list, check_point_dic, 
                                                    get_coordinate_system)

#
#
class NodesInMemory(Mapping):
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
        #self._labels: List[Union[str,int]] = []
        self._labels : array = array('I', [])
        self._number : array = array('I', [])
        self._x: array = array('f', [])
        self._y: array = array('f', [])
        self._z: array = array('f', [])
    #
    # ---------------------------------
    #
    def __setitem__(self, node_name: int,
                    coordinates: Union[List[float], Dict[str, float]]) -> None:
        """
        """
        try:
            self._labels.index(node_name)
            raise Exception('    *** warning point {:} already exist'
                            .format(node_name))
        except ValueError:
            coordinates = self._get_coordinates(coordinates)
            self._labels.append(node_name)
            self._number.append(len(self._labels))
            self._x.append(coordinates[0])
            self._y.append(coordinates[1])
            self._z.append(coordinates[2])            

    def __getitem__(self, node_name: int) -> Tuple:
        """
        node_name : node number
        """
        try:
            _index = self._labels.index(node_name)
            return self.system(x=self._x[_index], 
                               y=self._y[_index], 
                               z=self._z[_index],
                               name=node_name, 
                               number=self._number[_index], 
                               index=self._number[_index]-1)
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
    def renumbering(self, new_numbers:List[int]):
        """ """
        indexes = [self._labels.index(node_name) 
                   for node_name in new_numbers]
        #for x, index in enumerate(indexes):
        #    self._number[index] = x+1
        self._renumber(indexes)
        #print('')
    #
    def _renumber(self, indexes:List[int]):
        """
        """
        #lst_size = len(self._labels)
        #_index = [self._number.index(i+1)
        #          for i in range(lst_size)]
        #
        self._number = [self._number[indx] for indx in indexes]
        self._labels = [self._labels[indx] for indx in indexes]
        self._x = [self._x[indx] for indx in indexes]
        self._y = [self._y[indx] for indx in indexes]
        self._z = [self._z[indx] for indx in indexes]
        #print('-->')
    #
    #def update_number(self, node_name:int, value:int):
    #    """ """
    #    _index1 = self._labels.index(node_name)
    #    temp1 = self._number[_index1]
    #    _index2 = self._number.index(value)
    #    #memb = self._labels[_index2]
    #    self._number[_index1] = value
    #    self._number[_index2] = temp1
    #    #1/0
    #    #print('---', self._number[_index1], self._number[_index2])
    #
    def update_item(self, node_name:int, item:str, value:Union[float,int]):
        """ """
        1/0
        _items = {"number":self._number}
        try:
            _index = self._labels.index(node_name)
            _items[item][_index] = value
        except ValueError:
            raise IndexError('   *** node {:} does not exist'.format(node_name))        
        #print('---')
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
        raise IOError('   error coordinate not found')
    #
    def get_new_point(self, coordinates):
        """ """
        #create a new point
        while True:
            #node_name = "pnt_{:}".format(str(next(self.get_number())))
            node_name = next(self.get_number())
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
            n = max(self._labels) + 1
        except ValueError:
            n = start
        #
        while True:
            yield n
            n += 1
#
