# 
# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
from array import array
import logging
from math import isclose, dist
#from collections.abc import Mapping
#
# package imports
from steelpy.ufo.utils.node import NodePoint, NodeBasic
#
#
#
class NodesIM(NodeBasic):
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
                 '_system',  '_sets', '_component', '_boundary']
    
    def __init__(self, component: str, #boundary, 
                 system:str = 'cartesian') -> None:
        """
        """
        super().__init__(name=component, system=system)
        self._labels : list = []
        #
        self._x: array = array('f', [])
        self._y: array = array('f', [])
        self._z: array = array('f', [])
        self._number : array = array('I', [])
        #
        self._component = component
        #self._boundary = boundary
    #
    #
    #
    # ----------------------------------
    #
    def __setitem__(self, name: int|str,
                    coordinates: list[float]|dict[str,float]) -> None:
        """
        """
        try:
            self._labels.index(name)
            raise Exception(f'*** warning point {name} already exist')
        except ValueError:
            coordinates = self.get_coordinates(coordinates)
            self._labels.append(name)
            self._number.append(len(self._labels))
            self._x.append(coordinates[0])
            self._y.append(coordinates[1])
            self._z.append(coordinates[2])
            #
            #print('-- node')
            #self._boundary[node_name] = [0,0,0,0,0,0]

    def __getitem__(self, node_name: int) -> tuple:
        """
        node_name : node number
        """
        try:
            index = self._labels.index(node_name)
            node = NodePoint(name=node_name,
                             component=self._component, 
                             number=self._number[index], 
                             coord_system=self._system, 
                             x=self._x[index], 
                             y=self._y[index], 
                             z=self._z[index],
                             r=None, theta=None, phi=None,
                             title=node_name, 
                             index=self._number[index]-1,
                             boundary=None)
            return node.system()
        except ValueError:
            raise IndexError('   *** node {:} does not exist'.format(node_name))
    
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
    #
    # ----------------------------------
    #
    def renumbering(self, new_numbers:list[int]):
        """ """
        indexes = [self._labels.index(node_name) 
                   for node_name in new_numbers]
        #for x, index in enumerate(indexes):
        #    self._number[index] = x+1
        self._renumber(indexes)
        #print('')
    #
    def _renumber(self, indexes:list[int]):
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
    def get_maxmin(self):
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
    #  
    #
    def _isclose(self, item:str, value:float, key:str = '*',
                 rel_tol:float=1e-6, abs_tol:float=0.0)-> tuple:
        """ """
        if item in ['x']:
            val = self._x
            
        elif item in ['y']:
            val = self._y
        
        elif item in ['z']:
            val = self._z
        
        else:
            raise IOError(f' coord:{item} not valid')
        
        #
        idx = [x for x, item in enumerate(val)
                if isclose(value, item, abs_tol=abs_tol, rel_tol=rel_tol)]
        #
        items = [[self._number[i], self._number[i], self._system, 
                  self._x[i], self._y[i],  self._z[i], self._labels[i]]
                 for i in idx]
        return items
    #
    # ----------------------------------
#
