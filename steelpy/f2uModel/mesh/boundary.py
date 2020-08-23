# 
# Copyright (c) 2009-2020 fem2ufo
#

# Python stdlib imports
from array import array
from collections.abc import Mapping
import re
from typing import NamedTuple, Tuple, Union, List, Dict


# package imports

# -----------------------
#
class BoundaryItem(NamedTuple):
    """
    """
    x: float
    y: float
    z: float
    rx: float
    ry: float
    rz: float
    number: int
    name: Union[str, int]

    def __str__(self) -> str:
        return "{: 14.5f} {: 14.5f} {: 14.5f} {: 14.5f} {: 14.5f} {: 14.5f}".format(self.x, self.y, self.z,
                                                                                    self.rx, self.ry, self.rz)


#
class BoundaryNodes(Mapping):
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
                 '_labels', '_number']

    def __init__(self) -> None:
        """
        """
        self._number : array = array('I', [])
        self._x : array = array('f', [])
        self._y : array = array('f', [])
        self._z : array = array('f', [])
        self._rx : array = array('f', [])
        self._ry : array = array('f', [])
        self._rz : array = array('f', [])
        #
        self._labels: List[Union[int,str]] = [] #  = array('I', [])
        #self._type : List[Union[str, int]] = [] # array('I', [])
        #
        # fix_type = ['release', 'gap', 'prescribed', 'dependence',
        #            'link', 'master', 'rigid', 'constrain']        

    def __setitem__(self, node_number: int,
                    value:Union[List, Tuple, Dict, str]) -> None:
        """
        """
        try:
            # TODO : update data
            self._labels.index(node_number)
            raise Warning('    *** warning node {:} already exist'.format(node_number))
        except ValueError:
            self._labels.append(node_number)
            self._number.append(self._labels.index(node_number))
            #
            if isinstance(value, str):
                if re.match(r"\b(fix(ed)?)\b", value, re.IGNORECASE):
                    #self._type.append('fixed')
                    value = [1,1,1,1,1,1]
                elif re.match(r"\b(pinn(ed)?|roll)\b", value, re.IGNORECASE):
                    #self._type.append('pinned')
                    value = [1,1,1,1,0,0]
                else:
                    raise IOError("boundary type {:} not implemented".format(value))
                # update
                self._x.append(value[0])
                self._y.append(value[1])
                self._z.append(value[2])
                self._rx.append(value[3])
                self._ry.append(value[4])
                self._rz.append(value[5])
            else:
                if isinstance( value, (list, tuple) ):
                    self._x.append( value[ 0 ] )
                    self._y.append( value[ 1 ] )
                    self._z.append( value[ 2 ] )
                    self._rx.append( value[ 3 ] )
                    self._ry.append( value[ 4 ] )
                    self._rz.append( value[ 5 ] )
                elif isinstance( value, dict ):
                    self._x.append( value[ 'x' ] )
                    self._y.append( value[ 'y' ] )
                    self._z.append( value[ 'z' ] )
                    self._rx.append( value[ 'rx' ] )
                    self._ry.append( value[ 'ry' ] )
                    self._rz.append( value[ 'rz' ] )
                else:
                    print( '   *** error node input format not recognized' )
                #
                #self._type.append('user')
                # _type = _bound_type[_type]

    #
    def __getitem__(self, node_number: Union[None,str,int]) -> Union[Tuple,bool]:
        """
        """
        try:
            _index = self._labels.index(node_number)
            return BoundaryItem(x=self._x[_index], y=self._y[_index], z=self._z[_index],
                                rx=self._rx[_index], ry=self._ry[_index], rz=self._rz[_index],
                                number=self._number[_index], name=node_number)
        except ValueError:
            return False
            # raise IndexError

    #
    def __iter__(self):
        """
        """
        return iter(self._labels)
        # for node_number in self._labels:
        #    yield self.__getitem__(node_number)

    #
    def __delitem__(self, node_number: int) -> None:
        """
        """
        try:
            i = self._labels.index(node_number)
            self._labels.remove(node_number)
            self._number.pop(i)
            #self._type.pop(i)
            self._x.pop(i)
            self._y.pop(i)
            self._z.pop(i)
            self._rx.pop(i)
            self._ry.pop(i)
            self._rz.pop(i)
        except IndexError:
            #logging.warning('  ** boundary {:} does not exist'.format(node_number))
            raise Warning('  ** boundary {:} does not exist'.format(node_number))

    #
    def __contains__(self, value) -> bool:
        return value in self._labels

    def __len__(self) -> float:
        return len(self._labels)
#
#
class Boundaries:
    
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
