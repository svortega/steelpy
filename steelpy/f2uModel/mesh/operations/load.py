#
# Copyright (c) 2019-2021 steelpy
#

# Python stdlib imports
#from collections import defaultdict
#from collections.abc import Mapping
from typing import NamedTuple, Dict, List, Tuple, Union, Iterator
#import re


# package imports
from steelpy.process.math.vector import Vector
from steelpy.beam.static.singfunc import Trapezoidal, Point, Moment, SingFunction
#
#
class LineBeam(NamedTuple):
    """
    Local beam system
    """
    #qx1: float
    qy1: float
    qz1: float
    #
    #qx2: float
    qy2: float
    qz2: float
    #
    L1: float
    L2: float
    beam_length: float
    load_type: str = "line"
    #
    #
    def get_function(self):
        """ """
        output1 = Trapezoidal(self.qy1, self.qy2,
                              self.beam_length,
                              self.L1, self.L2)

        output2 = Trapezoidal(self.qz1, self.qz2,
                              self.beam_length,
                              self.L1, self.L2)

        return [output1, output2]
    #
    def __str__(self) -> str:
        """ """
        output = ""
        output += "Line      qy [N/m]    qz [N/m]    L [m  ]     End\n"
        output += "{:8s} {: 1.4E} {: 1.4E} {: 1.4E}  [1]\n".format(" ", self.qy1, self.qz1, self.L1)
        output += "{:8s} {: 1.4E} {: 1.4E} {: 1.4E}  [2]\n".format(" ", self.qy2, self.qz2, self.L2)
        return output
#
class PointBeam(NamedTuple):
    """
    """
    #fx: float
    fy: float
    fz: float
    L1:float
    name: Union[int,str]
    load_type: str = "point"
    #
    def __str__(self) -> str:
        """ """
        output = ""
        output += "Point     Fy [N  ]    Fz [N  ]    L1 [m  ]\n"
        output += "{:8s} {: 1.4E} {: 1.4E} {: 1.4E}\n".format(str(self.name), self.fy, self.fz, self.L1)
        return output
#
class MomentBeam(NamedTuple):
    """
    """
    #mx: float
    my: float
    mz: float
    L1:float
    name: Union[int,str]
    load_type: str = "moment"
    #
    def __str__(self) -> str:
        """ """
        output = ""
        output += "Moment    My [Nm ]    Mz [Nm ]    L1 [m  ]\n"
        output += "{:8s} {: 1.4E} {: 1.4E} {: 1.4E}\n".format(str(self.name), self.my, self.mz, self.L1)
        return output
#
#
#
class BeamLoad:

    __slots__ = ['_line', '_point', '_moment', '_torsion',
                 'beam_length']

    def __init__(self, L:float) -> None:
        """
        """
        self._point: List[ Tuple ] = [ ]
        self._line: List[ Tuple ] = [ ]
        self._moment: List[ Tuple ] = [ ]
        self._torsion: List[ Tuple ] = [ ]
        self.beam_length:float = L
    #
    @property
    def line(self):
        """ """
        return self._line

    @line.setter
    def line(self, load:List[int]):
        """ local system"""
        load.extend([ self.beam_length, "line" ])
        self._line.append(LineBeam._make(load))
    #
    def response(self, x: float,  E: float, Iy: float, Iz: float):
        """ """
        in_plane = Vector([0,0,0,0])
        out_plane = Vector([0,0,0,0])
        # line load
        for line in self._line:
            func = line.get_function()
            in_plane += func[0](x, E, Iy)
            out_plane += func[1](x, E, Iz)
        # point
        return [in_plane, out_plane]
    #
    def get_steps(self):
        """ """
        x_steps = [0,0.25,0.5,0.75,1]
        # line load
        for line in self._line:
            func = line.get_function()
            x_steps.extend(func[0].max_steps())
        # point
        #
        #if steps:
        x_steps = sorted(list(set(x_steps)))
        x_steps = [item for item in x_steps if item <= 1]
        #else:
        #    x_steps = [0,0.25,0.5,0.75,1]
        return x_steps
    #
    #
    #def get_response(self, x: float,  E: float, I: float):
    #    """ """
    #    plane = Vector([0,0,0,0])
    #    # line load
    #    for line in self._line:
    #        func = line.get_function()
    #        in_plane += func[0](x, E, Iy)
    #        out_plane += func[1](x, E, Iz)
    #    return plane
    #
    #
    # def __setitem__(self, load_name: Union[int, int],
    #                parameters: Union[List[float], Dict[str, float]]) -> None:
    #    """
    #    farg = [name, connectivity, material, section, type, group]
    #    """
    #    print('--')
    #
    #def __getitem__(self, load_name: Union[ str, int ]):
    #    """
    #    """
    #    return LoadType (load_name, cls=self )
    #
    ##
    #def __len__(self) -> float:
    #    """ """
    #    return len ( self._labels )
    #
    #def __iter__(self) -> Iterator:
    #    """ """
    #    for index, item in enumerate ( self._labels ):
    #        output = {"in_plane": SingFunction ( 0, 0 ),
    #                  "out_plane": SingFunction ( 0, 0 )}
    #        load = self._loads[ index ]
    #        if load.load_type == 'line':
    #            if load.qy1 or load.qy2:
    #                output[ "in_plane" ] = Trapezoidal ( load.qy1, load.qy2,
    #                                                     self.beam_length,
    #                                                     load.L1, load.L2 )
    #            if load.qz1 or load.qz2:
    #                output[ "out_plane" ] = Trapezoidal ( load.qz1, load.qz2,
    #                                                      self.beam_length,
    #                                                      load.L1, load.L2 )
    #            yield output
    #        elif load.load_type == 'point':
    #            if load.fy:
    #                output[ "in_plane" ] = Point ( load.fy,
    #                                               self.beam_length,
    #                                               load.L1 )
    #            if load.fz:
    #                output[ "out_plane" ] = Point ( load.fz,
    #                                                self.beam_length,
    #                                                load.L1 )
    #            yield output
    #        elif load.load_type == 'moment':
    #            if load.mz:
    #                output[ "in_plane" ] = Moment ( load.mz,
    #                                                self.beam_length,
    #                                                load.L1 )
    #            if load.my:
    #                output[ "out_plane" ] = Moment ( load.my,
    #                                                 self.beam_length,
    #                                                 load.L1 )
    #            yield output
    #        else:
    #            raise IOError ( "wrong type" )
    #    # return iter(self._labels)
    #
    #def __contains__(self, value) -> bool:
    #    return value in self._labels
    #
    ##
    #def __str__(self) -> str:
    #    """ """
    #    output = "------- Loading\n"
    #    output += "\n"
    #    for x, load_name in enumerate ( self._labels ):
    #        index = self._labels.index ( load_name )
    #        load = self._loads[ index ]
    #        # output += "Support {:} : {:}\n".format(support_name, self._fixity[x])
    #        output += load.__str__ ()
    #        output += "\n"
    #    return output
    ##


class LoadType:
    __slots__ = [ 'cls']

    def __init__(self, cls):
        """
        """
        self.cls = cls

    #
    #
    # def __getitem__(self, load_name: Union[str,int]):
    #    """
    #    """
    #    print('--')
    #
    @property
    def point(self):
        """ """
        index = self.cls._labels.index ( self._load_name )
        return self.cls._loads[ index ]
        # print('--')

    @point.setter
    def point(self, load: Union[ List, Dict ]):
        """ """
        load = self._get_point_load ( load )
        load.extend ( [ self._load_name, "point" ] )
        self.cls._loads.append ( PointBeam._make ( load ) )
        self.cls._labels.append ( self._load_name )

    #
    #
    @property
    def moment(self):
        """ """
        index = self.cls._labels.index ( self._load_name )
        return self.cls._loads[ index ]

    @moment.setter
    def moment(self, load: Union[ List, Dict ]):
        """ """
        load = self._get_moment_load ( load )
        load.extend ( [ self._load_name, "moment" ] )
        self.cls._loads.append ( MomentBeam._make ( load ) )
        self.cls._labels.append ( self._load_name )
        #

    #
    @property
    def line(self):
        """ """
        index = self.cls._labels.index ( self._load_name )
        return self.cls._loads[ index ]

    @line.setter
    def line(self, load:List[int]):
        """ local system"""
        #load = self._get_line_load(load)
        #load.extend ( [ self._load_name, "line" ] )
        self.cls._loads.append ( LineBeam._make ( load ) )
        self.cls._labels.append ( self._load_name )

    #
    #
    @property
    def torsion(self):
        """ """
        index = self.cls._labels.index ( self._load_name )
        load = self.cls._loads[ index ]
        print ( '--' )

    #