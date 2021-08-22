# 
# Copyright (c) 2019-2021 steelpy
# 

# Python stdlib imports
from collections import defaultdict
from collections.abc import Mapping
from typing import NamedTuple, Dict, List, Tuple, Union, Iterator
import re


# package imports
from steelpy.beam.static.singfunc import Trapezoidal, Point, Moment, SingFunction
#
#
#
class LineBeam(NamedTuple):
    """
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
    #
    name: Union[int,str]
    load_type: str = "line"
    #
    def __str__(self) -> str:
        """ """
        output = ""
        output += "Line      qy [N/m]    qz [N/m]    L [m  ]     End\n"
        output += "{:8s} {: 1.4E} {: 1.4E} {: 1.4E}  [1]\n".format(str(self.name), self.qy1, self.qz1, self.L1)
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
class Load(Mapping):
    
    __slots__ = ['_loads', '_labels', 'cls']
    
    def __init__(self, cls):
        """
        """
        self.cls = cls
        self._loads:List[Tuple] = []
        self._labels:List[Union[int,str]] = []
    #
    #def __setitem__(self, load_name: Union[int, int], 
    #                parameters: Union[List[float], Dict[str, float]]) -> None:
    #    """
    #    farg = [name, connectivity, material, section, type, group]
    #    """
    #    print('--')
    #
    def __getitem__(self, load_name: Union[str,int]):
        """
        """
        return LoadType(load_name, cls=self)
    #
    def __len__(self) -> float:
        """ """
        return len(self._labels)

    def __iter__(self) -> Iterator:
        """ """
        for index, item in enumerate(self._labels):
            output = {"in_plane":SingFunction(0, 0), 
                      "out_plane":SingFunction(0, 0)}
            load = self._loads[index]
            if load.load_type == 'line':
                if load.qy1 or load.qy2:
                    output["in_plane"] = Trapezoidal(load.qy1, load.qy2,
                                                     self.cls.beam_length,
                                                     load.L1, load.L2)
                if load.qz1 or load.qz2:
                    output["out_plane"] = Trapezoidal(load.qz1, load.qz2,
                                                      self.cls.beam_length,
                                                      load.L1, load.L2)
                yield output                
            elif load.load_type == 'point':
                if load.fy:
                    output["in_plane"] = Point(load.fy, 
                                               self.cls.beam_length,
                                               load.L1)
                if load.fz:
                    output["out_plane"] = Point(load.fz, 
                                                self.cls.beam_length,
                                                load.L1)
                yield output
            elif load.load_type == 'moment':
                if load.mz:
                    output["in_plane"] = Moment(load.mz, 
                                                self.cls.beam_length,
                                                load.L1)
                if load.my:
                    output["out_plane"] = Moment(load.my, 
                                                 self.cls.beam_length,
                                                 load.L1)
                yield output
            else:
                raise IOError("wrong type")
        #return iter(self._labels)

    def __contains__(self, value) -> bool:
        return value in self._labels
    #
    def __str__(self) -> str:
        """ """
        output = "------- Loading\n"
        output += "\n"
        for x, load_name in enumerate(self._labels):
            index = self._labels.index(load_name)
            load = self._loads[index]            
            #output += "Support {:} : {:}\n".format(support_name, self._fixity[x])
            output += load.__str__()
            output += "\n"
        return output    
#
class LoadType:
    __slots__ = ['cls', '_load_name']
    
    def __init__(self, load_name, cls):
        """
        """
        self.cls = cls
        self._load_name = load_name
    #
    #
    #def __getitem__(self, load_name: Union[str,int]):
    #    """
    #    """
    #    print('--')
    #
    @property
    def point(self):
        """ """
        index = self.cls._labels.index(self._load_name)
        return self.cls._loads[index]
        #print('--')
    
    @point.setter
    def point(self, load:Union[List,Dict]):
        """ """
        load = self._get_point_load(load)
        load.extend([self._load_name, "point"])
        self.cls._loads.append(PointBeam._make(load))
        self.cls._labels.append(self._load_name)
    #
    #
    @property
    def moment(self):
        """ """
        index = self.cls._labels.index(self._load_name)
        return self.cls._loads[index]
    
    @moment.setter
    def moment(self, load:Union[List,Dict]):
        """ """
        load = self._get_moment_load(load)
        load.extend([self._load_name, "moment"])
        self.cls._loads.append(MomentBeam._make(load))        
        self.cls._labels.append(self._load_name)  
    #
    #
    @property
    def line(self):
        """ """
        index = self.cls._labels.index(self._load_name)
        return self.cls._loads[index]
    
    @line.setter
    def line(self, load:Union[List,Dict]):
        """ """
        load = self._get_line_load(load)
        load.extend([self._load_name, "line"])
        self.cls._loads.append(LineBeam._make(load))          
        self.cls._labels.append(self._load_name)
    #
    #
    @property
    def torsion(self):
        """ """
        index = self.cls._labels.index(self._load_name)
        load = self.cls._loads[index]
        print('--')
    #
    #
    def _get_line_load(self, load):
        """ """
        if isinstance(load, (list, tuple)):
            try:
                load = check_list(load)
            except AttributeError:
                return load
        elif isinstance(load, dict):
            load = check_beam_dic(load)
        else:
            raise Exception('   *** Load input format not recognized')
        return load
    #
    def _get_point_load(self, load):
        """ """
        if isinstance(load, (list, tuple)):
            load = check_list(load)
        elif isinstance(load, dict):
            load = check_point_dic(load)
        else:
            raise Exception('   *** Load input format not recognized')
        return load
    #
    def _get_moment_load(self, load):
        """ """
        if isinstance(load, (list, tuple)):
            load = check_list(load)
        elif isinstance(load, dict):
            load = check_moment_dic(load)
        else:
            raise Exception('   *** Load input format not recognized')
        #load = [-load[0], load[1], load[2]]
        return load
#
#
class Combination(Mapping):
    __slots__ = ['_combination','_load', '_labels', 'cls']
    
    def __init__(self, cls):
        """
        """
        self.cls = cls
        self._combination:List[float] = []
        self._load:List[Union[int,str]] = []  
        self._labels:List[Union[int,str]] = []    
    #
    def __setitem__(self, load_name: Union[int, int], 
                    parameters: Union[List[float], Dict[str, float]]) -> None:
        """
        farg = [name, connectivity, material, section, type, group]
        """
        #print('--')
        self._labels.append(load_name)
        self._load.append(parameters[0])
        self._combination.append(parameters[1])
    #
    def __getitem__(self, load_name: Union[str,int]):
        """
        """
        print('--')
    #
    def __len__(self) -> float:
        """ """
        return len(self._labels)

    def __iter__(self) -> Iterator:
        """ """
        return iter(self._labels)

    def __contains__(self, value) -> bool:
        return value in self._labels
    #    
#
#
def get_value(data, label:str, steps:int):
    """ """
    new_data = []
    for x in range(steps):
        try:
            new_data.append(data[label][x])
        except IndexError:
            new_data.append(0)
    return new_data
#
def get_load_unit(item):
    """ """
    if item.units() == 'metre*second^-2*gram': # N 
        return "newton"
    elif item.units() == 'metre^2*second^-2*gram': # N*m
        return "newton*metre"
    elif item.units() == 'second^-2*gram': # N/m
        return "newton/metre"
    elif item.units() == 'metre': # m
        return "metre"
    else:
        raise IOError("units {:} not compatible"
                      .format(item.units()))
#
def check_list(data)->List[float]:
    """ 
    line = [qy1,qz1, qy2,qz2, L1,L2]
    """
    load = defaultdict(list)
    L = []
    lendata = len(data)
    for x in range(lendata):
        unit = get_load_unit(data[x])
        #
        if unit == 'metre':
            L.append(data[x].value)
        elif unit == 'newton':
            load['N'].append(data[x].convert(unit).value)
        elif unit == 'newton*metre':
            load['N*m'].append(data[x].convert(unit).value)
        else: # "newton/metre"
            load['N/m'].append(data[x].convert(unit).value)
    #
    if len(load) > 1:
        raise IOError("units not compatible")
    elif 'N' in load:
        new_data = get_value(load, label='N', steps=2)
        new_data.append(L[0])
    elif 'N*m' in load:
        new_data = get_value(load, label='N*m', steps=2)
        new_data.append(L[0])
    else: # 'N/m'
        new_data = fill_udl_list(load['N/m'], L)
    #
    return new_data
#
def check_point_dic(data)->List[float]:
    """ """
    new_data = [0,0,0]
    for key, item in data.items():
        #if re.match(r"\b(fx|axial)\b", str(key), re.IGNORECASE):
        #    new_data[0] = item.value
        if re.match(r"\b(py|fy|in(_)?plane)\b", str(key), re.IGNORECASE):
            new_data[0] = item.convert("newton").value
        elif re.match(r"\b(pz|fz|out(_)?plane)\b", str(key), re.IGNORECASE):
            new_data[1] = item.convert("newton").value
        elif re.match ( r"\b((l|d(istance)?)(_)?(1|start))\b", str ( key ), re.IGNORECASE ):
            new_data[2] = item.value
    return new_data
#
def check_moment_dic(data)->List[float]:
    """ """
    new_data = [0,0,0]
    for key, item in data.items():
        #if re.match(r"\b(mx|torsion)\b", str(key), re.IGNORECASE):
        #    new_data[0] = item.value
        if re.match(r"\b(my|in(_)?plane)\b", str(key), re.IGNORECASE):
            new_data[0] = item.convert("newton*metre").value
        elif re.match(r"\b(mz|out(_)?plane)\b", str(key), re.IGNORECASE):
            new_data[1] = item.convert("newton*metre").value
        elif re.match(r"\b((l|d(istance)?)(_)?(1|start))\b", str(key), re.IGNORECASE):
            new_data[2] = item.value
    return new_data
#
def check_beam_dic(data)->List[float]:
    """ 
    [qy1, qz1, qy2, qz2, L1, L2]
    """
    new_data = [0,0, None, None,0,0]
    for key, item in data.items():
        if re.match(r"\b((qy|in(_)?plane)(_)?(1|start)?)\b", str(key), re.IGNORECASE):
            new_data[0] = item.convert("newton/metre").value
        elif re.match(r"\b((qy|in(_)?plane)(_)?(2|end))\b", str(key), re.IGNORECASE):
            new_data[2] = item.convert("newton/metre").value
        elif re.match(r"\b((qz|out(_)?plane)(_)?(1|start)?)\b", str(key), re.IGNORECASE):
            new_data[1] = item.convert("newton/metre").value        
        elif re.match(r"\b((qz|out(_)?plane)(_)?(2|end))\b", str(key), re.IGNORECASE):
            new_data[3] = item.convert("newton/metre").value
        elif re.match(r"\b((l|d(_)?)(1|start))\b", key, re.IGNORECASE):
            new_data[4] = item.value
        elif re.match(r"\b((l|d(_)?)(2|end))\b", key, re.IGNORECASE):
            new_data[5] = item.value
    #
    if new_data[2] == None:
        new_data[2] = new_data[0]
    
    if new_data[3] == None:
        new_data[3] = new_data[1]
    #
    return new_data
#
def fill_udl_list(loading, L):
    """ """
    new_data = [0,0, None, None]
    for x in range(4):
        try:
            new_data[x] = loading[x]
        except IndexError:
            break
    #
    if new_data[2] == None:
        new_data[2] = new_data[0]
    
    if new_data[3] == None:
        new_data[3] = new_data[1]
    # update 
    Lnew = [0, 0]
    for x in range(2):
        try:
            Lnew[x] = L[x]
        except IndexError:
            break
    #
    new_data.extend(Lnew)
    return new_data
#