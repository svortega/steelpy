#
# Copyright (c) 2019 steelpy
#
# Python stdlib imports
from __future__ import annotations
from bisect import bisect_right
from dataclasses import dataclass
from typing import NamedTuple
#
# package imports
#
#from steelpy.utils.math.vector import Vector

# ---------------------------------------------------------------
#
# bload = BeamLoad(L, E) #, I)
# beload['title'] = ['point', P, L1]
# beload['title'] = ['moment', M, L1]
# beload['title'] = ['torsion', T, L1]
# beload['title'] = ['linear', q1, L1, q2, L2]
#
# bload.point('title', P, L1)
# bload.moment('title', M, L1)
# bload.torsion('title', T, L1)
# bload.linear('title', q1, L1, q2, L2)
# load at x distance from end 0 (left support)
# result = bload(x, I, load_name)
#
#
class BeamLoad:
    __slots__ = ['L', 'E', 'I', 'L1', '_load', 
                 '_type', '_labels']

    def __init__(self, L: float, E: float = 2.05e11) -> None:
        """
        L : beam lenght [m]
        E : Elastic modulus [Pa] (default steel)
        """
        self.L: float = L
        self.E:float = E
        #
        self._load = []
        #
        self._labels:list[str|int] = []
        self._type:list[str] = []
    #
    #
    def __setitem__(self, load_name:str|int,
                    load:list[str|float]) -> None:
        """
        """
        load_type = load[0]

        if re.match(r"\b(point)\b", load_type, re.IGNORECASE):
            self.point(load_name=load_name, P=load[1], L1=load[2])

        elif re.match(r"\b(moment)\b", load_type, re.IGNORECASE):
            self.moment(load_name=load_name, m=load[1], L1=load[2])

        elif re.match(r"\b(udl|linear)\b", load_type, re.IGNORECASE):
            self.linear(load_name, *load[1:])

        elif re.match(r"\b(torsion)\b", load_type, re.IGNORECASE):
            self.torsion(load_name=load_name, To=load[1], L1=load[2])

        else:
            raise IOError(f"load type {load_type} no supported")

    
    def __getitem__(self, load_name:str|int):
        """
        """
        index = [i for i, item in enumerate(self._labels)
                 if item == load_name]
        #
        #load = []
        #for idx in index:
        #    load += self._load[idx]
        return [self._load[idx] for idx in index]
    #
    #
    def point(self, load_name:str|int, P:float, L1:float):
        """
               |
             P |
        o------V------------------o
        |                          |
        +  L1  +                   +
        """
        self._labels.append(load_name)
        self._load.append(Point(P=P, L=self.L, L1=L1))
        self._type.append("point")

    def moment(self, load_name:str|int, M: float,  L1: float):
        """
               
               M 
        o------@------------------o
        |                         |
        +  L1  +                  +
        """
        self._labels.append(load_name)
        self._load.append(Moment(m=M, L=self.L, L1=L1))
        self._type.append("moment")

    def linear(self, load_name:str|int, q1: float, L1: float, 
               q2: float|None=None, L2: float|None = None):
        """
                        |
             q1         | q2
        o------|        |----------o
        |                          |
        +  L1  +        +    L2    +
        """
        if not q2:
            q2 = q1

        if not L2:
            L2 = L1
        
        self._labels.append(load_name)
        self._load.append(Trapezoidal(q1=q1, q2=q2,
                                      L=self.L, L1=L1, L2=L2))
        self._type.append("linear")

    def torsion(self, load_name:str|int, T:float, L1:float):
        """
               T 
        o------|------------------o
        |                         |
        +  L1  +                  +
        """
        self._labels.append(load_name)
        self._load.append(TorsionPoint(To=T, L=self.L, L1=L1))
        self._type.append("torsion")
    #
    def __call__(self, x: float, I: float, 
                 load_list:str|int|list|None=None):
        """
        x : distance from end 1 [m]
        I : Moment of intertia [m^4]
        load_list : list of load to be included in calculation (default use all)

        return: 
        [V, M, w, theta]
        """
        self.I = I
        # [V, M, w, theta]
        load = Vector([0,0,0,0])
        if not load_list:
            for item in self._load:
                load += item(x=x, E=self.E, I=self.I)
        else:
            for load_name in load_list:
                new_load = self.__getitem__(load_name)
                for item in new_load:
                    load += item(x=x, E=self.E, I=self.I)
            #1/0
        return load
#
#
@dataclass
class SingularFunction:

    def function_n(self, step: float, n: int) -> float:
        """<x-L>^n"""
        if n < 0:
            return 0
        elif step < 0:
            return 0
        elif n == 0:
            return 1
        else:
            return step**n
    #
    def fun_sincos(self, fun, step: float, C: float) -> float:
        """<x-L>^n"""
        if step < 0:
            return 0
        else:
            return fun(C * step)
    #
    def fun_evar(self, step: float, n: int):
        """
        """
        return self.function_n(step, n)

#
#
#
# ---------------------------------------------------------------
#
#
class Interpolate:
    def __init__(self, x_list, y_list):
        if any(y - x <= 0 for x, y in zip(x_list, x_list[1:])):
            raise ValueError("x_list must be in strictly ascending order!")
        self.x_list = x_list
        self.y_list = y_list
        intervals = zip(x_list, x_list[1:], y_list, y_list[1:])
        self.slopes = [(y2 - y1) / (x2 - x1) for x1, x2, y1, y2 in intervals]

    def __call__(self, x):
        if not (self.x_list[0] <= x <= self.x_list[-1]):
            raise ValueError("x out of bounds!")
        if x == self.x_list[-1]:
            return self.y_list[-1]
        i = bisect_right(self.x_list, x) - 1
        return self.y_list[i] + self.slopes[i] * (x - self.x_list[i])
#