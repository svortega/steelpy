#
# Copyright (c) 2019-2021 steelpy
#

# Python stdlib imports
from array import array
#from collections.abc import Mapping
from typing import NamedTuple, Dict, List, Tuple, Union
from dataclasses import dataclass
import re

# package imports
from steelpy.beam.static.support import Response, ReactPlane, FuncRes
from steelpy.f2uModel.plot.simplot import SimPlot
from steelpy.process.math.vector import Vector
#from steelpy.f2uModel.load.operations.actions import Actions
from steelpy.f2uModel.sections.process.stress import BeamStress
from steelpy.process.units.main import Units
#
#
#
class BeamPlane(NamedTuple):
    inplane:List[float]
    outplane:List[float]
    x:List[float]
    load_type:str
    #
    def _get_data(self, plane:str):
        """ """
        # if re.match(r"\b(my|in(_)?plane)\b", str(key), re.IGNORECASE):
        if re.match(r"\b(z|out(_)?plane)\b", plane, re.IGNORECASE):
            y = [item/1 for item in self.outplane]
        else:
            y = [item/1 for item in self.inplane]
        return y
    #
    def plot(self, plane:Union[str,None]=None):
        """ """
        #plane
        y = self._get_data(plane)
        data = list(zip(self.x, y))
        yhigh = max(y) * 1.20
        ylow  = min(y) * 1.20
        plt = SimPlot()
        y_axis = (ylow, yhigh)
        if sum(y_axis) == 0:
           y_axis = "automatic" 
        #try:
        test = plt.plotLine(data, title=plane.upper(), yaxis=y_axis,
                            xlab="Beam Length [m]", ylab=self.load_type.capitalize(), 
                            smooth=0, width=2, color="blue")
        #except ValueError:
        #    test = plt.plotLine(data, title=plane.upper(), #yaxis=(ylow, yhigh),
        #                        xlab="Beam Length [m]", ylab=self.load_type.capitalize(), 
        #                        smooth=0, width=2, color="blue")            
        #
        #plt.mainloop()
        test.mainloop()
    #
    def __call__(self, plane:Union[str,None]=None, 
                 units=None):
        """ """
        x = Vector(self.x)
        if plane:
            y = Vector(self._get_data(plane))
            return x, y
        else:
            y = Vector(self._get_data("in_plane"))
            z = Vector(self._get_data("out_plane"))
            return x, y, z
#
#
@dataclass
class BeamResponse:
    __slots__ = ['_loads', '_labels', '_result', 'cls', '_response']

    def __init__(self, cls):
        """
        """
        self.cls = cls
        self._response = ResponseResult()
    #
    def _get_steps(self, steps:int):
        """ """
        L = self.cls.beam_length
        step = L/steps
        x_steps = [x*step for x in range(steps+1)]
        load_step = []
        for loads in self.cls._load:
            for plane in ["in_plane", "out_plane"]:
                load = loads[plane]
                load_step.append(load.L1)
                # udl
                try:
                    res = load.L - load.L2
                    load_step.append(res)
                except AttributeError:
                    try:
                        1/load.L1 
                        load_step.append(load.L1 - 0.010)
                        load_step.append(load.L1 + 0.010)
                        #load_step.extend([load.L1 - x*0.010 for x in range(10)])
                        #load_step.extend([load.L1 + x*0.010 for x in range(10)])
                    except ZeroDivisionError:
                        continue
        #
        x_steps.extend(load_step)
        x_steps = sorted(list(set(x_steps)))
        x_steps = [item for item in x_steps if item <= L]
        return x_steps
    #
    def _get_properties(self):
        """Get material and section data"""
        sect = self.cls.section._get_properties()
        I_plane = {"in_plane":sect.Iy, "out_plane":sect.Iz}
        #
        mat = self.cls.material
        E = mat.E.convert("pascal").value
        L = self.cls.beam_length
        return L, I_plane, E    
    #
    def _get_load(self, steps:int):
        """
        steps : inst = beam length steps
        plane : str = in_plane/out_plane
        """
        L, I_plane, E = self._get_properties()
        loadres = {"in_plane":None, "out_plane":None}
        x_steps = self._get_steps(steps)
        no_steps = len(x_steps)
        for plane in loadres.keys():
            V=[0]*no_steps; M=[0]*no_steps 
            theta=[0]*no_steps; w=[0]*no_steps
            for loads in self.cls._load:
                load = loads[plane]
                for x, step in enumerate(x_steps):
                    V[x] += load.V(step)
                    M[x] += load.M(step)
                    theta[x] += load.theta(step, E, I_plane[plane])
                    w[x] += load.w(step, E, I_plane[plane])
            loadres[plane] = FuncRes(V, M, theta, w, x_steps)
        return ReactPlane(loadres["in_plane"], loadres["out_plane"])
    #
    def _solve_response(self, steps:int, E:float, I:float, 
                        load, reaction):
        """ """
        result = Response(E=E, I=I)
        # TODO : fix negative sign
        result.reacctions(V0=-reaction.R.convert("newton").value,
                          M0=-reaction.M.convert("newton*metre").value,
                          theta0=reaction.theta.value,
                          w0=reaction.w.value)
        #
        #result = reaction.response
        #
        V=[]; M=[]; theta=[]; w=[]
        x_steps = load.x
        for index, step in enumerate(x_steps):
            result.load(FV=load.V[index], FM=load.M[index],
                        Ftheta=load.theta[index], Fw=load.w[index])
            V.append(result.V(step))
            M.append(result.M(step))
            theta.append(result.theta(step))
            w.append(result.w(step))
        #print('-->')
        return FuncRes(V, M, theta, w, x_steps)
    #
    def __call__(self, steps:int=10):
        """ """
        try:
            self._result
        except AttributeError:
            L, I_plane, E = self._get_properties()
            load = self._get_load(steps=steps)
            #
            reaction = self.cls.support[1].in_plane
            in_plane = self._solve_response(steps=steps, 
                                            E= E,
                                            I=I_plane['in_plane'],
                                            load=load.in_plane,
                                            reaction=reaction)
            #
            reaction = self.cls.support[1].out_plane
            out_plane = self._solve_response(steps=steps,
                                             E= E,
                                             I=I_plane['out_plane'],                                             
                                             load=load.out_plane,
                                             reaction=reaction)
            # output
            #self._result = ReactPlane(in_plane, out_plane)
            #
            name = [self.cls.beam_name] * len(in_plane.x)
            self._response.in_plane = [name, *in_plane]
            self._response.out_plane = [name, *out_plane]
            #print('<-- mat' )
            self._result = self._response[self.cls.beam_name]
        return self._result
    #
    def shear(self, steps:int=10)-> Tuple:
        """ """
        _result = self.__call__(steps=steps)
        return BeamPlane(_result.Vy, _result.Vz,
                         _result.x,"shear")
    #
    def bending(self, steps:int=10)-> Tuple:
        """ """
        _result = self.__call__(steps=steps)
        return BeamPlane(_result.My, _result.Mz,
                         _result.x,"bending")
    #
    def slope(self, steps:int=10)-> Tuple:
        """ """
        _result = self.__call__(steps=steps)
        return BeamPlane(_result.theta_y, _result.theta_z,
                         _result.x,"slope")
    #
    def deflection(self, steps:int=10)-> Tuple:
        """ """
        _result = self.__call__(steps=steps)
        return BeamPlane(_result.wy, _result.wz,
                         _result.x,"deflection")
    #
    def stress(self, steps:int=10):
        """ """
        _result = self.__call__(steps=steps)
        stress = []
        for x in self._result.x:
            #res = self.cls.section.stress(self._result.forces(x))
            stress.append(self.cls.section.stress(_result.forces(x)))
        #stress = list(zip(self._result.x, stress))
        #print('-->')
        return stress
#
#
#
class Actions(NamedTuple):
    """
    Force & bending moments
    """
    x:Union[float, Units]
    Fx:Union[float, Units]
    Fy:Union[float, Units]
    Fz:Union[float, Units]
    Mx:Union[float, Units]
    My:Union[float, Units]
    Mz:Union[float, Units]
#
#
class ReSummary(NamedTuple):
    x:List[float]
    Px:List[float]
    Vy:List[float]
    Vz:List[float]
    Mx:List[float]
    My:List[float]
    Mz:List[float]
    theta_x:List[float]
    theta_y:List[float]
    theta_z:List[float]
    wx:List[float]
    wy:List[float]
    wz:List[float]
    #
    def forces(self, x):
        """ """
        #print('-->')
        index = self.x.index(x)
        return Actions(x, self.Px[index], self.Vy[index], self.Vz[index],
                       self.Mx[index], self.My[index], self.Mz[index])
#
#@dataclass
class ResponseResult:
    """
    Response
    x     : distance from end 1
    P     : Axial load
    V     : Shear
    M     : Moment
    T     : Torque
    phi   : Angle of twist
    theta : Slope
    w     : Deflection
    """
    __slots__ = ['_x', '_labels',
                 '_Px', '_Vy', '_Vz',
                 '_Mx', '_My', '_Mz', 
                 '_theta_x', '_theta_y', '_theta_z',
                 '_wx', '_wy', '_wz']
    
    def __init__(self):
        """
        """
        self._labels:List = []
        self._x:List = []
        self._Px = array('f', [])
        self._Vy = array('f', [])
        self._Vz = array('f', [])
        self._Mx = array('f', [])
        self._My = array('f', [])
        self._Mz = array('f', [])
        #
        self._theta_x = array('f', [])
        self._theta_y = array('f', [])
        self._theta_z = array('f', [])        
        #
        self._wx = array('f', [])
        self._wy = array('f', [])
        self._wz = array('f', [])
    #
    def __setitem__(self, name: int,
                    response: Union[List[float], Dict[str, float]]) -> None:
        """
        [x, Px, Vy, Vz, Mx, My, Mz, theta_x, theta_y, theta_z, wx, wy, wz]
        """
        #try:
        #    self._labels.index(name)
        #    raise Exception('    *** warning name {:} already exist'
        #                    .format(name))
        #except ValueError:
        self._labels.append(name)
        self._x.append(response[0])
        self._Px.append(response[1])
        self._Vy.append(response[2])
        self._Vz.append(response[3])
        self._Mx.append(response[4])
        self._My.append(response[5])
        self._Mz.append(response[6])
        #
        self._theta_x.append(response[7])
        self._theta_y.append(response[8])
        self._theta_z.append(response[9])
        #
        self._wx.append(response[10])
        self._wy.append(response[11])
        self._wz.append(response[12])
    #
    def __getitem__(self, name:int):
        """
        """
        d = {x:self._x[x] for x, item in enumerate(self._labels)
             if item == name}
        #items = {k: d[k] for k in sorted(d, key=d.get)}
        #d = [x for x, item in enumerate(self._labels)
        #     if item == name]
        res = []
        for index in sorted(d, key=d.get):
        #for index in d:
            #print(self._x[index])
            res.append([self._x[index], self._Px[index], self._Vy[index], self._Vz[index],
                        self._Mx[index], self._My[index], self._Mz[index],
                        self._theta_x[index], self._theta_y[index], self._theta_z[index],
                        self._wx[index], self._wy[index], self._wz[index]])
        #
        res = list(zip(*res))
        #print('<-- mat' )
        return ReSummary(*res)
    #
    @property
    def in_plane(self):
        """ """
        res =  FuncRes(V, M, theta, w, x_steps)
        print('<-- mat' )
    
    @in_plane.setter
    def in_plane(self, response:List[List[float]]):
        """ """
        #for x, name in enumerate(response[0]):
        try:
            self._labels.index(response[0][0])
            for x, name in enumerate(response[0]):
                index = self._x.index(response[5][x])
                self._Vy[index] = response[1][x]
                self._My[index] = response[2][x]
                self._theta_y[index] = response[3][x]
                self._wy[index] = response[4][x]
        except ValueError:
            for x, name in enumerate(response[0]):
                self.__setitem__(name, [response[5][x], 
                                        0, response[1][x], 0,
                                        0, response[2][x], 0,
                                        0, response[3][x], 0,
                                        0, response[4][x], 0])
        
    #
    @property
    def out_plane(self):
        """ """
        print('<-- mat' )
    
    @out_plane.setter
    def out_plane(self, response:List[float]):
        """ """
        try:
            self._labels.index(response[0][0])
            for x, name in enumerate(response[0]):
                index = self._x.index(response[5][x])
                self._Vz[index] = response[1][x]
                self._Mz[index] = response[2][x]
                self._theta_z[index] = response[3][x]
                self._wz[index] = response[4][x]
            #print('<-- mat' )
        except ValueError:
            for x, name in enumerate(response[0]):
                self.__setitem__(name, [response[5][x], 
                                        0, 0, response[1][x],
                                        0, 0, response[2][x],
                                        0, 0, response[3][x],
                                        0, 0, response[4][x]])
    #      
#
#
#
#