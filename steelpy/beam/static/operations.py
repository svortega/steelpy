#
# Copyright (c) 2019-2020 steelpy
#

# Python stdlib imports
#from collections.abc import Mapping
from typing import NamedTuple, Dict, List, Tuple, Union
from dataclasses import dataclass

# package imports
from steelpy.beam.static.support import Response, FuncRes, ReactPlane
from steelpy.f2uModel.plot.simplot import SimPlot


#
#
#
class BeamPlane(NamedTuple):
    inplane:List[float]
    outplane:List[float]
    x:List[float]
    load_type:str
    def plot(self,  plane:Union[str,None]=None):
        """ """
        #plane
        y = [item/1 for item in self.inplane]
        yhigh = max(y) * 1.20
        ylow  = min(y) * 1.20
        data = list(zip(self.x, y))
        plt = SimPlot()
        test = plt.plotLine(data, title='test curve', #yaxis=(ylow, yhigh),
                     xlab="Beam Length ", ylab=self.load_type, smooth=0)
        #plt.mainloop()
        test.mainloop()
#
#
@dataclass
class BeamOps:
    __slots__ = [ '_loads', '_labels', 'cls', '_result']

    def __init__(self, cls):
        """
        """
        self.cls = cls
    #
    def _get_load(self, steps:int):
        """
        steps : inst = beam length steps
        plane : str = in_plane/out_plane
        """
        L = self.cls._beam_length
        sect = self.cls.section._get_properties()
        I_plane = {"in_plane":sect.Iy, "out_plane":sect.Iy}
        mat = self.cls.material
        E = mat.E.convert ( "pascal" ).value
        loadres = {"in_plane":None, "out_plane":None}
        step = L/steps
        x_steps = [x*step for x in range(steps+1)]
        for plane in loadres.keys():
            V=[]; M=[]; theta=[]; w=[]
            for load in self.cls._load:
                for x in x_steps:
                    V.append(load[plane].V(x))
                    M.append(load[plane].M(x))
                    theta.append(load[plane].theta(x, E, I_plane[plane]))
                    w.append(load[plane].w(x, E, I_plane[plane]))
                loadres[plane] = FuncRes(V, M, theta, w, x_steps)
        return ReactPlane(loadres["in_plane"], loadres["out_plane"])
    #
    def _solve_response(self, steps:int, E:float, I:float, 
                        load, reaction):
        """ """
        result = Response(E=E, I=I)
        result.reacctions(V0=reaction.R.convert("newton").value,
                          M0=reaction.M.convert("newton*metre").value,
                          theta0=reaction.theta.value,
                          w0=reaction.w.value)
        #
        V=[]; M=[]; theta=[]; w=[]
        x_steps = load.x
        for index, x in enumerate(x_steps):
            result.load(FV=load.V[index], FM=load.M[index],
                             Ftheta=load.theta[index], Fw=load.w[index])
            V.append(result.V(x))
            M.append(result.M(x))
            theta.append(result.theta(x))
            w.append(result.w(x))
        #print('-->')
        return FuncRes( V, M, theta, w, x_steps )
    #
    def solve(self, steps:int=10):
        """ """
        try:
            self._result
        except AttributeError:
            sect = self.cls.section._get_properties()
            mat = self.cls.material
            E = mat.E.convert("pascal").value
            #
            load = self._get_load(steps=steps)
            #
            reaction = self.cls._reactions[0].in_plane
            in_plane = self._solve_response(steps=steps, 
                                            E= E,
                                            I=sect.Iy,
                                            load=load.in_plane,
                                            reaction=reaction)
            #
            reaction = self.cls._reactions[0].out_plane
            out_plane = self._solve_response(steps=steps,
                                             E= E,
                                             I=sect.Iz,                                             
                                             load=load.out_plane,
                                             reaction=reaction)
            # output
            self._result = ReactPlane(in_plane, out_plane)
    #
    def shear(self, steps:int=10)-> Tuple:
        """ """
        self.solve(steps=steps)
        return BeamPlane(self._result.in_plane.V,
                         self._result.out_plane.V,
                         self._result.in_plane.x, "shear")
    #
    def bending(self, steps:int=10)-> Tuple:
        """ """
        self.solve(steps=steps)
        return BeamPlane(self._result.in_plane.M,
                         self._result.out_plane.M,
                         self._result.in_plane.x, "bending")
    #
    def slope(self, steps:int=10)-> Tuple:
        """ """
        self.solve(steps=steps)
        return BeamPlane(self._result.in_plane.theta,
                         self._result.out_plane.theta,
                         self._result.in_plane.x, "slope")
    #
    def deflection(self, steps:int=10)-> Tuple:
        """ """
        self.solve(steps=steps)
        return BeamPlane(self._result.in_plane.w, 
                         self._result.out_plane.w,
                         self._result.in_plane.x, "deflection" )
#
#
#
#
#
#