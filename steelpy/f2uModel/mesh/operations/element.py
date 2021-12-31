#
# Copyright (c) 2009-2021 fem2ufo
#
# Python stdlib imports
#from itertools import chain
#from collections import Counter
#from collections import defaultdict
from dataclasses import dataclass
from math import factorial, cosh, sinh
#import functools
#import re
from typing import NamedTuple, Dict, Union, Tuple, List
#from math import  isclose, dist
#
# package imports
from steelpy.f2uModel.mesh.operations.load import BeamLoad

@dataclass
class BeamBasic:
    __slots__ = ['_support', '_load',  '_response',
                 'L', 'E', 'I', 'k', '_end_force']

    #
    def __init__(self, L: float,  E: float, Iy: float, Iz: float) -> None:
        """
        L : beam lengh
        """
        self.L:float = L
        self.E:float = E
        self.I:List[float] = [Iy, Iz]
        self.k:List[float] = [0, 0]
        self._load = BeamLoad(L=self.L)
        #self._response = Response( E=E, I=self.I )
    #
    def spring(self, k1: float, k2: float):
        """Get material and section data"""
        self.k = [k1, k2]

    #
    def supports(self, fixity):
        """ solve support 1 according determined boundary conditions"""
        supp_1 = fixity[ 0 ]
        try:
            supp_2 = fixity[ 1 ]
        except IndexError:
            supp_2 = [0,0,0,0,0,0]
        #
        self._support = []
        for x in range(2):
            boun1 = get_fixity(supp_1[0+x], supp_1[1+x], supp_1[4+x])
            boun2 = get_fixity(supp_2[0+x], supp_2[1+x], supp_2[4+x])
            function = self._support_func(boun1, boun2)
            self._support.append(function(L=self.L, E=self.E, I=self.I[x],
                                          k1=self.k[0], k2=self.k[1]))
            #self._support[-1](L=self.L, E=self.E, I=self.I[x],
            #                  k1=self.k[0], k2=self.k[1])
    #
    @property
    def load(self):
        """ beam load"""
        return self._load
    #
    @property
    def reactions(self):
        """ """
        load_total = self.load.response(self.L, self.E,
                                        self.I[0], self.I[1])
        res = []
        for x, support in enumerate(self._support):
            support.load(*load_total[x])
            res.append(support.reactions)
        return res
    #
    @property
    def end_force(self):
        """ []"""
        return self._end_force

    @end_force.setter
    def end_force(self, value):
        """ """
        self._end_force = value
    #
    def response(self, R1):
        """ 
        R1 = [FV,FM,Fw,Ftheta]
        """
        x_steps = self.load.get_steps()
        step_load = []
        for step in x_steps:
            step_load.append(self.load.response(step*self.L, self.E,
                                                self.I[0], self.I[1]))        
        #
        FuncRes = []
        res = []
        for x in range(2):
            res.append(Response(E=self.E, I=self.I[x]))
            res[-1].reacctions(*R1[x])
            #
            V=[]; M=[]; theta=[]; w=[]
            for index, xstep in enumerate(x_steps):
                step = xstep*self.L
                res[-1].load(*step_load[index][x])
                V.append(res[-1].V(step))
                M.append(res[-1].M(step))
                w.append(res[-1].w(step))
                theta.append(res[-1].theta(step))
            #print('-->')
            FuncRes.append([V, M, w, theta, x_steps])
        return FuncRes
    #
    def _support_func(self, supp1, supp2):
        """ """
        if supp1 == "pinned":
            if supp2 == "pinned":
                return PinnedPinned
            elif supp2 == "fixed":
                return PinnedFixed
            elif supp2 == "guided":
                return PinnedGuided
            elif supp2 == "spring":
                return PinnedSpring
            else:
                raise IOError ( "unstable" )

        elif supp1 == "fixed":
            if supp2 == "pinned":
                return FixedPinned
            elif supp2 == "fixed":
                return FixedFixed
            elif supp2 == "guided":
                return FixedGuided
            elif supp2 == "free":
                return FixedFree
            elif supp2 == "spring":
                return FixedSpring
            else:
                raise IOError ( "boundary no supported" )

        elif supp1 == "free":
            if supp2 == "fixed":
                return FreeFixed
            elif supp2 == "spring":
                return FreeSpring
            else:
                raise IOError ( "unstable" )

        elif supp1 == "guided":
            if supp2 == "pinned":
                return GuidedPinned
            elif supp2 == "fixed":
                return GuidedFixed
            elif supp2 == "spring":
                return GuidedSpring
            else:
                raise IOError ( "unstable" )

        elif supp1 == "spring":
            if supp2 == "pinned":
                return SpringPinned
            elif supp2 == "fixed":
                return SpringFixed
            elif supp2 == "guided":
                return SpringGuided
            elif supp2 == "free":
                return SpringFree
            elif supp2 == "spring":
                return SpringSpring
            else:
                raise IOError ( "boundary no supported" )
        else:
            raise IOError ( "boundary no supported" )

    #
    def __call__(self, x):
        """ """
        res = []
        for support in self._support:
            res.append(support.__call__( x ))
        return res
#
#
def get_fixity(axial, fix_1, fix_2):
    """ """
    try:
        1/fix_1
        try:
            1/fix_2
            return "fixed"
        except ZeroDivisionError:
            try:
                1/axial
                return "pinned"
            except ZeroDivisionError:
                return "guided"
    except ZeroDivisionError:
        return "free"
#
#
class ReacResult(NamedTuple):
    """Reactions"""
    R:float
    M:float
    theta:float
    w:float

    def __str__(self) -> str:
        """ """
        output = ""
        output += "React  [N  ]  {: 1.4E}\n".format(self.R)
        output += "Moment [N*m]  {: 1.4E}\n".format(self.M)
        output += "theta  [rad]  {: 1.4E}\n".format(self.theta)
        output += "delta  [m  ]  {: 1.4E}\n".format(self.w)          
        return output
#
#
@dataclass
class Response:
    """Table 11.2, Part A
    Simple Beams with arbitrary loadings
    General expresions.
    Pilkey"""

    __slots__ = [ 'E', 'I', 'FV', 'FM', 'Ftheta', 'Fw',
                  'V0', 'M0', 'theta0', 'w0' ]

    def __init__(self, E: float, I: float) -> None:
        """
        """
        self.E = E
        self.I = I

    #
    def load(self, FV: float, FM: float,
             Fw: float, Ftheta: float) -> None:
        """ """
        self.FV = FV
        self.FM = FM
        self.Ftheta = Ftheta
        self.Fw = Fw

    #
    def reacctions(self, V0: float, M0: float,
                   w0: float, theta0: float) -> None:
        """ """
        self.V0 = V0
        self.M0 = M0
        self.w0 = w0
        self.theta0 = theta0
    
    #
    def V(self, x: float) -> float:
        """ Shear force"""
        return self.V0 + self.FV

    #
    def M(self, x: float) -> float:
        """ Bending moment"""
        return self.M0 + self.V0 * x + self.FM
    #
    def T(self, x: float) -> float:
        """Torsion """
        step = x - self.L1
        return  self.T0  + self.To * self.function_n(step, 0)
    #
    def theta(self, x: float) -> float:
        """ Slope"""
        return (self.theta0 + self.V0 * x ** 2 / (2 * self.E * self.I)
                + self.M0 * x / (self.E * self.I) 
                + self.Ftheta)
    #
    def phi(self, x: float) -> float:
        """ angle of rotation at a distance x from the left end (radians) """
        # F1, F2, F3, F4 = self.F(x)
        return self.phi0 + self.T0 * x / (self.G * self.J) + self.Fphi
    #
    def w(self, x: float) -> float:
        """ Deflection"""
        return (self.w0 - self.theta0 * x
                - self.V0 * x ** 3 / (factorial ( 3 ) * self.E * self.I)
                - self.M0 * x ** 2 / (2 * self.E * self.I) 
                + self.Fw)
    #
    #
    def function_n(self, step:float,  n:int) -> float:
        """ <x-L>^n """
        if n < 0:
            return 0
        elif step < 0:
            return 0
        elif n == 0:
            return 1
        else:
            return step**n
#
#
@dataclass
class SuppBasic:
    __slots__ = [ 'L', 'E', 'I', '_response',
                  'FV', 'FM', 'Ftheta', 'Fw',
                  'V1', 'M1', 'w1', 'theta1' ]

    def __init__(self, L: float, E: float, I: float) -> None:
        """
        """
        self.L: float = L
        self.E: float = E
        self.I: float = I
        #
        # self.FV:float = 0
        # self.FM:float = 0
        # self.Ftheta:float = 0
        # self.Fw:float = 0
        # spring stiffness
        # self.k1 = k1
        # self.k2 = k2
        self._response = Response( E=self.E, I=self.I )

    #
    def load(self, FV: float, FM: float, Fw: float, Ftheta: float):
        """ """
        self.FV: float = FV
        self.FM: float = FM
        self.Fw: float = Fw
        self.Ftheta: float = Ftheta
        # set response function
        self._response.load( FV=FV, FM=FM, Fw=Fw, Ftheta=Ftheta)
    #
    #def torsion(self):
    #    """ """
    #    pass
    #
    @property
    def V0(self) -> float:
        """ Shear """
        return 0
        #

    @property
    def M0(self) -> float:
        """ Bending moment"""
        return 0
    #
    @property
    def T0(self) -> float:
        """ Torsional moment"""
        return 0
    #
    @property
    def phi0(self):
        """ Angle of rotation"""
        return 0
    #
    @property
    def theta0(self):
        """ Slope"""
        return 0

    #
    @property
    def w0(self) -> float:
        """ Deflection"""
        return 0

    #
    #
    def _spring_1(self, k1: float):
        """ """
        #
        A7 = 1 + k1 * self.L / (3 * self.E * self.I)
        A8 = 1 / self.L + k1 / (2 * self.E * self.I)
        A9 = 1 + k1 * self.L / (4 * self.E * self.I)
        A10 = 1 + k1 * self.L / (self.E * self.I)
        A11 = 3 * self.E * self.I / self.L ** 3 + 3 * k1 / self.L ** 2
        return A7, A8, A9, A10, A11

    #
    def _spring_2(self, k2: float):
        """ """
        A1 = 1 / self.L - k2 / (2 * self.E * self.I)
        A2 = k2 * self.Ftheta - self.FM
        A3 = 1 - k2 * self.L / (3 * self.E * self.I)
        A4 = 1 - k2 * self.L / (4 * self.E * self.I)
        A5 = 1 - k2 * self.L / (2 * self.E * self.I)
        A6 = 1 - k2 * self.L / (self.E * self.I)
        return A1, A2, A3, A4, A5, A6

    #
    #
    @property
    def respose(self):
        """ """
        return self._response

    #
    @property
    def reactions(self):
        """
        L : Length of the beam
        """
        self.V1, self.M1, self.w1, self.theta1 = self.__call__( self.L )
        return [ [ self.V0, self.M0, self.w0, self.theta0],
                 [ self.V1, self.M1, self.w1, self.theta1] ]

    #
    def __call__(self, x: float) -> List[ float ]:
        """
        x : distance from end 1
        """
        # update load
        self._response.load(FV=self.FV, FM=self.FM,
                            Fw=self.Fw, Ftheta=self.Ftheta)
        self._response.reacctions(V0=self.V0, M0=self.M0,
                                  w0=self.w0, theta0=self.theta0)
        V = self._response.V( x )
        M = self._response.M( x )
        theta = self._response.theta( x )
        w = self._response.w( x )
        return [ V, M, w, theta]
    #
    # Torsion
    def F(self, x: float):
        """ x : distance from the left end"""
        #a = self.a.value
        beta_x = self.beta * x
        F1 = cosh(beta_x)
        F2 = sinh(beta_x)
        F3 = cosh(beta_x) - 1.0
        F4 = sinh(beta_x) - beta_x
        return [F1, F2, F3, F4]
#
#
# Pinned
#
@dataclass
class PinnedPinned ( SuppBasic ):
    """ """
    __slots__ = [ 'E', 'I', 'FV', 'FM', 'Ftheta', 'Fw' ]

    def __init__(self, L: float, E: float, I: float,
                 k1: float = 0, k2: float = 0) -> float:
        """
        """
        super().__init__ ( L, E, I )

    #
    @property
    def V0(self) -> float:
        """ Shear force"""
        return -self.FM / self.L

    #
    @property
    def theta0(self) -> float:
        """ Slope"""
        return self.Fw / self.L + self.FM * self.L / (6 * self.E * self.I)
    #
    # Torsion
    @property
    def phi0(self) -> float:
        """ angle of rotation (radians) """
        phiA_0 = 0
        phiA_I = (self.To/(self.Cw*self.E*self.beta**2)
                  * (1 - self.L1/self.L - Ca2/C2))
        phiA_II = 0
        phiA_III = 0
        return [phiA_0, phiA_I, phiA_II, phiA_III]
    #
    @property
    def T0(self) -> float:
        """ Torsional moment"""
        return -self.To * (1 - self.L1/self.L)
#
@dataclass
class PinnedFixed ( SuppBasic ):
    """ """
    __slots__ = [ 'E', 'I', 'FV', 'FM', 'Ftheta', 'Fw' ]

    def __init__(self, L: float, E: float, I: float,
                 k1: float = 0, k2: float = 0) -> None:
        """
        """
        super().__init__ ( L, E, I )

    #
    @property
    def V0(self) -> float:
        """ Shear force"""
        return (-3 * self.E * self.I / self.L ** 3 * self.Fw
                - 3 * self.E * self.I / self.L ** 2 * self.Ftheta)

    #
    @property
    def theta0(self) -> float:
        """ Slope"""
        return 3 * self.Fw / (2 * self.L) + 0.50 * self.Ftheta
    #
    # torsion
    def phi0(self) -> float:
        """ angle of rotation (radians)
        Left end free to wrap but not twist, right end fixed"""
        phiA_0 = 0
        phiA_I = (self.To/(self.Cw*self.E*self.beta**2)
                  * (C3*Ca4 - C4*Ca3)/(C1*C4 - C2*C3))
        phiA_II = 0
        phiA_III = 0
        return [phiA_0, phiA_I, phiA_II, phiA_III]
    #
    @property
    def T0(self) -> float:
        """ Torsional moment"""
        return -self.To * (C1*Ca4 - C2*Ca3)/(C1*C4 - C2*C3)

#
@dataclass
class PinnedGuided ( SuppBasic ):
    """ """
    __slots__ = [ 'E', 'I', 'FV', 'FM', 'Ftheta', 'Fw' ]

    def __init__(self, L: float, E: float, I: float,
                 k1: float = 0, k2: float = 0) -> float:
        """
        """
        super().__init__ ( L, E, I )

    #
    @property
    def V0(self) -> float:
        """ Shear force"""
        return -self.FV

    #
    @property
    def theta0(self) -> float:
        """ Slope"""
        return L ** 2 / (2 * self.E * self.I) * self.FV - self.Ftheta

    # torsion
    def phi0(self) -> float:
        """ angle of rotation (radians) """
        phiA_0 = (self.To/(self.Cw*self.E*self.beta**3)
                  * (C3*Ca2 / C1) -Ca4)
        phiA_I = 0
        phiA_II = (-self.To/(self.Cw*self.E*self.beta)
                  * (Ca2/C1))
        phiA_III = 0
        return [phiA_0, phiA_I, phiA_II, phiA_III]
    #
    @property
    def T0(self) -> float:
        """ Torsional moment"""
        return 0
#
@dataclass
class PinnedSpring ( SuppBasic ):
    """ """
    __slots__ = [ 'E', 'I', 'FV', 'FM', 'Ftheta', 'Fw', 'k2', 'k1' ]

    def __init__(self, L: float, E: float, I: float,
                 k1: float, k2: float) -> None:
        """
        """
        super().__init__ ( L, E, I )
        # spring stiffness
        self.k1 = k1
        self.k2 = k2

    #
    @property
    def theta0(self) -> float:
        """ Bending moment"""
        A1, A2, A3, A4, A5, A6 = self._spring_2 ( self.k2 )
        return (A1 * self.Fw - self.L * A2 / (6 * self.E * self.I)) / A3

    #
    @property
    def V0(self) -> float:
        """ Shear """
        A1, A2, A3, A4, A5, A6 = self._spring_2 ( self.k2 )
        return ((self.k2 / self.L ** 2) * self.Fw + A2 / self.L) / A3


#
# Fixed
#
@dataclass
class FixedPinned ( SuppBasic ):
    """ """
    __slots__ = [ 'E', 'I', 'FV', 'FM', 'Ftheta', 'Fw' ]

    def __init__(self, L: float, E: float, I: float,
                 k1: float = 0, k2: float = 0) -> None:
        """
        """
        super().__init__ ( L, E, I )

    #
    @property
    def V0(self) -> float:
        """ Shear force"""
        return -3 * self.E * self.I / self.L ** 3 * self.Fw - 3 / (2 * self.L) * self.FM

    #
    @property
    def M0(self) -> float:
        """ Bending moment"""
        return 3 * self.E * self.I / self.L ** 2 * self.Fw + 0.50 * self.FM

    #
    # torsion
    def phi0(self) -> float:
        """ angle of rotation (radians) """
        phiA_0 = 0
        phiA_I = 0
        phiA_II = (self.To/(self.Cw*self.E*self.beta)
                  * (self.beta*self.L*A2 - self.beta*self.L1*C2)/(C1*C4 - C2*C3))
        phiA_III = -self.To / (self.Cw*self.E) * (A2-self.beta*self.L1*C1)/(C1*C4 - C2*C3)
        return [phiA_0, phiA_I, phiA_II, phiA_III]
    #
    @property
    def T0(self) -> float:
        """ Torsional moment"""
        return -self.To - self.To * (C1*Ca4 - C2*Ca3)/(C1*C4 - C2*C3)
#
@dataclass
class FixedFixed( SuppBasic ):
    """ """
    __slots__ = [ 'L','E', 'I', 'FV', 'FM', 'Ftheta', 'Fw', '_response']

    def __init__(self, L: float, E: float, I: float,
                 k1: float = 0, k2: float = 0) -> None:
        """
        """
        super().__init__(L, E, I )

    #
    @property
    def V0(self) -> float:
        """ Shear force"""
        return (-12 * self.E * self.I / self.L ** 3 * self.Fw
                - 6 * self.E * self.I / self.L ** 2 * self.Ftheta)

    #
    @property
    def M0(self) -> float:
        """ Bending moment"""
        return (6 * self.E * self.I / self.L ** 2 * self.Fw
                + 2 * self.E * self.I / self.L * self.Ftheta)

    #
    # torsion
    def phi0(self) -> float:
        """ angle of rotation (radians) """
        phiA_0 = 0
        phiA_I = (self.To/(self.Cw*self.E*self.beta)
                  * (C3*Ca4 - C4*Ca3) / (C2*C4 - C3**2))
        phiA_II = 0
        return [phiA_0, phiA_I, phiA_II]
    #
    @property
    def T0(self) -> float:
        """ Torsional moment"""
        return - self.To * (C2*Ca4 - C3*Ca3) / (C2*C4 - C3**2)
#
@dataclass
class FixedFree ( SuppBasic ):
    """ """
    __slots__ = [ 'E', 'I', 'FV', 'FM', 'Ftheta', 'Fw' ]

    def __init__(self, L: float, E: float, I: float,
                 k1: float = 0, k2: float = 0) -> None:
        """
        """
        super().__init__ ( L, E, I )

    #
    @property
    def V0(self) -> float:
        """ Shear force"""
        return -self.FV

    #
    @property
    def M0(self) -> float:
        """ Bending moment"""
        return self.L * self.FV - self.FM

    #
    # torsion
    def phi0(self) -> float:
        """ angle of rotation (radians) """
        phiA_0 = 0
        phiA_I = (self.To/(self.Cw*self.E*self.beta)
                  * ((A2-C2) / C1 ))
        phiA_II = 0
        return [phiA_0, phiA_I, phiA_II]
    #
    @property
    def T0(self) -> float:
        """ Torsional moment"""
        return -self.To
#
@dataclass
class FixedGuided ( SuppBasic ):
    """ """
    __slots__ = [ 'E', 'I', 'FV', 'FM', 'Ftheta', 'Fw' ]

    def __init__(self, L: float, E: float, I: float,
                 k1: float = 0, k2: float = 0) -> None:
        """
        """
        super().__init__ ( L, E, I )

    #
    @property
    def V0(self) -> float:
        """ Shear force"""
        return -self.FV

    #
    @property
    def M0(self) -> float:
        """ Bending moment"""
        return - self.E * self.I / self.L * self.Ftheta + 0.50 * self.FV * self.L


#
@dataclass
class FixedSpring ( SuppBasic ):
    """ """
    __slots__ = [ 'E', 'I', 'FV', 'FM', 'Ftheta', 'Fw', 'k2', 'k1' ]

    def __init__(self, L: float, E: float, I: float,
                 k1: float, k2: float) -> None:
        """
        """
        super().__init__ ( L, E, I )
        # spring stiffness
        self.k1 = k1
        self.k2 = k2

    #
    @property
    def V0(self) -> float:
        """ Shear """
        A1, A2, A3, A4, A5, A6 = self._spring_2 ( self.k2 )
        return ((-3 * self.E * self.I * A5 / self.L ** 3) * self.Fw + 3 * A2 / (2 * self.L)) / A4

    #
    @property
    def M0(self) -> float:
        """ Bending moment"""
        A1, A2, A3, A4, A5, A6 = self._spring_2 ( self.k2 )
        return ((3 * self.E * self.I / self.L ** 2) * A6 * self.Fw - A2 / 2) / A4


#
# free
#
@dataclass
class FreeFixed ( SuppBasic ):
    """ """
    __slots__ = [ 'E', 'I', 'FV', 'FM', 'Ftheta', 'Fw' ]

    def __init__(self, L: float, E: float, I: float,
                 k1: float = 0, k2: float = 0) -> None:
        """
        """
        super().__init__ ( L, E, I )

    #
    @property
    def theta0(self) -> float:
        """ Bending moment"""
        return -self.Ftheta

    #
    @property
    def w0(self) -> float:
        """ Shear force"""
        return -self.Fw - self.L * self.Ftheta

    #
    # torsion
    def phi0(self) -> float:
        """ angle of rotation (radians) """
        phiA_0 = (self.To/(self.Cw*self.E*self.beta**3)
                  * (C2*Ca3 / C1 - Ca4))
        phiA_I = (self.To/(self.Cw*self.E*self.beta**2)
                  * (Ca3 / C1 ))
        phiA_II = 0
        return [phiA_0, phiA_I, phiA_II]
    #
    @property
    def T0(self) -> float:
        """ Torsional moment"""
        return 0
#
@dataclass
class FreeSpring ( SuppBasic ):
    """ """
    __slots__ = [ 'E', 'I', 'FV', 'FM', 'Ftheta', 'Fw', 'k2', 'k1' ]

    def __init__(self, L: float, E: float, I: float,
                 k1: float, k2: float) -> None:
        """
        """
        super().__init__ ( L, E, I )
        # spring stiffness
        self.k1 = k1
        self.k2 = k2

    #
    @property
    def theta0(self) -> float:
        """ Slope"""
        A1, A2, A3, A4, A5, A6 = self._spring_2 ( self.k2 )
        return -self.Fw - self.L * A2 / self.k2

    #
    @property
    def w0(self) -> float:
        """ Deflection"""
        A1, A2, A3, A4, A5, A6 = self._spring_2 ( self.k2 )
        return - A2 / self.k2


#
# guide
#
@dataclass
class GuidedPinned ( SuppBasic ):
    """ """
    __slots__ = [ 'E', 'I', 'FV', 'FM', 'Ftheta', 'Fw' ]

    def __init__(self, L: float, E: float, I: float,
                 k1: float = 0, k2: float = 0) -> None:
        """
        """
        super().__init__ ( L, E, I )

    #
    @property
    def M0(self) -> float:
        """ Bending moment"""
        return - self.FM

    #
    @property
    def w0(self) -> float:
        """ Shear force"""
        return - self.Fw - 0.50 * self.L ** 2 / (2 * self.E * self.I) * self.FM


#
@dataclass
class GuidedFixed ( SuppBasic ):
    """ """
    __slots__ = [ 'E', 'I', 'FV', 'FM', 'Ftheta', 'Fw' ]

    def __init__(self, L: float, E: float, I: float,
                 k1: float = 0, k2: float = 0) -> None:
        """
        """
        super().__init__ ( L, E, I )

    #
    @property
    def M0(self) -> float:
        """ Bending moment"""
        return - self.E * self.I / self.L * self.Ftheta

    #
    @property
    def w0(self) -> float:
        """ Shear force"""
        return - self.Fw - 0.50 * self.Ftheta * self.L


#
@dataclass
class GuidedSpring ( SuppBasic ):
    """ """
    __slots__ = [ 'E', 'I', 'FV', 'FM', 'Ftheta', 'Fw', 'k2', 'k1' ]

    def __init__(self, L: float, E: float, I: float,
                 k1: float, k2: float) -> None:
        """
        """
        super().__init__ ( L, E, I )
        # spring stiffness
        self.k1 = k1
        self.k2 = k2

    #
    @property
    def theta0(self) -> float:
        """ Slope"""
        A1, A2, A3, A4, A5, A6 = self._spring_2 ( self.k2 )
        return (-A6 * self.Fw + A2 * (self.L ** 2 / (2 * self.E * self.I))) / A6

    #
    @property
    def M0(self) -> float:
        """ Bending moment"""
        A1, A2, A3, A4, A5, A6 = self._spring_2 ( self.k2 )
        return A2 / A6


#
# Spring
#
@dataclass
class SpringPinned ( SuppBasic ):
    """ """
    __slots__ = [ 'E', 'I', 'FV', 'FM', 'Ftheta', 'Fw', 'k1' ]

    def __init__(self, L: float, E: float, I: float,
                 k1: float, k2: float = 0) -> None:
        """
        """
        super().__init__ ( L, E, I )
        # spring stiffness
        self.k1 = k1
        # self.k2 = k2

    #
    @property
    def theta0(self) -> float:
        """ Slope"""
        A7, A8, A9, A10, A11 = self._spring_1 ( self.k1 )
        return (self.Fw / self.L + L * self.FM / (6 * self.E * self.I)) / A7

    #
    @property
    def V0(self) -> float:
        """ Shear """
        A7, A8, A9, A10, A11 = self._spring_1 ( self.k1 )
        return -(self.k1 / self.L ** 2 * self.Fw + self.FM * A8) / A7

    #
    @property
    def M0(self) -> float:
        """ Bending moment"""
        return self.k1 * self.theta0


#
@dataclass
class SpringFixed ( SuppBasic ):
    """ """
    __slots__ = [ 'E', 'I', 'FV', 'FM', 'Ftheta', 'Fw', 'k1' ]

    def __init__(self, L: float, E: float, I: float,
                 k1: float, k2: float = 0) -> None:
        """
        """
        super().__init__ ( L, E, I )
        # spring stiffness
        self.k1 = k1
        # self.k2 = k2

    #
    @property
    def theta0(self) -> float:
        """ Slope"""
        A7, A8, A9, A10, A11 = self._spring_1 ( self.k1 )
        return (3 * self.Fw / (2 * self.L) + self.Ftheta / 2) / A9

    #
    @property
    def V0(self) -> float:
        """ Shear """
        A7, A8, A9, A10, A11 = self._spring_1 ( self.k1 )
        return -(A11 * self.Fw + (3 * self.E * self.I / self.L) * self.Ftheta * A8) / A9

    #
    @property
    def M0(self) -> float:
        """ Bending moment"""
        return self.k1 * self.theta0


#
@dataclass
class SpringFree ( SuppBasic ):
    """ """
    __slots__ = [ 'E', 'I', 'FV', 'FM', 'Ftheta', 'Fw', 'k1' ]

    def __init__(self, L: float, E: float, I: float,
                 k1: float, k2: float = 0) -> None:
        """
        """
        super().__init__ ( L, E, I )
        # spring stiffness
        self.k1 = k1
        # self.k2 = k2

    #
    def theta(self, L: float) -> float:
        """ Slope"""
        # A7,A8,A9,A10,A11 = self._spring_1(self.k1, L)
        return 1 / self.k1 * (self.FM - L * self.FV)

    #
    def V(self, L: float) -> float:
        """ Shear """
        # A7,A8,A9,A10,A11 = self._spring_1(self.k1, L)
        return - self.FV

    #
    def M(self, L: float) -> float:
        """ Bending moment"""
        return self.k1 * self.theta ( L )
    #


@dataclass
class SpringGuided ( SuppBasic ):
    """ """
    __slots__ = [ 'E', 'I', 'FV', 'FM', 'Ftheta', 'Fw', 'k1' ]

    def __init__(self, L: float, E: float, I: float,
                 k1: float, k2: float = 0) -> None:
        """
        """
        super().__init__ ( L, E, I )
        # spring stiffness
        self.k1 = k1
        # self.k2 = k2

    #
    @property
    def theta0(self) -> float:
        """ Slope"""
        A7, A8, A9, A10, A11 = self._spring_1 ( self.k1 )
        return (-self.Ftheta + (self.L ** 2 / (2 * self.E * self.I) * self.FV)) / A10

    #
    @property
    def V0(self) -> float:
        """ Shear """
        # A7,A8,A9,A10,A11 = self._spring_1(self.k1, L)
        return - self.FV

    #
    @property
    def M0(self) -> float:
        """ Bending moment"""
        return self.k1 * self.theta0


#
@dataclass
class SpringSpring ( SuppBasic ):
    """ """
    __slots__ = [ 'E', 'I', 'FV', 'FM', 'Ftheta', 'Fw',
                  'k2', 'k1' ]

    def __init__(self, L: float, E: float, I: float,
                 k1: float, k2: float) -> None:
        """
        """
        super().__init__ ( L, E, I )
        # spring stiffness
        self.k1 = k1
        self.k2 = k2

    #
    @property
    def theta0(self) -> float:
        """ Slope"""
        A1, A2, A3, A4, A5, A6 = self._spring_2 ( self.k2 )
        A12, A13 = self._spring_12 ( self.k1, self.k2, self.L )
        return (A5 / L * self.Fw - self.L * A2 / (6 * self.E * self.I)) / A13

    #
    @property
    def V0(self) -> float:
        """ Shear """
        A1, A2, A3, A4, A5, A6 = self._spring_2 ( self.k2 )
        A7, A8, A9, A10, A11 = self._spring_1 ( self.k1 )
        A12, A13 = self._spring_12 ( self.k1, self.k2, self.L )
        return (A12 * self.Fw + A2 * A8) / A13

    #
    @property
    def M0(self) -> float:
        """ Bending moment"""
        return self.k1 * self.theta0

    #
    def _spring_12(self, k1: float, k2: float, L: float):
        """ """
        A12 = (k2 - k1) / L ** 2 + k2 * k1 / (self.E * self.I * L)
        A13 = (1 + k1 * L / (3 * self.E * self.I)
               - k2 * L / (3 * self.E * self.I)
               - k1 * k2 * L ** 2 / (12 * (self.E * self.I) ** 2))
        return A12, A13
    #


#
class Parameters ( NamedTuple ):
    """ """
    A1: float
    A2: float
    A3: float
    A4: float
    A5: float
    A6: float
    #
    A7: float
    A8: float
    A9: float
    A10: float
    A11: float
    #
    A12: float
    A13: float
#
#