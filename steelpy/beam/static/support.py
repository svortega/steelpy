# 
# Copyright (c) 2019-2021 steelpy
# 

# Python stdlib imports
from collections.abc import Mapping
from dataclasses import dataclass
from math import factorial
from typing import NamedTuple, Dict, List, Tuple, Union, Iterator


# package imports
from steelpy.process.units.main import Units


#
#@dataclass
class ReacRes(NamedTuple):
    R:Units
    M:Units
    theta:Units
    w:Units
    #
    #R_outplane: Units
    #M_outplane:Units
    #theta_outplane:Units
    #w_outplane:Units
    def __str__(self) -> str:
        """ """
        output = ""
        output += "React  [N  ]  {: 1.4E}\n".format(self.R.value)
        output += "Moment [N*m]  {: 1.4E}\n".format(self.M.value)
        output += "theta  [rad]  {: 1.4E}\n".format(self.theta.value)
        output += "delta  [m  ]  {: 1.4E}\n".format(self.w.value)          
        return output
#
@dataclass
class ReacInput:
    R:Units
    M:Units
    theta:Units
    #R_outplane: Units
    #M_outplane:Units
    #theta_outplane:Units
#
@dataclass
class ReactPlane:
    in_plane:Union[Tuple,ReacInput]
    out_plane:Union[Tuple,ReacInput]
    #
    @property
    def y(self):
        return self.in_plane
    
    @property
    def z(self):
        return self.out_plane
    #
    def __str__(self) -> str:
        """ """
        output = ""
        output += "        Unit   In Plane    Out Plane\n"
        output += "React  [N  ]  {: 1.4E} {: 1.4E}\n".format(self.in_plane.R.convert("newton").value,
                                                            self.out_plane.R.convert("newton").value)
        output += "Moment [Nm ]  {: 1.4E} {: 1.4E}\n".format(self.in_plane.M.convert("newton*metre").value,
                                                             self.out_plane.M.convert("newton*metre").value)
        output += "theta  [rad]  {: 1.4E} {: 1.4E}\n".format(self.in_plane.theta.value,
                                                             self.out_plane.theta.value)
        output += "delta  [m  ]  {: 1.4E} {: 1.4E}\n".format(self.in_plane.w.value,
                                            self.out_plane.w.value)
        return output
#       
class FuncRes(NamedTuple):
    V :List[float]
    M :List[float]
    theta :List[float]
    w :List[float]
    x :List[float]
    #
    @property
    def FV(self) -> float:
        return sum(self.V)

    @property
    def FM(self) -> float:
        return sum(self.M)

    @property
    def Ftheta(self) -> float:
        return sum(self.theta)

    @property
    def Fw(self) -> float:
        return sum(self.w)
#
#
class Support(Mapping):
    __slots__ = ['cls', '_fixity', '_reactions', 
                 '_supcoord', '_labels', '_k']
    
    def __init__(self, cls) -> None:
        """
        """
        self.cls = cls
        #
        self._fixity:List = []
        self._reactions:List = []
        self._supcoord:List = []
        self._labels:List = []
        #
        self._k:List = []
    #
    def __setitem__(self, support_name: Union[int, int], 
                    parameters: Union[List[float], Dict[str, float]]) -> None:
        """
        Support type : pinned/fixed/free/guided/spring
        Support position along the beam length (optional)
        """
        try:
            self._labels.index(support_name)
            raise Exception('    *** warning support {:} already exist'
                            .format(support_name))
        except ValueError:
            self._labels.append(support_name)
            self._reactions.append(-1)
            self._k.append({"in_plane":0, "out_plane":0})
            if isinstance(parameters, (list, tuple)):
                self._fixity.append(parameters[0])
                try:
                    item = parameters[1]
                    if parameters[0] == 'spring':
                        item = item.convert("newton/metre").value
                        self._k[-1]["in_plane"] = item
                        try:
                            item2 = parameters[2].convert("newton/metre").value
                            self._k[-1]["out_plane"] = item2
                        except IndexError:
                            self._k[-1]["out_plane"] = item
                    else:
                        self._supcoord.append(item.value)                      
                except IndexError:
                    self._supcoord.append(-1)
            elif isinstance( parameters, dict ):
                print('fix here')
            elif isinstance(parameters, str):
                self._fixity.append(parameters)
                self._supcoord.append(-1)
            else:
                raise Exception('   *** Support input format not recognized')
        #print('--')
    #
    def __getitem__(self, support_name: Union[str,int]) -> List[Tuple]:
        """
        """
        index = self._labels.index(support_name)
        if self._reactions[index] == -1:
            if self._fixity[index] == "user":
                units = Units ()
                default = ReacInput(R= 0*units.N,
                                    M= 0*units.N*units.m,
                                    theta= 0*units.radians,
                                    w= 0*units.m)
                self._reactions[index] = ReactPlane(default, default)
            else:
                self.__call__()
        return self._reactions[index]
    #
    def _get_load(self):
        """ """
        L, I_plane, E = self._get_properties()
        loadres = {"in_plane":None, "out_plane":None}
        for plane in loadres.keys():
            V=[]; M=[]; theta=[]; w=[]; x=[]
            try:
                for load in self.cls.load:
                    V.append(load[plane].V(L))
                    M.append(load[plane].M(L))
                    theta.append(load[plane].theta(L, E, I_plane[plane]))
                    w.append(load[plane].w(L, E, I_plane[plane]))
                loadres[plane] = FuncRes(V, M, theta, w, x)
            except AttributeError:
                continue
        return loadres
    #
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
    def _get_support_1(self, load):
        """ solve support 1 according determined boundary conditions"""
        L, I_plane, E = self._get_properties()
        #
        supp_1 = self._fixity[0]
        try:
            supp_2 = self._fixity[1]
        except IndexError:
            supp_2 = "free"
            self._labels.append(2)
            self._reactions.append(-1)
            self._supcoord.append(-1)
        suppfun = self._support_func(supp_1, supp_2)
        #
        supp = {"in_plane":None, "out_plane":None}
        for key, item in load.items():
            if not item:
                supp[key] = SuppBasic(E=E, I=I_plane[key],
                                      FV=0, FM=0, Ftheta=0, Fw=0)
                continue
            try:
                supp[key] = suppfun(E=E, I=I_plane[key],
                                    FV=item.FV, FM=item.FM,
                                    Ftheta=item.Ftheta, Fw=item.Fw,
                                    k=self._k[0][key])
            except TypeError: # Two springs data
                supp[key] = suppfun(E=E, I=I_plane[key],
                                    FV=item.FV, FM=item.FM,
                                    Ftheta=item.Ftheta, Fw=item.Fw,
                                    k1=self._k[0][key], k2=self._k[1][key])                
        return supp
    #
    def _get_support_2(self, load, support_1):
        """ solve support n according determined boundary conditions"""
        L, I_plane, E = self._get_properties()
        #
        supp = {"in_plane":None, "out_plane":None}
        for key, item in load.items():
            if not item:
                supp[key] = SuppBasic(E=E, I=I_plane[key],
                                      FV=0, FM=0, Ftheta=0, Fw=0)
                continue
            supp[key] = Response(E=E, I=I_plane[key])
            supp[key].load(FV=item.FV, FM=item.FM,
                           Ftheta=item.Ftheta, Fw=item.Fw)
            supp[key].reacctions(support_1[key].V(L), support_1[key].M(L), 
                                 support_1[key].theta(L), support_1[key].w(L))
        return supp
    #
    def __call__(self):
        """ """
        if len(self._reactions) > 2:
            raise RuntimeError("max 2 supports currently available")
        #
        load = self._get_load()
        L = self.cls.beam_length
        #
        supports = [self._get_support_1(load)]
        #
        supports.append(self._get_support_2(load, supports[0]))
        #
        units = Units()
        for index, item in enumerate(self._labels):
            #fixity = self.cls._fixity[index]
            in_plane = ReacRes(R= supports[index]["in_plane"].V(L)*units.N,
                               M= supports[index]["in_plane"].M(L)*units.N*units.m,
                               theta= supports[index]["in_plane"].theta(L)*units.radians,
                               w= supports[index]["in_plane"].w(L)*units.m)
            #
            out_plane = ReacRes(R= supports[index]["out_plane"].V(L)*units.N,
                                M= supports[index]["out_plane"].M(L)*units.N*units.m,
                                theta= supports[index]["out_plane"].theta(L)*units.radians,
                                w= supports[index]["out_plane"].w(L) *units.m) 
            #
            self._reactions[index] = ReactPlane(in_plane, out_plane)
        #print('--')
        return self._reactions
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
                raise IOError("unstable")
        
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
                raise IOError("boundary no supported")
        
        elif supp1 == "free":
            if supp2 == "fixed":
                return FreeFixed
            elif supp2 == "spring":
                return FreeSpring
            else:
                raise IOError("unstable")
        
        elif supp1 == "guided":
            if supp2 == "pinned":
                return GuidedPinned
            elif supp2 == "fixed":
                return GuidedFixed
            elif supp2 == "spring":
                return GuidedSpring
            else:
                raise IOError("unstable")
        
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
                raise IOError("boundary no supported")
        else:
            raise IOError("boundary no supported")
    #
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
    def __str__(self) -> str:
        """ """
        output = "------- Reactions Results\n"
        output += "\n"
        _reactions = self.__call__()
        for x, support_name in enumerate(self._labels):
            try:
                output += "Support {:} : {:}\n".format(support_name, self._fixity[x])
                output += _reactions[x].__str__()
                output += "\n"
            except IndexError:
                continue
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
                  'V0', 'M0', 'theta0', 'w0']

    def __init__(self, E: float, I: float) -> None:
        """
        """
        self.E = E
        self.I = I
    #
    def load(self, FV: float, FM: float, 
             Ftheta: float, Fw: float) -> None:
        """ """
        self.FV = FV
        self.FM = FM
        self.Ftheta = Ftheta
        self.Fw = Fw
    #
    def reacctions(self,V0:float, M0:float, 
                   theta0:float, w0:float) -> None:
        """ """
        self.V0=V0
        self.M0=M0
        self.theta0=theta0
        self.w0=w0
    #
    def V(self, x:float) -> float:
        """ Shear force"""
        return self.V0 + self.FV
    #
    def M(self, x:float) -> float:
        """ Bending moment"""
        return self.M0 + self.V0 * x + self.FM
    #
    def theta(self, x:float) -> float:
        """ Slope"""
        return (self.theta0 + self.V0* x**2/(2*self.E*self.I)
                + self.M0 * x/(self.E*self.I) + self.Ftheta)
    #
    def w(self, x:float) -> float:
        """ Deflection"""
        return (self.w0 - self.theta0 * x
                - self.V0 * x**3 / (factorial(3) * self.E * self.I)
                - self.M0 * x**2 / (2 * self.E * self.I) + self.Fw)
#
#
@dataclass
class SuppBasic:
    __slots__ = ['E', 'I', 'FV', 'FM', 'Ftheta', 'Fw']

    def __init__(self, E: float, I: float,
                 FV: float, FM: float, Ftheta: float, Fw: float) -> None:
        """
        """
        self.E:float = E
        self.I:float = I
        self.FV:float = FV
        self.FM:float = FM
        self.Ftheta:float = Ftheta
        self.Fw:float = Fw
        # spring stiffness
        #self.k1 = k1
        #self.k2 = k2
    #
    def V(self, L:float) -> float:
        """ Shear """
        return 0     
    #
    def M(self, L:float) -> float:
        """ Bending moment"""
        return 0
    #
    def theta(self, L:float):
        """ Slope"""
        return 0
    #
    def w(self, L:float) -> float:
        """ Deflection"""
        return 0
    #
    #
    def _spring_1(self, k1:float, L:float):
        """ """
        #
        A7 = 1 + k1 * L / (3*self.E*self.I)
        A8 = 1/L + k1 / (2*self.E*self.I)
        A9 = 1 + k1 * L / (4*self.E*self.I)
        A10 = 1 + k1 * L / (self.E*self.I)
        A11 = 3*self.E*self.I /L**3  + 3*k1 /L**2
        return  A7,A8,A9,A10,A11
    #
    def _spring_2(self, k2:float, L:float):
        """ """
        A1 = 1/L - k2 / (2*self.E*self.I)
        A2 = k2 * self.Ftheta - self.FM
        A3 = 1 - k2 * L / (3*self.E*self.I)
        A4 = 1 - k2 * L / (4*self.E*self.I)
        A5 = 1 - k2 * L / (2*self.E*self.I)
        A6 = 1 - k2 * L / (self.E*self.I)
        return A1, A2, A3, A4, A5, A6
    #    
#
# Pinned
#
@dataclass
class PinnedPinned(SuppBasic):
    """ """
    __slots__ = [ 'E', 'I', 'FV', 'FM', 'Ftheta', 'Fw' ]

    def __init__(self, E: float, I: float,
                 FV: float, FM: float, Ftheta: float, Fw: float,
                 k:float=0) -> float:
        """
        """
        SuppBasic.__init__(self, E, I, FV, FM, Ftheta, Fw)
    #
    #@property
    def V(self, L:float) -> float:
        """ Shear force"""
        return -self.FM/L
    #
    #@property
    def theta(self, L:float) -> float:
        """ Slope"""
        return self.Fw/L + self.FM * L/(6*self.E*self.I)
#
@dataclass
class PinnedFixed(SuppBasic):
    """ """
    __slots__ = [ 'E', 'I', 'FV', 'FM', 'Ftheta', 'Fw' ]

    def __init__(self, E: float, I: float,
                 FV: float, FM: float, Ftheta: float, Fw: float,
                 k:float=0) -> None:
        """
        """
        SuppBasic.__init__(self, E, I, FV, FM, Ftheta, Fw)
    #
    #@property
    def V(self, L:float) -> float:
        """ Shear force"""
        return (-3*self.E*self.I/L**3 * self.Fw
                - 3*self.E*self.I/L**2  * self.Ftheta)
    #
    #@property
    def theta(self, L:float) -> float:
        """ Slope"""
        return 3*self.Fw/(2*L) + 0.50 * self.Ftheta
#
@dataclass
class PinnedGuided(SuppBasic):
    """ """
    __slots__ = [ 'E', 'I', 'FV', 'FM', 'Ftheta', 'Fw' ]

    def __init__(self, E: float, I: float,
                 FV: float, FM: float, Ftheta: float, Fw: float,
                 k:float=0) -> float:
        """
        """
        SuppBasic.__init__(self, E, I, FV, FM, Ftheta, Fw)
    #
    #@property
    def V(self, L:float) -> float:
        """ Shear force"""
        return -self.FV
    #
    #@property
    def theta(self, L:float) -> float:
        """ Slope"""
        return L**2/(2*self.E*self.I) * self.FV - self.Ftheta
#
@dataclass
class PinnedSpring(SuppBasic):
    """ """
    __slots__ = ['E', 'I', 'FV', 'FM', 'Ftheta', 'Fw', 'k2', 'k1']

    def __init__(self, E: float, I: float,
                 FV: float, FM: float, Ftheta: float, Fw: float,
                 k1:float, k2:float) -> None:
        """
        """
        SuppBasic.__init__(self, E, I, FV, FM, Ftheta, Fw)
        # spring stiffness
        self.k1 = k1
        self.k2 = k2
    #
    def theta(self, L:float) -> float:
        """ Bending moment"""
        A1, A2, A3, A4, A5, A6 = self._spring_2(self.k2, L)
        return (A1*self.Fw - L*A2/(6*self.E*self.I))/A3
    
    #
    def V(self, L:float) -> float:
        """ Shear """
        A1, A2, A3, A4, A5, A6 = self._spring_2(self.k2, L)
        return ((self.k2/L**2)* self.Fw + A2/L)/A3
#
# Fixed
#
@dataclass
class FixedPinned(SuppBasic):
    """ """
    __slots__ = [ 'E', 'I', 'FV', 'FM', 'Ftheta', 'Fw' ]

    def __init__(self, E: float, I: float,
                 FV: float, FM: float, Ftheta: float, Fw: float,
                 k:float=0) -> None:
        """
        """
        SuppBasic.__init__(self, E, I, FV, FM, Ftheta, Fw)
    #
    def V(self, L:float) -> float:
        """ Shear force"""
        return -3*self.E*self.I/L**3 * self.Fw - 3/(2*L) * self.FM
    #
    def M(self, L:float) -> float:
        """ Bending moment"""
        return 3*self.E*self.I/L**2 * self.Fw + 0.50*self.FM
#
@dataclass
class FixedFixed(SuppBasic):
    """ """
    __slots__ = [ 'E', 'I', 'FV', 'FM', 'Ftheta', 'Fw' ]

    def __init__(self, E: float, I: float,
                 FV: float, FM: float, Ftheta: float, Fw: float,
                 k:float=0) -> None:
        """
        """
        SuppBasic.__init__(self, E, I, FV, FM, Ftheta, Fw)

    #
    def V(self, L:float) -> float:
        """ Shear force"""
        return (-12 * self.E * self.I / L**3 * self.Fw
                - 6 * self.E * self.I / L**2 * self.Ftheta)

    #
    def M(self, L:float) -> float:
        """ Bending moment"""
        return (6 * self.E * self.I / L**2 * self.Fw
                + 2 * self.E * self.I / L * self.Ftheta)

#
@dataclass
class FixedFree(SuppBasic):
    """ """
    __slots__ = [ 'E', 'I', 'FV', 'FM', 'Ftheta', 'Fw' ]

    def __init__(self, E: float, I: float,
                 FV: float, FM: float, Ftheta: float, Fw: float,
                 k:float=0) -> None:
        """
        """
        SuppBasic.__init__( self, E, I, FV, FM, Ftheta, Fw)

    #
    def V(self, L:float) -> float:
        """ Shear force"""
        return -self.FV
    #
    def M(self, L:float) -> float:
        """ Bending moment"""
        return L * self.FV - self.FM
#
@dataclass
class FixedGuided(SuppBasic):
    """ """
    __slots__ = [ 'E', 'I', 'FV', 'FM', 'Ftheta', 'Fw' ]

    def __init__(self, E: float, I: float,
                 FV: float, FM: float, Ftheta: float, Fw: float,
                 k:float=0) -> None:
        """
        """
        SuppBasic.__init__( self, E, I, FV, FM, Ftheta, Fw)
    #
    def V(self, L:float) -> float:
        """ Shear force"""
        return -self.FV
    #
    def M(self, L:float) -> float:
        """ Bending moment"""
        return - self.E*self.I/L * self.Ftheta + 0.50 * self.FV * L
#
@dataclass
class FixedSpring(SuppBasic):
    """ """
    __slots__ = ['E', 'I', 'FV', 'FM', 'Ftheta', 'Fw', 'k2', 'k1']

    def __init__(self, E: float, I: float,
                 FV: float, FM: float, Ftheta: float, Fw: float,
                 k1:float, k2:float) -> None:
        """
        """
        SuppBasic.__init__(self, E, I, FV, FM, Ftheta, Fw)
        # spring stiffness
        self.k1 = k1
        self.k2 = k2
    #
    def V(self, L:float) -> float:
        """ Shear """
        A1, A2, A3, A4, A5, A6 = self._spring_2(self.k2, L)
        return ((-3*self.E*self.I*A5 / L**3) * self.Fw + 3*A2 / (2*L))/A4
    
    #
    def M(self, L:float) -> float:
        """ Bending moment"""
        A1, A2, A3, A4, A5, A6 = self._spring_2(self.k2, L)
        return ((3*self.E*self.I / L**2 )* A6 * self.Fw - A2/2) / A4

#
# free
#
@dataclass
class FreeFixed(SuppBasic):
    """ """
    __slots__ = [ 'E', 'I', 'FV', 'FM', 'Ftheta', 'Fw' ]

    def __init__(self, E: float, I: float,
                 FV: float, FM: float, Ftheta: float, Fw: float,
                 k:float=0) -> None:
        """
        """
        SuppBasic.__init__( self, E, I, FV, FM, Ftheta, Fw)

    #
    def theta(self, L:float) -> float:
        """ Bending moment"""
        return -self.Ftheta
    
    #
    def w(self, L:float) -> float:
        """ Shear force"""
        return -self.Fw - L * self.Ftheta    
#
@dataclass
class FreeSpring(SuppBasic):
    """ """
    __slots__ = ['E', 'I', 'FV', 'FM', 'Ftheta', 'Fw', 'k2', 'k1']

    def __init__(self, E: float, I: float,
                 FV: float, FM: float, Ftheta: float, Fw: float,
                 k1:float, k2:float) -> None:
        """
        """
        SuppBasic.__init__(self, E, I, FV, FM, Ftheta, Fw)
        # spring stiffness
        self.k1 = k1
        self.k2 = k2
    #
    def theta(self, L:float) -> float:
        """ Slope"""
        A1, A2, A3, A4, A5, A6 = self._spring_2(self.k2, L)
        return -self.Fw - L*A2/self.k2
    
    #
    def w(self, L:float) -> float:
        """ Deflection"""
        A1, A2, A3, A4, A5, A6 = self._spring_2(self.k2, L)
        return - A2/self.k2
#
# guide
#
@dataclass
class GuidedPinned(SuppBasic):
    """ """
    __slots__ = [ 'E', 'I', 'FV', 'FM', 'Ftheta', 'Fw' ]

    def __init__(self, E: float, I: float,
                 FV: float, FM: float, Ftheta: float, Fw: float,
                 k:float=0) -> None:
        """
        """
        SuppBasic.__init__( self, E, I, FV, FM, Ftheta, Fw)
    #
    def M(self, L:float) -> float:
        """ Bending moment"""
        return - self.FM
    
    #
    def w(self, L:float) -> float:
        """ Shear force"""
        return - self.Fw -  0.50  * L**2 / (2*self.E*self.I) * self.FM
#
@dataclass
class GuidedFixed(SuppBasic):
    """ """
    __slots__ = [ 'E', 'I', 'FV', 'FM', 'Ftheta', 'Fw' ]

    def __init__(self, E: float, I: float,
                 FV: float, FM: float, Ftheta: float, Fw: float,
                 k:float=0) -> None:
        """
        """
        SuppBasic.__init__( self, E, I, FV, FM, Ftheta, Fw)
    #
    def M(self, L:float) -> float:
        """ Bending moment"""
        return - self.E*self.I/L * self.Ftheta
    
    #
    def w(self, L:float) -> float:
        """ Shear force"""
        return - self.Fw -  0.50 * self.Ftheta * L  
#
@dataclass
class GuidedSpring(SuppBasic):
    """ """
    __slots__ = ['E', 'I', 'FV', 'FM', 'Ftheta', 'Fw', 'k2', 'k1']

    def __init__(self, E: float, I: float,
                 FV: float, FM: float, Ftheta: float, Fw: float,
                 k1:float, k2:float) -> None:
        """
        """
        SuppBasic.__init__(self, E, I, FV, FM, Ftheta, Fw)
        # spring stiffness
        self.k1 = k1
        self.k2 = k2
    #
    def theta(self, L:float) -> float:
        """ Slope"""
        A1, A2, A3, A4, A5, A6 = self._spring_2(self.k2, L)
        return (-A6 * self.Fw + A2* (L**2 /(2*self.E*self.I))) / A6
    
    #
    def M(self, L:float) -> float:
        """ Bending moment"""
        A1, A2, A3, A4, A5, A6 = self._spring_2(self.k2, L)
        return A2 / A6
#
# Spring
#
@dataclass
class SpringPinned(SuppBasic):
    """ """
    __slots__ = ['E', 'I', 'FV', 'FM', 'Ftheta', 'Fw', 'k1']

    def __init__(self, E: float, I: float,
                 FV: float, FM: float, Ftheta: float, Fw: float,
                 k:float) -> None:
        """
        """
        SuppBasic.__init__(self, E, I, FV, FM, Ftheta, Fw)
        # spring stiffness
        self.k1 = k
        #self.k2 = k2    
    #
    def theta(self, L:float) -> float:
        """ Slope"""
        A7,A8,A9,A10,A11 = self._spring_1(self.k1, L)
        return (self.Fw/L + L*self.FM/(6*self.E*self.I))/A7
    
    #
    def V(self, L:float) -> float:
        """ Shear """
        A7,A8,A9,A10,A11 = self._spring_1(self.k1, L)
        return -(self.k1/L**2 * self.Fw + self.FM * A8)/A7
    #
    def M(self, L:float) -> float:
        """ Bending moment"""
        return self.k1 * self.theta(L)
#
@dataclass
class SpringFixed(SuppBasic):
    """ """
    __slots__ = ['E', 'I', 'FV', 'FM', 'Ftheta', 'Fw', 'k1']

    def __init__(self, E: float, I: float,
                 FV: float, FM: float, Ftheta: float, Fw: float,
                 k:float) -> None:
        """
        """
        SuppBasic.__init__(self, E, I, FV, FM, Ftheta, Fw)
        # spring stiffness
        self.k1 = k
        #self.k2 = k2
    #
    def theta(self, L:float) -> float:
        """ Slope"""
        A7,A8,A9,A10,A11 = self._spring_1(self.k1, L)
        return (3*self.Fw/(2*L) + self.Ftheta/2) / A9
    
    #
    def V(self, L:float) -> float:
        """ Shear """
        A7,A8,A9,A10,A11 = self._spring_1(self.k1, L)
        return -(A11 * self.Fw + (3*self.E*self.I/L) * self.Ftheta * A8)/A9
    #
    def M(self, L:float) -> float:
        """ Bending moment"""
        return self.k1 * self.theta(L)    
#
@dataclass
class SpringFree(SuppBasic):
    """ """
    __slots__ = ['E', 'I', 'FV', 'FM', 'Ftheta', 'Fw', 'k1']

    def __init__(self, E: float, I: float,
                 FV: float, FM: float, Ftheta: float, Fw: float,
                 k:float) -> None:
        """
        """
        SuppBasic.__init__(self, E, I, FV, FM, Ftheta, Fw)
        # spring stiffness
        self.k1 = k
        #self.k2 = k2    
    #
    def theta(self, L:float) -> float:
        """ Slope"""
        #A7,A8,A9,A10,A11 = self._spring_1(self.k1, L)
        return 1/self.k1 * (self.FM - L*self.FV)
    
    #
    def V(self, L:float) -> float:
        """ Shear """
        #A7,A8,A9,A10,A11 = self._spring_1(self.k1, L)
        return - self.FV
    #
    def M(self, L:float) -> float:
        """ Bending moment"""
        return self.k1 * self.theta(L)    
#
@dataclass
class SpringGuided(SuppBasic):
    """ """
    __slots__ = ['E', 'I', 'FV', 'FM', 'Ftheta', 'Fw', 'k1']

    def __init__(self, E: float, I: float,
                 FV: float, FM: float, Ftheta: float, Fw: float,
                 k:float) -> None:
        """
        """
        SuppBasic.__init__(self, E, I, FV, FM, Ftheta, Fw)
        # spring stiffness
        self.k1 = k
        #self.k2 = k2
    #
    def theta(self, L:float) -> float:
        """ Slope"""
        A7,A8,A9,A10,A11 = self._spring_1(self.k1, L)
        return (-self.Ftheta + (L**2/(2*self.E*self.I) * self.FV ))/ A10
    
    #
    def V(self, L:float) -> float:
        """ Shear """
        #A7,A8,A9,A10,A11 = self._spring_1(self.k1, L)
        return - self.FV
    #
    def M(self, L:float) -> float:
        """ Bending moment"""
        return self.k1 * self.theta(L)    
#
@dataclass
class SpringSpring(SuppBasic):
    """ """
    __slots__ = ['E', 'I', 'FV', 'FM', 'Ftheta', 'Fw',
                 'k2', 'k1']

    def __init__(self, E: float, I: float,
                 FV: float, FM: float, Ftheta: float, Fw: float,
                 k1:float, k2:float) -> None:
        """
        """
        SuppBasic.__init__(self, E, I, FV, FM, Ftheta, Fw)
        # spring stiffness
        self.k1 = k1
        self.k2 = k2
    #
    def theta(self, L:float) -> float:
        """ Slope"""
        A1, A2, A3, A4, A5, A6 = self._spring_2(self.k2, L)
        A12, A13 = self._spring_12(self.k1, self.k2, L)
        return (A5/L * self.Fw - L*A2 / (6*self.E*self.I))/A13
    
    #
    def V(self, L:float) -> float:
        """ Shear """
        A1, A2, A3, A4, A5, A6 = self._spring_2(self.k2, L)
        A7,A8,A9,A10,A11 = self._spring_1(self.k1, L)
        A12, A13 = self._spring_12(self.k1, self.k2, L)
        return (A12 * self.Fw + A2 * A8) / A13
    #
    def M(self, L:float) -> float:
        """ Bending moment"""
        return self.k1 * self.theta(L)
    #
    def _spring_12(self, k1:float, k2:float, L:float):
        """ """
        A12 = (k2 - k1)/L**2 +  k2*k1 / (self.E*self.I*L)
        A13 = (1 + k1 * L / (3*self.E*self.I) 
               - k2 * L / (3*self.E*self.I) 
               -  k1 * k2 * L**2 / (12* (self.E*self.I)**2))
        return  A12, A13    
#
#
class Parameters(NamedTuple):
    """ """
    A1:float
    A2:float
    A3:float
    A4:float
    A5:float
    A6:float
    #
    A7:float
    A8:float
    A9:float
    A10:float
    A11:float
    #
    A12:float
    A13:float    
#
#
