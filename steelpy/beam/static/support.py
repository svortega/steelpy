# 
# Copyright (c) 2019-2020 steelpy
# 

# Python stdlib imports
#from collections.abc import Mapping
from dataclasses import dataclass
from math import factorial
from typing import NamedTuple, Dict, List, Tuple, Union


# package imports
from steelpy.process.units.main import Units


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
#
#
class Support:
    __slots__ = ['cls']
    
    def __init__(self, cls) -> None:
        """
        """
        self.cls = cls      
    #
    def __setitem__(self, support_name: Union[int, int], 
                    parameters: Union[List[float], Dict[str, float]]) -> None:
        """
        Support type : pinned/fixed/free/guided
        Support position along the beam length (optional)
        """
        try:
            self.cls._labels.index(support_name)
            raise Exception('    *** warning support {:} already exist'
                            .format(support_name))
        except ValueError:
            self.cls._labels.append(support_name)
            self.cls._reactions.append(-1)
            if isinstance( parameters, (list, tuple) ):
                self.cls._fixity.append(parameters[0])
                try:
                    self.cls._supcoord.append(parameters[1].value)
                except IndexError:
                    self.cls._supcoord.append(-1)
            elif isinstance( parameters, dict ):
                print('fix here')
            elif isinstance ( parameters, str ):
                self.cls._fixity.append(parameters)
                self.cls._supcoord.append( -1 )
            else:
                raise Exception ( '   *** Support input format not recognized' )

        #print('--')
    #
    def __getitem__(self, support_name: Union[str,int]) -> List[Tuple]:
        """
        """
        index = self.cls._labels.index(support_name)
        if self.cls._reactions[ index ] == -1:
            if self.cls._fixity[index] == "user":
                units = Units ()
                default = ReacInput(R= 0*units.N,
                                    M= 0*units.N*units.m,
                                    theta= 0*units.radians,
                                    w= 0*units.m)
                self.cls._reactions[index] = ReactPlane(default, default)
            else:
                self.solve()
        return self.cls._reactions[index]
    #
    def _get_load(self):
        """ """
        sect = self.cls.section._get_properties()
        I_plane = {"in_plane":sect.Iy, "out_plane":sect.Iy}
        mat = self.cls.material
        E = mat.E.convert("pascal").value
        L = self.cls._beam_length
        loadres = {"in_plane":None, "out_plane":None}
        for plane in loadres.keys():
            V=[]; M=[]; theta=[]; w=[]; x=[]
            try:
                for load in self.cls._load:
                    V.append(load[plane].V(L))
                    M.append(load[plane].M(L))
                    theta.append(load[plane].theta(L, E, I_plane[plane]))
                    w.append(load[plane].w(L, E, I_plane[plane]))
                loadres[plane] = FuncRes(V, M, theta, w, x)
            except AttributeError:
                continue
        return loadres
    #
    def _get_support_1(self, load):
        """ solve support 1 according determined boundary conditions"""
        sect = self.cls.section._get_properties()
        I_plane = {"in_plane":sect.Iy, "out_plane":sect.Iy}
        mat = self.cls.material
        E = mat.E.convert("pascal").value
        #
        supp_1=self.cls._fixity[0]
        try:
            supp_2=self.cls._fixity[1]
        except IndexError:
            supp_2="free"
            self.cls._labels.append(2)
            self.cls._reactions.append(-1)
            self.cls._supcoord.append(-1)
        suppfun = self._supp_func(supp_1, supp_2)
        #
        supp = {"in_plane":None, "out_plane":None}
        for key, item in load.items():
            if not item:
                supp[key] = SuppBasic(E=E, I=I_plane[key],
                                      FV=0, FM=0, Ftheta=0, Fw=0)
                continue
            supp[key] = suppfun(E=E, I=I_plane[key],
                                FV=item.FV, FM=item.FM,
                                Ftheta=item.Ftheta, Fw=item.Fw)
        return supp
    #
    def _get_support_2(self, load, support_1):
        """ solve support n according determined boundary conditions"""
        sect = self.cls.section._get_properties()
        I_plane = {"in_plane":sect.Iy, "out_plane":sect.Iy}
        mat = self.cls.material
        E = mat.E.convert("pascal").value
        L = self.cls._beam_length
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
    def solve(self):
        """ """
        if len(self.cls._reactions) > 2:
            raise RuntimeError("max 2 supports currently available")
        #
        load = self._get_load()
        L = self.cls._beam_length
        #
        supports = [self._get_support_1(load)]
        #
        supports.append(self._get_support_2(load, supports[0]))
        print("---")
        #print(in_plane.FM)
        #print( in_plane.Fw )
        #
        print( "--- R1" )
        print( "R =", supports[0]["in_plane"].V(L) )
        print( "M =", supports[0]["in_plane"].M(L) )
        print( "theta =", supports[0]["in_plane"].theta(L) )
        print( "delta =", supports[0]["in_plane"].w(L) )
        print( "--- R2" )
        print( "R =",supports[1]["in_plane"].V(L))
        print( "M =", supports[1]["in_plane"].M(L))
        print( "theta =", supports[1]["in_plane"].theta(L))
        print( "delta =", supports[1]["in_plane"].w(L))
        #
        units = Units()
        for index, item in enumerate(self.cls._labels):
            #fixity = self.cls._fixity[index]
            in_plane = ReacRes(R= supports[index]["in_plane"].V(L)*units.N,
                               M= supports[index]["in_plane"].M(L)*units.N*units.m,
                               theta= supports[index]["in_plane"].theta(L)*units.radians,
                               w= supports[index]["in_plane"].w(L)*units.m)
            out_plane = ReacRes(R= supports[index]["out_plane"].V(L)*units.N,
                                M= supports[index]["out_plane"].M(L)*units.N*units.m,
                                theta= supports[index]["out_plane"].theta(L)*units.radians,
                                w= supports[index]["out_plane"].w(L) *units.m) 
            self.cls._reactions[ index ] = ReactPlane(in_plane, out_plane)
        #print('--')
    #
    #
    def _supp_func(self, supp1, supp2):
        """ """
        if supp1 == "pinned":
            if supp2 == "pinned":
                return PinnedPinned
            elif supp2 == "fixed":
                return PinnedFixed
            elif supp2 == "guided":
                return PinnedGuided
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
            else:
                raise IOError("boundary no supported")
        
        elif supp1 == "free":
            if supp2 == "fixed":
                return FreeFixed
            else:
                raise IOError("unstable")
        
        elif supp1 == "guided":
            if supp2 == "pinned":
                return GuidedPinned
            elif supp2 == "fixed":
                return GuidedFixed
            else:
                raise IOError("unstable")
        
        else:
            raise IOError("boundary no supported")
#
#
#
@dataclass
class Response:
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
        return -1*(self.V0 + self.FV)
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
    __slots__ = [ 'E', 'I', 'FV', 'FM', 'Ftheta', 'Fw']

    def __init__(self, E: float, I: float,
                 FV: float, FM: float, Ftheta: float, Fw: float) -> None:
        """
        """
        #self.L:float = L
        self.E:float = E
        self.I:float = I
        self.FV:float = FV
        self.FM:float = FM
        self.Ftheta:float = Ftheta
        self.Fw:float = Fw
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
# Pinned
#
@dataclass
class PinnedPinned(SuppBasic):
    """ """
    __slots__ = [ 'E', 'I', 'FV', 'FM', 'Ftheta', 'Fw' ]

    def __init__(self, E: float, I: float,
                 FV: float, FM: float, Ftheta: float, Fw: float) -> float:
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
                 FV: float, FM: float, Ftheta: float, Fw: float) -> None:
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
                 FV: float, FM: float, Ftheta: float, Fw: float) -> float:
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
# Fixed
#
@dataclass
class FixedPinned(SuppBasic):
    """ """
    __slots__ = [ 'E', 'I', 'FV', 'FM', 'Ftheta', 'Fw' ]

    def __init__(self, E: float, I: float,
                 FV: float, FM: float, Ftheta: float, Fw: float) -> None:
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
                 FV: float, FM: float, Ftheta: float, Fw: float) -> None:
        """
        """
        SuppBasic.__init__( self, E, I, FV, FM, Ftheta, Fw )

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
                 FV: float, FM: float, Ftheta: float, Fw: float) -> None:
        """
        """
        SuppBasic.__init__( self, E, I, FV, FM, Ftheta, Fw )

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
                 FV: float, FM: float, Ftheta: float, Fw: float) -> None:
        """
        """
        self.L = L
        SuppBasic.__init__( self, E, I, FV, FM, Ftheta, Fw )
    #
    def V(self, L:float) -> float:
        """ Shear force"""
        return -self.FV
    #
    def M(self, L:float) -> float:
        """ Bending moment"""
        return - self.E*self.I/L * self.Ftheta + 0.50 * self.FV * L
#
# free
#
@dataclass
class FreeFixed(SuppBasic):
    """ """
    __slots__ = [ 'E', 'I', 'FV', 'FM', 'Ftheta', 'Fw' ]

    def __init__(self, E: float, I: float,
                 FV: float, FM: float, Ftheta: float, Fw: float) -> None:
        """
        """
        SuppBasic.__init__( self, E, I, FV, FM, Ftheta, Fw )

    #
    def theta(self, L:float) -> float:
        """ Bending moment"""
        return -self.Ftheta
    
    #
    def w(self, L:float) -> float:
        """ Shear force"""
        return -self.Fw - L * self.Ftheta    
#
# guide
#
@dataclass
class GuidedPinned(SuppBasic):
    """ """
    __slots__ = [ 'E', 'I', 'FV', 'FM', 'Ftheta', 'Fw' ]

    def __init__(self, E: float, I: float,
                 FV: float, FM: float, Ftheta: float, Fw: float) -> None:
        """
        """
        self.L = L
        SuppBasic.__init__( self, E, I, FV, FM, Ftheta, Fw )
    #
    def M(self, L:float) -> float:
        """ Bending moment"""
        return - self.FM
    
    #
    def w(self, L:float) -> float:
        """ Shear force"""
        return - self.Fw -  0.50  * L**2 / (2*self.E*self.I) * self.FM
#
#
@dataclass
class GuidedFixed(SuppBasic):
    """ """
    __slots__ = [ 'E', 'I', 'FV', 'FM', 'Ftheta', 'Fw' ]

    def __init__(self, E: float, I: float,
                 FV: float, FM: float, Ftheta: float, Fw: float) -> None:
        """
        """
        self.L = L
        SuppBasic.__init__( self, E, I, FV, FM, Ftheta, Fw )
    #
    def M(self, L:float) -> float:
        """ Bending moment"""
        return - self.E*self.I/L * self.Ftheta
    
    #
    def w(self, L:float) -> float:
        """ Shear force"""
        return - self.Fw -  0.50 * self.Ftheta * L  
#
#
