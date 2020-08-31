#
# Copyright (c) 2009-2020 fem2ufo
# 

# Python stdlib imports
from array import array
from collections.abc import Mapping
#from dataclasses import dataclass
from typing import NamedTuple, Tuple, List, Union, Iterable, Dict  


# package imports
from steelpy.trave3D.processor.operations import trns_3Dv
from steelpy.trave3D.preprocessor.assemble import Rmatrix
from steelpy.f2uModel.load.operations import (NodeLoadMaster, 
                                              check_list_units, 
                                              check_list_number,
                                              check_beam_dic)

#
# ---------------------------------
#
class LineBeam(NamedTuple):
    """
    """
    qx1: float
    qy1: float
    qz1: float
    #
    qx2: float
    qy2: float
    qz2: float
    #
    L1: float
    L2: float
    #
    number: int
    load_name: str
    system:str
    load_complex:int
#
class PointBeam(NamedTuple):
    """
    """
    fx: float
    fy: float
    fz: float
    mx: float
    my: float
    mz: float
    distance:float
    number: int
    load_name: str
    system:str
    load_complex:int
#
class BeamDistributed(Mapping):
    """
    """
    __slots__ = ['_type', '_labels', '_index', '_complex',
                 '_L1', '_qx1', '_qy1', '_qz1', 
                 '_L2', '_qx2', '_qy2', '_qz2',
                 '_system', '_system_flag', '_title']

    def __init__(self) -> None:
        """
        """
        # ens 1
        self._L1: array = array("f", [])
        self._qx1: array = array("f", [])
        self._qy1: array = array("f", [])
        self._qz1: array = array("f", [])
        # end 2
        self._L2: array = array("f", [])
        self._qx2: array = array("f", [])
        self._qy2: array = array("f", [])
        self._qz2: array = array("f", [])
        #
        self._labels: List[Union[str, int]] = []
        self._title: List[str] = []
        self._index: int
        self._complex: array = array("I", [])
        # 0-global/ 1-local
        self._system_flag:int = 0
        self._system: array = array("I", [])
    #
    def __setitem__(self, element_name: Union[int, str], 
                    udl: Union[List[float], Dict[str,float]]) -> None:
        """
        """
        self._labels.append(element_name)
        self._title.append(element_name)
        self._index = len(self._labels)-1
        self._system.append(self._system_flag)
        self._complex.append(0)
        #
        # update inputs
        udl = self._get_line_load(udl)
        # end 1
        self._qx1.append(udl[0])
        self._qy1.append(udl[1])
        self._qz1.append(udl[2])
        # end 2
        self._qx2.append(udl[3])
        self._qy2.append(udl[4])
        self._qz2.append(udl[5])
        # distance from ends
        self._L1.append(udl[6])
        self._L2.append(udl[7])
    #
    def __getitem__(self, element_name: Union[int, str]) -> List[Tuple]:
        """
        """
        _index_list: List = [x for x, _item in enumerate(self._labels)
                             if _item == element_name]
        
        _udl_list: List[Tuple] = []
        for _index in _index_list:
            _udl_list.append(LineBeam(self._qx1[_index], self._qy1[_index], self._qz1[_index],
                                      self._qx2[_index], self._qy2[_index], self._qz2[_index],
                                      self._L1[_index], self._L2[_index],
                                      self._labels[_index], self._title[_index],
                                      self._system[_index], self._complex[_index]))
        return _udl_list
    #
    @property
    def name(self) -> str:
        """
        """
        return self._title[self._index]
    
    @name.setter
    def name(self, load_name:str) -> None:
        """
        """
        try:
            self._title[self._index] = load_name
        except AttributeError:
            #self.load_name = load_name
            raise IndexError("load name not found")    
    #
    #
    @property
    def coordinate_system(self):
        if self._system_flag != 0:
            return "local"
        return "global"
    
    @coordinate_system.setter
    def coordinate_system(self, system:Union[str,int]):
        """
        Coordinate system for load : global or local (member)
        """
        self._system_flag = 0
        if system in ['local', 'member', 1]:
            self._system_flag = 1
        
    #
    #
    def __len__(self) -> float:
        return len(self._labels)

    def __delitem__(self, element_name: int) -> None:
        """
        """
        _index_list: List = [x for x, _item in enumerate(self._labels)
                             if _item == element_name]

        _index_list.sort(reverse=True)
        for _index in _index_list:
            self._L1.pop(_index)
            self._qx1.pop(_index)
            self._qy1.pop(_index)
            self._qz1.pop(_index)
            #
            self._L2.pop(_index)
            self._qx2.pop(_index)
            self._qy2.pop(_index)
            self._qz2.pop(_index)
            #
            self._labels.pop(_index)
            self._title.pop(_index)
            self._complex.pop(_index)

    def __contains__(self, value) -> bool:
        return value in self._labels

    def __iter__(self) -> Iterable:
        """
        """
        items = list(set(self._labels))
        return iter(items)
    #
    #
    def get_nodal_load(self, elements, nodes, materials, sections) -> List:
        """
        """
        items = [ ]
        set_items = set(self._labels)
        if not set_items:
            return []
        #
        for item in set_items:
            _index = self._labels.index(item)
            element = elements[item]
            length = element.length(nodes)
            material = materials[element.material]
            emod = material.E 
            gmod = material.G 
            section = sections[element.section]
            area = section.area 
            iy = section.Iz 
            iz = section.Iy
            #
            #elementR = element.R
            #unit_vector = element.unit_vector
            #
            nloads = self.__getitem__(item)
            for nl in nloads:
                # local nodal loading
                nload = [nl[0], nl[1], nl[2], 0, 0, 0,
                         nl[3], nl[4], nl[5], 0, 0, 0,]
                L1 = nl.L1
                L2 = nl.L2
                #
                try: # local system
                    1/self._system[_index]
                    #print('local system')                   
                except ZeroDivisionError:
                    # global system
                    #print('global system')
                    univec = element.unit_vector(nodes)
                    R = Rmatrix(*univec, element.beta)
                    nload = trns_3Dv(nload, R)                  
                #
                lnload = [0]*12
                #
                if nload[0] or nload[6]:
                    h = length - L1 - L2
                    axial, xc = self.trapezoid(a=-nload[0], b=-nload[6], h=h) # N/m neg
                    lnload[0], lnload[6] = self.distributed_axial(axial, length,
                                                                  l1=L1, l2=L2)                
                #
                # This is in local system
                if nload[1] or nload[7]:
                    try:
                        1/(abs(L1) + abs(L2))
                        eq_nload = self.distributed_load3(w1=nload[1],
                                                          w2=nload[7],
                                                          L=length,
                                                          l1=L1, l2=L2)
                    except ZeroDivisionError:
                        eq_nload = self.distributed_load(pa=nload[1],
                                                         pb=nload[7],
                                                         L=length,
                                                         EI=emod*iy,
                                                         GAs=gmod*area)
                    # node 1
                    lnload[1] += eq_nload[0]
                    lnload[5] -= eq_nload[1]
                    # node 2
                    lnload[7]  += eq_nload[2]
                    lnload[11] -= eq_nload[3]
                #node_load_z 
                if nload[2] or nload[8]:
                    try:
                        1/(abs(L1) + abs(L2))
                        eq_nload = self.distributed_load3(w1=nload[2],
                                                          w2=nload[8],
                                                          L=length,
                                                          l1=L1, l2=L2)
                    except ZeroDivisionError:
                        eq_nload = self.distributed_load(pa=nload[2],
                                                         pb=nload[8],
                                                         L=length,
                                                         EI=emod*iz,
                                                         GAs=gmod*area)
                    # node 1
                    lnload[2] += eq_nload[0]
                    lnload[4] += eq_nload[1]
                    # node 2
                    lnload[8]  += eq_nload[2]
                    lnload[10] += eq_nload[3]
                #
                items.append([*lnload, element.name])
        return  items
    #
    def distributed_load(self, pa:float, pb:float, 
                         L:float, EI:float, GAs:float) -> List[float]:
        """
        pa:
        pb:
        L:
        EI:
        GAs:
        """
        Pa = pa * (L**4 / (30 * EI)  - L**2 / (3 * GAs))
        Pb = pb * (L**4 / (120 * EI) - L**2 / (6 * GAs))
        Fw = Pa + Pb
        Ftheta = -(3*pa + pb) * L**3/(24*EI)
        Fv = -(pa+pb) * L/2.0
        Fm = -(2*pa+pb) * L**2 / 6.0
        nodal_load = [Fw, Ftheta, Fv, Fm]
        return nodal_load_vector(nodal_load, L, EI, GAs)
    #
    def distributed_load3(self, w1:float, w2:float, 
                          L:float, l1:float, l2:float):
        """
        Case 4 from Matrix Analysis of Framed Structures [Aslam Kassimali]
        """
        Fa = (w1*(L - l1)**3/(20*L**3) * ((7*L + 8*l1) - l2*(3*L + 2*l1)/(L - l1)
                                          * (1 + l2/(L-l1) + l2**2/(L-l1)**2) 
                                          + 2*l2**4/(L-l1)**3)
              + w2*(L-l1)**3/(20*L**3) * ((3*L + 2*l1) * (1 + l2/(L-l1) + l2**2/(L-l1)**2)
                                          - l2**3/(L-l1)**2 * (2 + (15*L - 8*l2)/(L-l1))))
        #
        Ma = (w1*(L - l1)**3/(60*L**2) * (3*(L + 4*l1) - l2*(2*L + 3*l1)/(L - l1)
                                          * (1 + l2/(L-l1) + l2**2/(L-l1)**2)
                                          + 3*l2**4/(L-l1)**3)
              + w2*(L-l1)**3/(60*L**2) * ((2*L + 3*l1) * (1 + l2/(L-l1) + l2**2/(L-l1)**2)
                                          - 3*l2**3/(L-l1)**2 * (1 + (5*L - 4*l2)/(L-l1))))
        #
        Fb = (w1+w2)/2.0 * (L-l1-l2) - Fa
        Mb = ((L-l1-l2)/6.0 * (w1*(-2*L + 2*l1 - l2) 
                               - w2*(L-l1+2*l2)) + Fa*L - Ma)
        return [Fa, Ma, Fb, Mb]
    #
    def distributed_load2(self, w, L, l1=0, l2=0):
        """
        """
        Fa = w*L/2*(1 - l1/L**4*(2*L**3 - 2*l1**2*L + l1**3) - l2**3/L**4*(2*L-l2))
        Ma = w*L**2/12 * (1 - l1**2/L**4 *(6*L**2 - 8*l1*L + 3*l1**2) - l2**3/L**4*(4*L - 3*l2))
        Fb = w*L/2*(1 - l1**3/L**4*(6*L - l1) - l2/L**4*(2*L**3 - 2*l2**2 + l2**3))
        Mb = -w*L**2/12 * (1 - l1**3/L**4 *(4*L - 3*l1) - l2**2/L**4*(6*L**2 - 8*l2*L + 3*l2**2))
        return [Fa, Ma, Fb, Mb]
    #
    def distributed_moment(self, ca:float, cb:float, 
                           L:float, EI:float, GAs:float) -> List[float]:
        """
        ca:
        cb:
        L:
        EI:
        GAs:
        """
        Fw = (3*ca + cb) * L**3 / (24*EI)
        Ftheta = -(2*ca + cb) * L**2 / (6*EI)
        Fv = 0
        Fm = -(ca + cb) * L/2.0
        nodal_load = [Fw, Ftheta, Fv, Fm]
        return nodal_load_vector(nodal_load, L, EI, GAs)        
    #
    def distributed_axial(self, w:float, L:float, 
                          l1:float, l2:float) -> List[float]:
        """ 
        Case 6 from Matrix Analysis of Framed Structures [Aslam Kassimali]
        """
        Fa = w/(2*L) * (L-l1-l2) * (L-l1+l2)
        Fb = w/(2*L) * (L-l1-l2) * (L+l1-l2)
        return [Fa, Fb]
    #
    def trapezoid(self, a, b, h):
        """
        """
        area = h * (a+b)/2.0
        ceq = area / h
        
        if a < b:
            xb = h/3.0 * (2*a + b) / (a+b)
            xb = h - xb
        else:
            xb = h/3.0 * (2*b + a) / (a+b)
        return [ceq, xb]
    #
    def _get_line_load(self, load):
        """ """
        if isinstance(load, (list, tuple)):
            try:
                load = check_list_units(load)
            except AttributeError:
                load = check_list_number(load, steps=8)
        elif isinstance(load, dict):
            load = check_beam_dic(load)
        else:
            raise Exception('   *** Load input format not recognized')
        return load
#
#
#
class BeamPoint(NodeLoadMaster):
    __slots__ =  ['_title', '_labels', '_index', '_complex',
                  '_fx', '_fy', '_fz', '_mx', '_my', '_mz']
    
    def __init__(self) -> None:
        """
        """
        super().__init__()
    #
    #
    def __setitem__(self, element_name:Union[int, str], 
                    point_load: Union[List[float], Dict[str,float]]) -> None:
        """
        """
        self._labels.append(element_name)
        self._title.append(element_name)
        self._system.append(self._system_flag)
        self._complex.append(0)
        #
        point_load = self._get_beam_load(point_load)
        self._fx.append(point_load[0])
        self._fy.append(point_load[1])
        self._fz.append(point_load[2])
        self._mx.append(point_load[3])
        self._my.append(point_load[4])
        self._mz.append(point_load[5])
        self._distance.append(point_load[6])
        #
        self._index = len(self._labels) - 1
        #print('--')
    
    def __getitem__(self, element_name:Union[int, str])-> List[Tuple]:
        """
        """
        _index_list: List = [x for x, _item in enumerate(self._labels)
                             if _item == element_name]
        #
        _points: List = []
        for _index in _index_list:
            _points.append(PointBeam(self._fx[_index], self._fy[_index], self._fz[_index],
                                     self._mx[_index], self._my[_index], self._mz[_index],
                                     self._distance[_index], self._labels[_index], 
                                     self._title[_index], self._system[_index], 
                                     self._complex[_index]))
        return _points    
    #
    @property
    def coordinate_system(self):
        if self._system_flag != 0:
            return "local"
        return "global"
    
    @coordinate_system.setter
    def coordinate_system(self, system:Union[str,int]):
        """
        Coordinate system for load : global or local (member)
        """
        self._system_flag = 0
        if system in ['local', 'member', 1]:
            self._system_flag = 1
        
    #    
    #
    def get_nodal_load(self, elements, nodes, materials, sections) -> List:
        """
        """
        items = [ ]
        set_items = set(self._labels)
        if not set_items:
            return []
        for item in set_items:
            _index = self._labels.index(item)
            element = elements[item]
            length = element.length(nodes)
            material = materials[element.material]
            emod = material.E 
            gmod = material.G 
            section = sections[element.section]
            area = section.area 
            iy = section.Iz 
            iz = section.Iy 
            #
            # local nodal loading
            ploads = self.__getitem__(item)
            for pload in ploads:
                # local nodal loading
                Pdistance = pload.distance
                lnload = [0]*12
                try: # local system
                    1/self._system[_index]
                    pload = list(pload[:6])
                    #print('local system')
                except ZeroDivisionError:
                    # global system
                    #print('global system')
                    univec = element.unit_vector(nodes)
                    R = Rmatrix(*univec, element.beta)
                    pload = trns_3Dv(pload[:6], R)
                #
                # This is in local system
                # X axial
                if pload[0]:
                    lnload[0], lnload[6] = self.axial_load(-pload[0], 
                                                           length,
                                                           Pdistance)
                # Y
                if pload[1]:
                    eq_nload = self.point_load(W=pload[1], # negative
                                               L=length,
                                               lp=Pdistance, 
                                               EI=emod*iy, 
                                               GAs=gmod*area)                  
                    # node 1
                    lnload[1] += eq_nload[0]
                    lnload[5] -= eq_nload[1]
                    # node 2
                    lnload[7]  += eq_nload[2]
                    lnload[11] -= eq_nload[3]
                # Z
                if pload[2]:
                    eq_nload = self.point_load(W=pload[2], 
                                               L=length,
                                               lp=Pdistance, 
                                               EI=emod*iz, 
                                               GAs=gmod*area)
                    #                    
                    # node 1
                    lnload[2] += eq_nload[0]
                    lnload[4] += eq_nload[1]
                    # node 2
                    lnload[8]  += eq_nload[2]
                    lnload[10] += eq_nload[3]
                # Mx torsion
                if pload[3]:
                    lnload[3], lnload[9] = self.torsion_load(pload[3], 
                                                             length, 
                                                             Pdistance)
                # My
                if pload[4]:
                    eq_nload = self.point_moment(C=pload[4], # negative?
                                                 L=length,
                                                 lp=Pdistance, 
                                                 EI=emod*iz, 
                                                 GAs=gmod*area)
                    #
                    #eq_nload = self.point_moment2(M=pload[4], 
                    #                             L=length,
                    #                             l1=Pdistance)
                    # node 1
                    lnload[2] += eq_nload[0]
                    lnload[4] += eq_nload[1]
                    # node 2
                    lnload[8]  += eq_nload[2]
                    lnload[10] += eq_nload[3]
                # Mz
                if pload[5]:
                    eq_nload = self.point_moment(C=-pload[5], 
                                                 L=length,
                                                 lp=Pdistance, 
                                                 EI=emod*iy, 
                                                 GAs=gmod*area)
                    #
                    #eq_nload = self.point_moment2(M=pload[5], 
                    #                              L=length,
                    #                              l1=Pdistance)
                    # node 1
                    lnload[1] += eq_nload[0]
                    lnload[5] -= eq_nload[1]
                    # node 2
                    lnload[7]  += eq_nload[2]
                    lnload[11] -= eq_nload[3]
            #
            items.append([*lnload, element.name])
        return  items    
    #
    def point_load(self, W:float, L:float, lp:float, 
                   EI:float, GAs:float) -> List[float]:
        """
        """
        l1 = L - lp
        Fw = W*(l1**3/(6*EI) - l1/GAs)
        Ftheta = -W*l1**2 / (2*EI)
        Fv = -W
        Fm = -W*l1
        nodal_load = [Fw, Ftheta, Fv, Fm]
        return nodal_load_vector(nodal_load, L, EI, GAs)
    #
    def point_moment(self, C:float, L:float, lp:float, 
                     EI:float, GAs:float) -> List[float]:
        """
        """
        l1 = L - lp
        Fw = C * l1**2 / (2*EI)
        Ftheta = -C * l1 / EI
        Fv = 0
        Fm = -C
        nodal_load = [Fw, Ftheta, Fv, Fm]
        return nodal_load_vector(nodal_load, L, EI, GAs)        
    #
    def point_moment2(self, M:float, L:float, l1:float) -> List[float]:
        """
        Case 2 from Matrix Analysis of Framed Structures [Aslam Kassimali]
        """
        l2 = L - l1
        Fa = -6 * M * l1*l2 / L**3
        Ma = M*l2/L**2 * (l2-2*l1)
        Fb = 6 * M * l1*l2 / L**3
        Mb = M*l1/L**2 * (l1-2*l2)
        return [Fa, Ma, Fb, Mb]
        
    #
    def axial_load(self, W:float, L:float, l1:float) -> List[float]:
        """
        Case 5 from Matrix Analysis of Framed Structures [Aslam Kassimali]
        """
        l2 = L - l1
        Fa = W*l2/L
        Fb = W*l1/L
        return [Fa, Fb]
    #
    def torsion_load(self, Mt:float, L:float, l1:float) -> List[float]:
        """
        Case 7 from Matrix Analysis of Framed Structures [Aslam Kassimali]
        """
        l2 = L - l1
        Fa = Mt*l2/L
        Fb = Mt*l1/L
        return [Fa, Fb]
    #
    def point_load2(self, W:float, L:float, l1:float) -> List[float]:
        """
        Case 1 from Matrix Analysis of Framed Structures [Aslam Kassimali]
        """
        l2 = L - l1
        Fa = W*l2**2/L**3 / (3*l1 + l2)
        Ma = W*l1*l2**2/L**2
        Fb = W*l1**2/L**3 / (l1 + 3*l2)
        Mb = -W*l1**2 * l2/L**2
        return [Fa, Ma, Fb, Mb]
#
#
#@staticmethod
def nodal_load_vector(nodal_load:List[float], L:float, 
                      EI:float, GAs:float) -> List[float]:
    """
    return [point load end 1, moment end 1, point load end 2, moment end 2]
    """
    Fw, Ftheta, fv, fm = nodal_load
    delta = 2*L * (L**2/(12*EI) + 1/GAs)
    Va = -(2*Fw + L * Ftheta)/delta
    Ma = (Fw*L/EI + 2* Ftheta*(L**2/(6*EI) - 1/GAs)) * EI / delta
    Vb = -fv + (2*Fw + L * Ftheta)/delta
    Mb = -fm + (Fw*L/EI + 2* Ftheta *(L**2/(3*EI) + 1/GAs))* EI / delta
    return [ Va, -Ma, Vb, -Mb ]
#
