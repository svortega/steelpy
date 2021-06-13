#
# Copyright (c) 2009-2021 fem2ufo
# 

# Python stdlib imports
from array import array
from collections.abc import Mapping
#from collections import defaultdict
#from dataclasses import dataclass
from typing import NamedTuple, Tuple, List, Union, Iterable, Dict  


# package imports
from steelpy.trave3D.processor.operations import trns_3Dv
from steelpy.trave3D.preprocessor.assemble import Rmatrix
from steelpy.f2uModel.load.operations.operations import(check_list_units,
                                                        check_list_number,
                                                        check_beam_dic)
from steelpy.process.math.vector import Vector
#
# ---------------------------------
#
class LineBeam(NamedTuple):
    """
    """
    # FIXME: start 1 intead zero
    qx0: float
    qy0: float
    qz0: float
    #
    qx1: float
    qy1: float
    qz1: float
    #
    L0: float
    L1: float
    #
    number: int
    load_name: str
    system:int
    load_complex:int
    #
    @property
    def coordinate_system(self):
        if self.system != 0:
            return "local"
        return "global"
    #
    def node_equivalent(self, elements, materials, sections):
        """ """
        #nodal_load = []
        eq_ln = self.line2node(elements, materials, sections)
        global_nodal_load = local2global(eq_ln, elements)
        local_nodal_load = beam_eq(eq_ln)
        return global_nodal_load, local_nodal_load
    #
    def line2node(self, element, material, section) -> List:
        """
        """
        length = element.length
        emod = material.E.convert('pascal').value
        gmod = material.G.convert('pascal').value
        area = section.area
        iy = section.Iz
        iz = section.Iy
        #
        # local nodal loading
        nload = [self.qx0, self.qy0, self.qz0, 0, 0, 0,
                 self.qx1, self.qy1, self.qz1, 0, 0, 0,]
        L1 = self.L0
        L2 = self.L1
        #
        try:  # local system
            1 / self.system
            # print('local system')
        except ZeroDivisionError:
            # global system
            # print('global system')
            univec = element.unit_vector
            R = Rmatrix ( *univec, element.beta )
            nload = trns_3Dv ( nload, R )
            #
        lnload = [0] * 12
        #
        if nload[0] or nload[6]:
            h = length - L1 - L2
            axial, xc = trapezoid(a=-nload[0], b=-nload[6], h=h)  # N/m neg
            lnload[0], lnload[6] = distributed_axial(axial, length,
                                                     l1=L1, l2=L2)
        #
        # This is in local system
        if nload[1] or nload[7]:
            try:
                1 / (abs(L1) + abs(L2))
                eq_nload = distributed_load3(w1=nload[1],
                                             w2=nload[7],
                                             L=length,
                                             l1=L1, l2=L2)
            except ZeroDivisionError:
                eq_nload = distributed_load(pa=nload[1],
                                            pb=nload[7],
                                            L=length,
                                            EI=emod * iy,
                                            GAs=gmod * area)
            # node 1
            lnload[1] += eq_nload[0]
            lnload[5] -= eq_nload[1]
            # node 2
            lnload[7]  += eq_nload[2]
            lnload[11] -= eq_nload[3]
        # node_load_z
        if nload[2] or nload[8]:
            try:
                1 / (abs(L1) + abs(L2))
                eq_nload = distributed_load3(w1=nload[2],
                                             w2=nload[8],
                                             L=length,
                                             l1=L1, l2=L2)
            except ZeroDivisionError:
                eq_nload = distributed_load(pa=nload[2],
                                            pb=nload[8],
                                            L=length,
                                            EI=emod * iz,
                                            GAs=gmod * area)
            # node 1
            lnload[2] += eq_nload[0]
            lnload[4] += eq_nload[1]
            # node 2
            lnload[8]  += eq_nload[2]
            lnload[10] += eq_nload[3]
        #
        return lnload
    #
    def __str__(self, units:str="si") -> str:
        """ """
        output  = (f"{str(self.number):12s} {self.L0: 1.3e} {self.qx0: 1.3e} "
                   f"{self.qy0: 1.3e} {self.qz0: 1.3e} "
                   f"{self.coordinate_system.upper()}\n")
        if (step := self.load_name) == "NULL":
            step = 12 * " "
        output += (f"{step[:12]:12s} {self.L1: 1.3e} {self.qx1: 1.3e} "
                   f"{self.qy1: 1.3e} {self.qz1: 1.3e} "
                   f"{self.load_complex}\n")
        return output
#
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
    L0:float
    number: int
    load_name: str
    system:int
    load_complex:int
    #
    @property
    def distance(self):
        """ """
        return self.L0
    #
    def node_equivalent(self, elements, materials, sections):
        """ """
        #nodal_load = []
        eq_ln = self.point2node(elements, materials, sections)
        global_nodal_load = local2global(eq_ln, elements)
        local_nodal_load = beam_eq(eq_ln)
        return global_nodal_load, local_nodal_load
    #
    def point2node(self, element, material, section) -> List:
        """
        """
        length = element.length
        emod = material.E.convert('pascal').value
        gmod = material.G.convert('pascal').value
        area = section.area
        iy = section.Iz
        iz = section.Iy
        #
        # local nodal loading
        #ploads = self.__getitem__ ( item )
        pload = [self.fx, self.fy, self.fz, self.mx, self.my, self.mz]
        #for pload in ploads:
        # local nodal loading
        Pdistance = self.L0
        lnload = [0] * 12
        try:  # local system
            1 / self.system
            #pload = list( pload[ :6 ] )
            # print('local system')
        except ZeroDivisionError:
            # global system
            # print('global system')
            univec = element.unit_vector
            R = Rmatrix( *univec, element.beta )
            pload = trns_3Dv (pload, R)
        #
        # This is in local system
        # X axial
        if pload[0]:
            lnload[0], lnload[6] = axial_load(-pload[0],
                                              length,
                                              Pdistance)
        # Y
        if pload[1]:
            eq_nload = point_load(W=pload[1],  # negative
                                  L=length,
                                  lp=Pdistance,
                                  EI=emod * iy,
                                  GAs=gmod * area)
            # node 1
            lnload[1] += eq_nload[0]
            lnload[5] -= eq_nload[1]
            # node 2
            lnload[7]  += eq_nload[2]
            lnload[11] -= eq_nload[3]
        # Z
        if pload[2]:
            eq_nload = point_load(W=pload[2],
                                  L=length,
                                  lp=Pdistance,
                                  EI=emod * iz,
                                  GAs=gmod * area)
            #
            # node 1
            lnload[2] += eq_nload[0]
            lnload[4] += eq_nload[1]
            # node 2
            lnload[8]  += eq_nload[2]
            lnload[10] += eq_nload[3]
        # Mx torsion
        if pload[3]:
            lnload[3], lnload[9] = torsion_load(pload[3],
                                                length,
                                                Pdistance)
        # My
        if pload[4]:
            eq_nload = point_moment(C=pload[4],  # negative?
                                    L=length,
                                    lp=Pdistance,
                                    EI=emod * iz,
                                    GAs=gmod * area)
            #
            # eq_nload = self.point_moment2(M=pload[4],
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
            eq_nload = point_moment(C=-pload[5],
                                    L=length,
                                    lp=Pdistance,
                                    EI=emod * iy,
                                    GAs=gmod * area)
            #
            # eq_nload = self.point_moment2(M=pload[5],
            #                              L=length,
            #                              l1=Pdistance)
            # node 1
            lnload[1] += eq_nload[0]
            lnload[5] -= eq_nload[1]
            # node 2
            lnload[7]  += eq_nload[2]
            lnload[11] -= eq_nload[3]
        #
        return lnload
    #
    @property
    def coordinate_system(self):
        if self.system != 0:
            return "local"
        return "global"    
    #
    def __str__(self, units:str="si") -> str:
        """ """
        output  = (f"{str(self.number):12s} {self.L0: 1.3e} {self.fx: 1.3e} "
                   f"{self.fy: 1.3e} {self.fy: 1.3e} "
                   f"{self.coordinate_system.upper()}\n")
        if (step := self.load_name) == "NULL":
            step = 12 * " "
        output += (f"{step[:12]:12s} {10*' '} {self.mx: 1.3e} "
                   f"{self.my: 1.3e} {self.mz: 1.3e} "
                   f"{self.load_complex}\n")
        return output    
#
#
class BeamDistMaster(Mapping):
    
    def __init__(self) -> None:
        """
        """
        self._labels: List[Union[str, int]] = []
        self._title: List[str] = []
        self._index: int
        self._complex: array = array("I", [])
        # 0-global/ 1-local
        self._system_flag:int = 0
        self._system: array = array("I", [])
    #
    def __len__(self) -> float:
        return len(self._labels)
    #
    def __contains__(self, value) -> bool:
        return value in self._labels
    #
    def __iter__(self) -> Iterable:
        """
        """
        items = list(set(self._labels))
        return iter(items)
    #
    def __str__(self) -> str:
        """ """
        print('---')    
#
# ---------------------------------
#
def local2global(eq_lnloads, element):
    """ """
    #nodal_load = [ ]
    #item = eq_lnloads
    # global nodal load
    univec = element.unit_vector
    R = Rmatrix(*univec, element.beta)
    gnload = trns_3Dv(eq_lnloads[:12], R)
    #
    # axial
    gnload[0] *= -1
    gnload[6] *= -1
    #
    # rotation X
    #gnload[3] *= -1
    #gnload[9] *= -1
    #
    # rotation Y
    gnload[4] *= -1
    gnload[10] *= -1
    #
    # rotation Z
    gnload[5] *= -1
    gnload[11] *= -1
    #
    # node 1
    #n_index = nodes[end_nodes[0]].index
    #nodal_load[n_index] += Vector(gnload[:6])
    #nodal_load.append(Vector(gnload[:6]))
    # node 2
    #n_index = nodes[end_nodes[1]].index
    #nodal_load[n_index] += Vector(gnload[6:])
    #nodal_load.append(Vector(gnload[6:]))
    #return nodal_load
    return Vector(gnload)
#
def beam_eq(eq_lnloads):
    """ """
    #member_nload = [] # zeros(nmb, 12)
    #for item in eq_lnloads:
    item = eq_lnloads
    #
    # axial
    item[0] *= -1
    item[6] *= -1
    #
    # rotation X
    #item[3] *= -1
    #item[9] *= -1
    #
    # rotation Y
    item[4] *= -1
    item[10] *= -1
    #
    item[5] *= -1
    item[11] *= -1
    #
    #try:
        #member_nload += Vector(item[:12])
    #except KeyError:
        #member_nload = Vector(item[:12])
    #return member_nload
    return Vector(item[:12])
#
# ---------------------------------
#
#
def line2nodeY(self, elements, materials, sections) -> List:
    """
    """
    #items = [ ]
    set_items = set(self._labels)
    if not set_items:
        return []
    #
    for item in set_items:
        #_index = self._labels.index(item)
        element = elements[item]
        material = materials[element.material]
        section = sections[element.section].properties
        #property = section.properties
        #
        #nodes = element.nodes
        #
        nloads = self.__getitem__(item)
        for nl in nloads:
            nodal_load, member_nload = nl.node_equivalent(element, material, section)
            print('-->')
#
def line2nodeX(self, elements, materials, sections) -> List:
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
        length = element.length
        material = materials[element.material]
        emod = material.E.convert('pascal').value
        gmod = material.G.convert('pascal').value 
        section = sections[element.section].properties
        #property = section.properties
        area = section.area
        iy = section.Iz
        iz = section.Iy
        #
        #nodes = element.nodes
        #
        nloads = self.__getitem__(item)
        for nl in nloads:
            # FIXME:
            #items = nl.node_equivalent(element, material, section)
            # local nodal loading
            nload = [nl[0], nl[1], nl[2], 0, 0, 0,
                     nl[3], nl[4], nl[5], 0, 0, 0,]
            L1 = nl.L0
            L2 = nl.L1
            #
            try: # local system
                1/self._system[_index]
                #print('local system')                   
            except ZeroDivisionError:
                # global system
                #print('global system')
                univec = element.unit_vector
                R = Rmatrix(*univec, element.beta)
                nload = trns_3Dv(nload, R)                  
            #
            lnload = [0]*12
            #
            if nload[0] or nload[6]:
                h = length - L1 - L2
                axial, xc = trapezoid(a=-nload[0], b=-nload[6], h=h) # N/m neg
                lnload[0], lnload[6] = distributed_axial(axial, length,
                                                         l1=L1, l2=L2)
            #
            # This is in local system
            if nload[1] or nload[7]:
                try:
                    1/(abs(L1) + abs(L2))
                    eq_nload = distributed_load3(w1=nload[1],
                                                 w2=nload[7],
                                                 L=length,
                                                 l1=L1, l2=L2)
                except ZeroDivisionError:
                    eq_nload = distributed_load(pa=nload[1],
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
                    eq_nload = distributed_load3(w1=nload[2],
                                                 w2=nload[8],
                                                 L=length,
                                                 l1=L1, l2=L2)
                except ZeroDivisionError:
                    eq_nload = distributed_load(pa=nload[2],
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
#
def distributed_load(pa:float, pb:float,
                     L:float, EI:float, GAs:float) -> List[float]:
    """
    pa:
    pb:
    L:
    EI:
    GAs:
    """
    Pa = pa * (L**4 / (30 * EI) - L**2 / (3 * GAs))
    Pb = pb * (L**4 / (120 * EI) - L**2 / (6 * GAs))
    Fw = Pa + Pb
    Ftheta = -(3*pa + pb) * L**3/(24*EI)
    Fv = -(pa+pb) * L/2.0
    Fm = -(2*pa+pb) * L**2 / 6.0
    nodal_load = [Fw, Ftheta, Fv, Fm]
    return nodal_load_vector(nodal_load, L, EI, GAs)
#
def distributed_load3(w1:float, w2:float,
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
def distributed_load2(w, L, l1=0, l2=0):
    """
    """
    Fa = w*L/2*(1 - l1/L**4*(2*L**3 - 2*l1**2*L + l1**3) - l2**3/L**4*(2*L-l2))
    Ma = w*L**2/12 * (1 - l1**2/L**4 *(6*L**2 - 8*l1*L + 3*l1**2) - l2**3/L**4*(4*L - 3*l2))
    Fb = w*L/2*(1 - l1**3/L**4*(6*L - l1) - l2/L**4*(2*L**3 - 2*l2**2 + l2**3))
    Mb = -w*L**2/12 * (1 - l1**3/L**4 *(4*L - 3*l1) - l2**2/L**4*(6*L**2 - 8*l2*L + 3*l2**2))
    return [Fa, Ma, Fb, Mb]
#
def distributed_moment(ca:float, cb:float,
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
def distributed_axial(w:float, L:float,
                      l1:float, l2:float) -> List[float]:
    """
    Case 6 from Matrix Analysis of Framed Structures [Aslam Kassimali]
    """
    Fa = w/(2*L) * (L-l1-l2) * (L-l1+l2)
    Fb = w/(2*L) * (L-l1-l2) * (L+l1-l2)
    return [Fa, Fb]
#
def trapezoid(a, b, h):
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
# ---------------------------------
#
def point2node(self, elements, materials, sections) -> List:
    """
    """
    items = [ ]
    set_items = set(self._labels)
    if not set_items:
        return []
    for item in set_items:
        ploads = self.__getitem__(item)
        if not ploads:
            continue
        #
        _index = self._labels.index(item)
        element = elements[item]
        length = element.length
        material = materials[element.material]
        emod = material.E.convert('pascal').value
        gmod = material.G.convert('pascal').value
        section = sections[element.section].properties
        area = section.area 
        iy = section.Iz 
        iz = section.Iy 
        #
        # local nodal loading
        ploads = self.__getitem__(item)
        for pload in ploads:
            # local nodal loading
            Pdistance = pload.L0
            lnload = [0]*12
            try: # local system
                1/self._system[_index]
                pload = list(pload[:6])
                #print('local system')
            except ZeroDivisionError:
                # global system
                #print('global system')
                univec = element.unit_vector
                R = Rmatrix(*univec, element.beta)
                pload = trns_3Dv(pload[:6], R)
            #
            # This is in local system
            # X axial
            if pload[0]:
                lnload[0], lnload[6] = axial_load(-pload[0], 
                                                  length,
                                                  Pdistance)
            # Y
            if pload[1]:
                eq_nload = point_load(W=pload[1], # negative
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
                eq_nload = point_load(W=pload[2], 
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
                lnload[3], lnload[9] = torsion_load(pload[3], 
                                                    length, 
                                                    Pdistance)
            # My
            if pload[4]:
                eq_nload = point_moment(C=pload[4], # negative?
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
                eq_nload = point_moment(C=-pload[5], 
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
#
def point_load(W:float, L:float, lp:float, 
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
def point_moment(C:float, L:float, lp:float, 
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
def point_moment2(M:float, L:float, l1:float) -> List[float]:
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
#
def axial_load(W:float, L:float, l1:float) -> List[float]:
    """
    Case 5 from Matrix Analysis of Framed Structures [Aslam Kassimali]
    """
    l2 = L - l1
    Fa = W*l2/L
    Fb = W*l1/L
    return [Fa, Fb]
#
def torsion_load(Mt:float, L:float, l1:float) -> List[float]:
    """
    Case 7 from Matrix Analysis of Framed Structures [Aslam Kassimali]
    """
    l2 = L - l1
    Fa = Mt*l2/L
    Fb = Mt*l1/L
    return [Fa, Fb]
#
def point_load2(W:float, L:float, l1:float) -> List[float]:
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
# ---------------------------------
#
def get_line_loadXX(load):
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
    #
    return load
#
#
def get_beam_point_load(data:List[float], steps:int=6)->List[float]:
    """ """
    new_data = [ ]
    # point load
    for x in range(1, steps+1):
        try:
            new_data.append(data[x])
        except IndexError:
            new_data.append(0.0)
    # L
    new_data.append(data[0])
    return new_data
#
#
