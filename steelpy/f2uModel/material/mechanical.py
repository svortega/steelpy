# 
# Copyright (c) 2009-2019 fem2ufo
# 
"""
defines mechanical material type: linear or non-linear material
"""

# Python stdlib imports
from array import array
from dataclasses import dataclass
from typing import NamedTuple, ClassVar, Dict, List, Tuple, Union, Iterable

# package imports
from steelpy.process.units.main import Units
import steelpy.f2uModel.material.print_report as print_material
 
#
@dataclass
class GetMaterial:
    """ Linear Material Data"""

    __slots__ = ['name', 'index', 'cls', 'type']

    def __init__(self, cls, material_name) -> None:
        """
        """
        self.index: int = cls._labels.index(material_name)
        self.cls: ClassVar = cls
        self.name: Union[int, str] = material_name
        self.type = "elastic"
    #
    @property
    def Fy(self):
        return self.cls._Fy[self.index] * f2u_units.Pa
    #
    @Fy.setter
    def Fy(self, item):
        self.cls._Fy[self.index] = item.value
    #
    @property
    def E(self):
        return self.cls._E[self.index] * f2u_units.Pa
    
    @E.setter
    def E(self, item):
        self.cls._E[self.index] = item.value
    #
    @property
    def G(self):
        return self.cls._G[self.index] * f2u_units.Pa
    
    @G.setter
    def G(self, item):
        self.cls._G[self.index] = item.value     
    #
    @property
    def Fu(self):
        try:
            1/self.cls._Fu[self.index]
            return self.cls._Fu[self.index] * f2u_units.Pa
        except ZeroDivisionError:
            return self.cls._Fy[self.index] / 0.75 * f2u_units.Pa
    
    @Fu.setter
    def Fu(self, item):
        """
        Fu : 
        """
        self.cls._Fu[self.index] = item.value
    #
    @property
    def density(self):
        return self.cls._density[self.index] * f2u_units.kg/f2u_units.m**3
    
    @density.setter
    def density(self, item):
        self.cls._density[self.index] = item.value
    #
    #
    @property
    def alpha(self):
        """"""
        return self.cls._alpha[self.index] * f2u_units.K
    #
    @property
    def poisson(self):
        """"""
        return self.cls._poisson[self.index]
    #
    def get_name(self, number, _grade=False):
        """
        """
        if not _grade:
            try:
                _grade = round(self.cls._grade[self.index])
            except TypeError:
                _grade = 'Zero'
        # in case material already named (i.e. sacs case)
        _name = '_' + self.name
        if '_MT' in _name:
            _name = ''
        #
        self.name[self.index] = 'MT' + str(number).zfill(3) + '_' + str(_grade) + _name
    #
    def equal(self, other, grade=None):
        """
        """
        if not grade:
            grade = other.grade
        # TODO : check if material type is needed
        if self.cls._type[self.index] == other.type \
        and self.cls._E[self.index] == other.E \
        and self.cls._grade[self.index] == grade \
        and self.cls._poisson[self.index] == other.poisson \
        and self.cls._density[self.index] == other.density :
            return True
        else:
            return False
    #
    def print_properties(self):
        """ """
        return print_material.print_plastic_material(self)
        #for line in text:
        #    print(line.rstrip())
    #
    def set_default(self):
        """ """
        self.cls._cls._default = self.name
#
#
#
class MaterialElastic:
    """
    Represents an istotrop linear elastic material used for FEM simulations
    
    Material
        |_ name
        |_ number
        |_ type
        |_ Emodulus
        |_ poisson
        |_ Fy
        |_ density
        |_ alpha
        |_ damping
        |_ stiffness [k1, k2, k3,..., kn]
        |_ spring    [sp[0],...]
        |_ Gmodulus (shear modulus)
    
    **Parameters**:  
      :number:  integer internal number 
      :name:  string node external name
      :coordinate:  list node coordinates 
    """
    __slots__ = ['_E', '_poisson', '_density', '_cls',
                 #'units', 'type', 'sets',
                 '_Fy', '_Fu', '_damping', '_alpha', 'f2u_units',
                  '_grade', '_G', '_labels', '_number']
    
    def __init__(self, cls)-> None:
        """
        """
        global f2u_units
        f2u_units = Units()
        #
        self._cls = cls
        #self.type = material_type
        #
        self._labels: List[Union[str,int]] = []
        self._number : List[int] = array('I', [])        
        self._grade: List[str] = []
        self._Fy: List[float] = array('f', [])
        self._Fu: List[float] = array('f', [])
        self._poisson: List[float] = array('f', [])
        self._density: List[float] = array('f', [])
        self._E: List[float] = array('f', [])
        self._G: List[float] = array('f', [])
        self._alpha: List[float] = array('f', [])          
    #
    def __setitem__(self, material_name: Union[int, str],
                    properties: Union[List[float], Dict[str, float]]) -> None:
        """
        """
        try:
            self._labels.index(material_name)
            raise Exception('    *** warning material {:} already exist'
                            .format(material_name))
        except ValueError:
            number = next(self.get_number())
            self._labels.append(material_name)
            self._number.append(number)
            self._grade.append(-1)
            # default
            self._Fy.append(280_000_000)    # Pa
            self._Fu.append(0)    # Pa
            self._E.append(205_000_000_000) # Pa
            self._G.append(77_200_000_000)  # Pa
            self._poisson.append(0.30)
            self._density.append(7.850E12) # kg/m^3
            self._alpha.append(1.2E-5)  # K
            #
            index = self._labels.index(material_name)
            self._Fy[index] = properties[0]
    #
    def __getitem__(self, material_name: int) -> Tuple:
        """
        """
        try:
            _index = self._labels.index(material_name)
            return GetMaterial(self, material_name)
        except ValueError:
            raise IndexError('   *** material {:} does not exist'.format(material_name))        
    #
    def get_number(self, start:int=1)-> Iterable[int]:
        """
        """
        try:
            n = max(self._number) + 1
        except ValueError:
            n = start
        #
        while True:
            yield n
            n += 1
    #
#
#
class Spring(NamedTuple):
    """
    """
    #number: int
    force: List
    displacement: List
#
@dataclass
class Curve:
    """
    """
    __slots__ = ['type', '_spring', 'cls']
    
    def __init__(self, cls):
        """
        """
        self.cls = cls
    #
    @property
    def spring(self):
        """
        """
        _force = [_spring[1] for _spring in self._spring]
        _disp = [_spring[0] for _spring in self._spring]
        
        return Spring(force=_force, 
                      displacement=_disp)
    
    @spring.setter
    def spring(self, data):
        """
        """
        #self._force.append(data[0])
        #self._displacement.append(data[1])
        self._spring = data
#