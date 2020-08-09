# 
# Copyright (c) 2009-2020 fem2ufo
#

# Python stdlib imports
from array import array
from collections.abc import Mapping
import re
from typing import NamedTuple, Tuple, List, Iterator, Dict, ClassVar


# package imports
from steelpy.f2uModel.properties.operations.operations import BasicProperty


#
class CDCMX:
    """
    """
    __slots__ = ('cd', 'cm', 'type',
                 'KC_number', 'profile')
    
    def __init__(self) -> None:
        """
        """
        self.cd = []
        self.cm = []
        self.KC_number = []
        self.type = None
    #
    def depth_profile(self, depth, cdcm):
        """
        """
        if self.type:
            self.profile.extend([depth, cdcm])
        else:
            self.type = 'profile'
            self.profile = [depth, cdcm]
#
class DragMassCoefficientXX:
    """
    type : specified, profile, diameter, rule, roughness/Reynolds number, roughness/KC number
    smooth  = [cd_Normal_x, cd_Longitudinal_y , cd_Longitudinal_z, cm_Normal_x, cm_Longitudinal_y , cm_Longitudinal_z]
    rough = [cd_Normal_x, cd_Longitudinal_y , cd_Longitudinal_z, cm_Normal_x, cm_Longitudinal_y , cm_Longitudinal_z]
    
    coefficient = [cd_Normal_x, cd_Longitudinal_y , cd_Longitudinal_z, cm_Normal_x, cm_Longitudinal_y , cm_Longitudinal_z]
    """
    __slots__ = ('number', 'name', 'items', 'sets', 'profile',
                 'type', 'zone', 'diameter', #'smooth', 'rough',
                 'Reynolds_number', 'KC_number', 'direction',
                 'concepts', 'coefficient')
    
    def __init__(self,  name, number):
        #
        self.number = number
        self.name = name
        self.coefficient = None
        #self.rough = None
        self.type = None
        #
        self.items = []
        self.sets = {}
    #     
    #
    #
    def depth_profile(self, depth, cdcm):
        """
        """
        #_coefficient = MorrisonCoefficient.set_cdcm(self, cdcm)
        
        if self.type:
            self.profile.extend([depth, cdcm])
        else:
            self.type = 'profile'
            self.profile = [depth, cdcm]
    #
    #
    def set_global_direction(self):
        """
        """
        pass
#
#
#
class Direction(NamedTuple):
    """
    Direction dependent coeffients
    """
    x:float
    y:float
    z:float    
#
class CdCm(NamedTuple):
    """
    Morison parameters CD & Cm
    """
    Cdx:float
    Cdy:float
    Cdz:float
    #
    Cmx:float
    Cmy:float
    Cmz:float
    #
    #number:int
    #sets:List
    #case:str
    #
    @property
    def drag(self):
        return Direction(self.Cdx, self.Cdy, self.Cdz)
    
    @property
    def mass(self):
        return Direction(self.Cmx, self.Cmy, self.Cmz)  
#
#
class CdCmParameters(Mapping):
    """
    """
    __slots__ = ('_labels', '_title',
                 '_cdx', '_cdy','_cdz',
                 '_cmx', '_cmy','_cmz',
                 '_case', '_sets')
    
    def __init__(self) -> None:
        """
        """
        self._labels = array('I', [])
        #
        self._cdx = array('f', [])
        self._cdy = array('f', [])
        self._cdz = array('f', [])
        #
        self._cmx = array('f', [])
        self._cmy = array('f', [])
        self._cmz = array('f', [])
        #
        self._title: List = []
        self._case: List = []
        self._sets: List = []
    
    def __getitem__(self, parameter_name) -> Tuple:
        """
        """
        try:
            _index = self._title.index(parameter_name)
            return CdCm(self._cdx[_index], self._cdy[_index], self._cdz[_index],
                        self._cmx[_index], self._cmy[_index], self._cmz[_index],
                        self._labels[_index], self._sets[_index], self._case[_index])
        except ValueError:
            raise KeyError(' cdcm {:} does not exist'.format(parameter_name))
    
    def __setitem__(self, parameter_name, cdcm) -> None:
        """
        """
        if not self.__contains__(parameter_name) :
            self._title.append(parameter_name)
            try:
                _number = max(self._labels) + 1
            except ValueError:
                _number = 1
            self._labels.append(_number)
            
            if isinstance(cdcm, (list, tuple)):
                if len(cdcm) == 2:
                    self._cdx.append(0)
                    self._cdy.append(cdcm[0])
                    self._cdz.append(cdcm[0])
                    #
                    self._cmx.append(0)
                    self._cmy.append(cdcm[1])
                    self._cmz.append(cdcm[1])                     
                elif len(cdcm) == 6:
                    self._cdx.append(cdcm[0])
                    self._cdy.append(cdcm[1])
                    self._cdz.append(cdcm[2])
                    #
                    self._cmx.append(cdcm[3])
                    self._cmy.append(cdcm[4])
                    self._cmz.append(cdcm[5])                    
                else:
                    raise IOError("cdcm data don't undesrtood")
                #
                self._sets.append(MeshGroupCases(label=_number, 
                                                 elements=[], nodes=[]))
                self._case.append('specified')
            else:
                raise IOError("cdcm data type no yet implemented")
        else:
            logging.warning(' ** cdcm {:} data updated'.format(parameter_name))
            _index = self._title.index(parameter_name)
            _sets = copy.copy(self.sets[_index])
            self.__delitem__(parameter_name)
            self.__setitem__(parameter_name, cdcm)
            self._sets[_index] = _sets
    #
    def get_item_name(self, number:int):
        """
        get item name by number
        """
        _index = self._labels.index(number)
        return self._title[_index]
    #
    def __delitem__(self, parameter_name) -> None:
        """
        """
        try:
            _index = self._title.index(parameter_name)
            self._title.remove(parameter_name)
            self._labels.pop(_index)
            self._cdx.pop(_index)
            self._cdy.pop(_index)
            self._cdz.pop(_index)
            #
            self._cmx.pop(_index)
            self._cmy.pop(_index)
            self._cmz.pop(_index)
            #
            self._sets.pop(_index)
            self._case.pop(_index)
        except IndexError:
            print('    *** warning cdcm property {:} does not exist'
                  .format(parameter_name))
            return
    
    def __len__(self) -> float:
        return len(self._labels)

    
    def __iter__(self) -> Iterator[Tuple]:
        return iter(self._title)
    
    def __contains__(self, value) -> bool:
        return value in self._title   
#
#
class CdCmCoefficients(Mapping):
    """
    """
    __slots__ =  ['_cdcm']
    
    def __init__(self):
        """
        """
        self._cdcm: Dict = {}
    
    def __getitem__(self, cdcm_name) -> ClassVar:
        """
        """
        return self._cdcm[cdcm_name]
    
    def __setitem__(self, cdcm_name, cdcm_type) -> None:
        """
        rule
        specified
        diametre
        """
        if re.match(r"\b((rule(\_)?)(api|iso))\b", cdcm_type, re.IGNORECASE):
            self._cdcm[cdcm_name] = DiametreFunction()
            self._cdcm[cdcm_name].rule = cdcm_type
        elif re.match(r"\b(specified)\b", cdcm_type, re.IGNORECASE):
            self._cdcm[cdcm_name] = SpecifiedFunction()
        elif re.match(r"\b(diamet(re|re))\b", cdcm_type, re.IGNORECASE):
            self._cdcm[cdcm_name] = DiametreFunction()
        elif re.match(r"\b(depth(\_)?profile)\b", cdcm_type, re.IGNORECASE):
            self._cdcm[cdcm_name] = DepthProfileFunction()
        else:
            raise IOError("CdCm type {:} not implemented".format(cdcm_type))         
    #
    #
    def __len__(self) -> float:
        return len(self._cdcm)

    
    def __iter__(self) -> Iterator[Tuple]:
        return iter(self._cdcm)
    
    def __contains__(self, value) -> bool:
        return value in self._cdcm     
#
#
class DiametreFunction(BasicProperty):
    __slots__ = ['_type', '_diameter', '_rule']
    
    def __init__(self):
        """
        """
        BasicProperty.__init__(self)
        self._type = 'diameter_function'
        self._diameter:Dict = {}
    
    def diameter_function(self, diameter, smooth, rough):
        """
        """
        self._diameter[diameter] = [smooth, rough]
    #
    @property
    def rule(self):
        """
        """
        return self._rule
    
    @rule.setter
    def rule(self, value:str):
        """
        """
        self._rule = value
#
class SpecifiedFunction(BasicProperty):
    __slots__ = ['_type', '_coefficients']
    
    def __init__(self):
        """
        """
        BasicProperty.__init__(self)
        self._type = 'specified_function'
    
    #@property
    def set_cdcm(self, cdcm):
        """
        """
        if type(cdcm) == list:
            self._coefficients = CdCm._make(cdcm)
        elif type(cdcm) == dict:
            self._coefficients = CdCm.Fixity(**cdcm)
        else:
            self._coefficients = CdCm(Cdx, Cdy, Cdz, Cmx, Cmy, Cmz)
#
class DepthProfileFunction(BasicProperty):
    __slots__ = ['_type', '_profile']
    
    def __init__(self):
        """
        """
        BasicProperty.__init__(self)
        self._type = 'depth_profile'
        self._profile:List = []
    #
    #
    def depth_profile(self, depth, cdcm):
        """
        """
        #_coefficient = MorrisonCoefficient.set_cdcm(self, cdcm)
        self._profile.extend([depth, cdcm])
   