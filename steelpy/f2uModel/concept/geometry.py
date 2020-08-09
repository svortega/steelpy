# 
# Copyright (c) 2009-2019 fem2ufo
# 


# Python stdlib imports
from collections.abc import Mapping
from typing import Dict, List, Tuple, ClassVar, Iterable, Union

# package imports



#
#
#  FE Classes
#
#
class SectionX:
    """
    Section property
    
    Section
        |_ name
        |_ number
        |_ type
        |_ shape
        |_ properties (object)
        |_ SFV shear factor vertical
        |_ SFH shear factor horizontal
        |_ elements [m1, m2, m3,..., mn]
        |_ length section's length
        |_ group 
        |_ step [section_1, section_2,...,section_n]
    
    **Parameters**:  
      :number:  integer internal number 
      :name:  string node external name
    """
    #
    __slots__ = ('number', 'name', 'type', 'shape',
                 'properties', 'SFV', 'SFH', 'length',
                 'elements', 'group', 'step', 'material')
    #
    def __init__(self, Name, Number, Type, Shape,
                 SFv=1.0, SFh=1.0):
        #
        self.number = Number
        self.name = Name
        self.type = Type
        self.shape = Shape
        self.SFV = SFv
        self.SFH = SFh
        self.elements = []
        #self.step = []
    #
    def equal(self, other):
        """
        """
        if other.shape == 'general section':
            
            if self.shape == 'general section':
                if self.properties.Area == other.properties.Area \
                and self.properties.Iip == other.properties.Iip \
                and self.properties.Iop == other.properties.Iop \
                and self.properties.J == other.properties.J \
                and self.properties.SCV == other.properties.SCV \
                and self.properties.SCH == other.properties.SCH :
                    return True
                else:
                    return False
            else:
                return False
        #
        elif self.shape == 'general section':
            
            if other.shape == 'general section':
                if self.properties.Area == other.properties.Area \
                and self.properties.Iip == other.properties.Iip \
                and self.properties.Iop == other.properties.Iop \
                and self.properties.J == other.properties.J \
                and self.properties.SCV == other.properties.SCV \
                and self.properties.SCH == other.properties.SCH :
                    return True
                else:
                    return False
            else:
                return False         
        #
        else:
            if self.properties.B1 == other.properties.B1 \
            and self.properties.B2 == other.properties.B2 \
            and self.properties.Bfb == other.properties.Bfb \
            and self.properties.Bft == other.properties.Bft \
            and self.properties.D == other.properties.D \
            and self.properties.Tfb == other.properties.Tfb \
            and self.properties.Tft == other.properties.Tft \
            and self.properties.Tw == other.properties.Tw :
                return True
            else:
                return False            
    #
    def get_name(self, length_unit, number=False):
        """
        """
        if not number:
            number = self.number
        
        if self.shape == 'general section':
            self.name = 'GB' + str(number).zfill(3)
        
        else:
            self.name = self.properties.get_name(length_unit, number)
#
#
class SupportX:
    """
    FE Fixity
    
    Boundary
        |_ name
        |_ number
        |_ type
        |_ constrain [x, y, z, mx, my, mz]
    
    **Parameters**:  
      :name : coupling, pile
      :type : free, fix, dependent, prescribed, supernode, master
      :constrain : 
            0- free
            1- fixed
            2- prescribed
            3- dependent
            4- supernode
    """
    
    __slots__ = ('number', 'name', 'type', 'releases',
                 'master', 'slave', 'nodes', 'elements',
                 'link', 'dependence')
    #
    fix_type = ['release', 'gap', 'prescribed', 'dependence',
                'link', 'master', 'rigid', 'constrain']
    #
    def __init__(self, Name, Number=None):
        """
        """
        self.number = Number
        self.name = Name
        self.type = 'free'
        self.slave = []
        self.nodes = []
        self.link = []
    #
    def set_releases(self, constrain):
        """
        """
        if type(constrain) == list:
            self.releases = geomop.Fixity._make(constrain)
        
        elif type(constrain) == dict:
            self.releases = geomop.Fixity(**constrain)
    #
    @property
    def fixed(self):
        """
        """
        constrain = [1, 1, 1, 1, 1, 1]
        self.releases = geomop.Fixity._make(constrain)
    #
    @property
    def pinned(self):
        """
        """
        constrain = [1, 1, 1, 0, 0, 0]
        self.releases = geomop.Fixity._make(constrain)
#
#
class HingesX(Mapping):
    __slots__ = ['_hinges', '_releases']
    
    def __init__(self) -> None:
        """
        """
        self._hinges:Dict = {}
        self._releases: ClassVar = Releases()
    
    def __setitem__(self, hinge_name: str, releases:List[int]) -> None:
        """
        """
        self._releases[hinge_name] = self._releases()
    
    def __getitem__(self, hinge_name: str)-> ClassVar:
        """
        """
        try:
            return self._hinges[hinge_name]
        except KeyError:
            raise KeyError('Invalid key')
    
