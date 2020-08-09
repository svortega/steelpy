# 
# Copyright (c) 2019-2020 steelpy
# 

# Python stdlib imports
from collections.abc import Mapping
from typing import NamedTuple, Dict, List, Tuple, Union


# package imports
from steelpy.f2uModel.sections.shapes.tubular import Tubular
from steelpy.f2uModel.sections.shapes.tee import Tee
from steelpy.f2uModel.sections.shapes.ibeam import Ibeam
from steelpy.f2uModel.sections.shapes.solid import Trapeziod, Circle, Rectangle

#
class Sections(Mapping):
    
    __slots__ = ['_sections', '_default']
    
    def __init__(self):
        """
        alpha : Thermal expansion ratio (in 1/K at 20C )
        """
        #  ----- Section -----
        self._sections : Dict = {}
        self._default :Union[str,None] = None
    #
    #
    def __setitem__(self, shape_name, shape_type):
        """
        """
        item = shape_type.lower()
        
        if item in ['ibeam', 'isection', 'i', 
                    'w', 'm', 's', 'hp',
                    'ub', 'uc', 'he', 'ipe']:
            self._sections[shape_name] = Ibeam(cls=self)
        elif item in ['tee', 't']:
            self._sections[shape_name] = Tee(cls=self)
        elif item in ['tubular', 'pipe']:
            self._sections[shape_name] = Tubular(cls=self)
        elif item in ['rectangle']:
            self._sections[shape_name] = Rectangle(cls=self)        
        else:
            raise Exception(" section item {:} not recognized".format(shape_type))
            #print('fix here')        
        #print('-->')
        self._sections[shape_name].name = shape_name
        self._sections[shape_name].number = len(self._sections)
    
    def __getitem__(self, shape_name):
        """
        """
        if shape_name in self._sections:
            return self._sections[shape_name]
        else:
            raise Exception(" section name {:} not found".format(shape_name))   
    #
    def get_item_by_number(self, shape_name):
        """
        """
        _items = {_item.number:key for key, _item in self._sections.items()}
        try:
            _name = _items[shape_name]
            return self.__getitem__(_name)
        except KeyError:
            raise KeyError('Invalid section number')
    #
    #
    def __delitem__(self, shape_name: str) -> None:
        """
        """
        #_number = self._sections[shape_name].number
        del self._sections[shape_name]
    
    def __len__(self):
        return len(self._sections)
    
    def __iter__(self):
        """
        """
        return iter(self._sections)    
    
    def __contains__(self, value):
        return value in self._sections
    #
    def get_properties(self):
        """
        """
        summary = {}
        for key, item in self._sections.items():
            summary[key] = item._get_properties()
        return summary