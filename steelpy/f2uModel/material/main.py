# 
# Copyright (c) 2009-2020 fem2ufo
# 


# Python stdlib imports
from collections.abc import Mapping
from typing import NamedTuple, Dict, List, Iterable, Union
#from array import array

# package imports
import steelpy.f2uModel.material.operations as operations
from steelpy.f2uModel.material.mechanical import MaterialElastic, Curve
#from steelpy.process.units.units import Units
#
#
"""This module stores material classes.
   -------------------------------
   Nonlinear materials:
   - Nonlinear elastic material
   - Bilinear elastoplastic material
     Yield criterion:
      - von Mises
      - Treca
      - Morh-Coulomb
      - Drucker-Prager
     Hardering rules:
      - Isotropic
      - Kinematic
      - Isotropic + Kinematic
   - Multilinear plastic material
   - Rigid-plastic material
   
   -----------------------------
"""
#
# 
#
#
class Materials(Mapping):
    """
    """
    __slots__ =  ['_type', '_labels', '_number',
                  '_elastic', '_curve', '_default']
    
    def __init__(self):
        """
        """
        self._default = None
        self._labels:List[Union[str,int]] = []
        self._type:List[Union[str,int]] = []
        self._number:List[int] = []
        #
        self._elastic = MaterialElastic(self)
        self._curve = Curve(self)
    
    def __setitem__(self, material_name:Union[str, int], 
                    material_type:str) -> None:
        """
        """
        try:
            self._labels.index(material_name)
            raise IOError('   error material {:} already exist'.format(material_name))
        except ValueError:
            _material_type = operations.find_material_type(material_type)
            self._labels.append(material_name)
            self._type.append(_material_type)
            number = next(self.get_number())
            self._number.append(number)
            #
            if 'curve' == _material_type :
                self._curve[material_name]
            elif 'elastic' == _material_type :
                self._elastic[material_name] = [280_000_000] # default
            else:
                raise IOError(' material type {:} not recognised'
                              .format(material_type))
        #
    
    def __getitem__(self, material_name:str):
        """
        """
        try:
            _index = self._labels.index(material_name)
            _material_type = self._type[_index]
        except IndexError:
            raise KeyError('Invalid material name : {:}'.format(material_name))
        #
        if 'curve' == _material_type :
            return self._curve[material_name]
        elif 'elastic' == _material_type :
            return self._elastic[material_name]
        else:
            raise IOError(' material type {:} not recognised'
                          .format(material_type))        
    
    def __delitem__(self, material_name: str) -> None:
        """
        """
        #material_number = self._materials[material_name].number
        del self._labels[material_name]
        #del self._labels[material_number]
    
    def __len__(self):
        return len(self._labels)
    
    def __iter__(self):
        """
        """
        return iter(self._labels)
    
    def __contains__(self, value):
        return value in self._labels
    #
    # Modify material
    #
    #
    #@property
    def get_item_by_number(self, material_name):
        """
        """
        #_items = {_item.number:key for key, _item in self._materials.items()}
        try:
            #material_name = _items[material_number]
            return self.__getitem__(material_name)
        except KeyError:
            raise KeyError('Invalid material name')
    #
    def get_material(self):
        """
        """
        summary = {}
        for key, mat in self._labels.items():
            summary[key] = mat._get_data()
        return summary
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