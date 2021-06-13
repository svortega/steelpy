# 
# Copyright (c) 2009-2021 fem2ufo
# 


# Python stdlib imports
from collections.abc import Mapping
from typing import List, Iterable, Union, Dict
#

# package imports
#import steelpy.f2uModel.material.operations as operations
from steelpy.f2uModel.material.mechanical import MaterialElastic, Curve
#
#
class MaterialInmemory(Mapping):
    """
    """
    __slots__ =  ['_type', '_labels', '_number',
                  '_elastic', '_curve', '_default']
    
    def __init__(self):
        """
        """
        self._default:Union[str,None] = None
        self._labels:List[Union[str,int]] = []
        self._type:List[Union[str,int]] = []
        self._number:List[int] = []
        #
        self._elastic = MaterialElastic()
        self._curve = Curve(self)
    
    def __setitem__(self, material_name: Union[str, int],
                    properties: List[Union[str, float]]) -> None:
        """
        """
        try:
            self._labels.index(material_name)
            raise IOError('   error material {:} already exist'.format(material_name))
        except ValueError:
            material_type = properties[0]
            self._labels.append(material_name)
            self._type.append(material_type)
            mat_number = next(self.get_number())
            self._number.append(mat_number)
            #
            if 'curve' == material_type :
                raise Exception('--> Mat type No ready')
                #self._curve[mat_number]
            elif 'elastic' == material_type :
                self._elastic[mat_number] = [material_name, *properties[1:]]
            else:
                raise IOError(' material type {:} not recognised'
                              .format(material_type))
    
    def __getitem__(self, material_name:Union[str, int]):
        """
        """
        try:
            index = self._labels.index(material_name)
            material_type = self._type[index]
            mat_number = self._number[index]
        except ValueError:
            raise KeyError('Invalid material name : {:}'.format(material_name))
        #
        if 'curve' == material_type :
            return self._curve[mat_number]
        elif 'elastic' == material_type :
            return self._elastic[mat_number]
        else:
            raise IOError(' material type {:} not recognised'
                          .format(material_type))
    
    def __delitem__(self, material_name: Union[str, int]) -> None:
        """
        """
        _index = self._labels.index(material_name)
        1/0
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
        1/0
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