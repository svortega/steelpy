# 
# Copyright (c) 2009 steelpy
#

# Python stdlib imports
from __future__ import annotations
from array import array
#from dataclasses import dataclass
#from collections import namedtuple
#import math
#import re
#

# package imports
from steelpy.sections.utils.shape.tubular import TubularBasic
from steelpy.sections.utils.shape.main import ShapeBasic


#
#
#
def find_tubular_dimensions(line_in: str) -> str:
    """
    """
    _key = {"diameter": r"\b(d(iamet(ro|er|re))?(y)?)\b",
            "thickness": r"\b((w(all)?(\_)?)?t(hickness|hk)?)\b"}

    keyWord, line_out, _match = search_line(line_in, _key)
    return keyWord


#
def get_dimension(self, dim: str, value: float):
    """
    """
    # Section Definition
    if dim == 'diameter':
        self.diameter = value
    elif dim == 'thickness':
        self.thickness = value
    else:
        raise IOError('   *** error Tubular geometry {:} not found'.format(dim))


#
def get_compactness(diameter: float, thickness: float) -> str:
    """
    """
    dt = diameter / thickness
    if dt > 80.0:
        compactness = 'slender'
        
    elif dt > 60.0:
        compactness = 'noncompact'
        
    else:
        compactness = 'compact'
        
    return compactness


#
#
# ----------------------------------------
#      Standard Section Profiles
# ----------------------------------------
#
#
class TubularIM(ShapeBasic):
    """ """
    __slots__ = ['_labels', '_number', '_title', '_d', '_tw']
    
    def __init__(self):
        """ """
        super().__init__()
        self._d: array = array('f', [])
        self._tw: array = array('f', [])        
        #self._process = TubularBasic()
    #
    #
    def __setitem__(self, shape_name: int|str, parameters: list) -> None:
        """
        parameters = [node1, node2, material, section, roll_angle]
        """
        try:
            self._labels.index(shape_name)
            raise Exception('element {:} already exist'.format(shape_name))
        except ValueError:
            self._labels.append(shape_name)
            self._title.append('NULL')
            mnumber = next(self.get_number())
            self._number.append(mnumber)
            #
            self._d.append(parameters.d)
            self._tw.append(parameters.tw)
            self._type.append('Tubular')
    #
    def __getitem__(self, shape_name: str | int):
        """
        """
        try:
            index = self._labels.index(shape_name)
        except ValueError:
            raise Exception(f" section name {shape_name} not found")
        #
        return TubularBasic(name=self._labels[index],
                            diameter=self._d[index],
                            thickness=self._tw[index])
    #
    @property
    def df(self):
        """ """
        df = DBframework()
        #stype = ['tubular' for _ in self._labels]
        data = {"Number": self._number,
                "Title": self._title,
                "Name": self._labels,
                "type": self._type,
                "d" : self._d,
                "tw" : self._tw}
        #print('-->')
        return df.DataFrame(data)
    
    @df.setter
    def df(self, df):
        """ """
        for item in df.name:
            mnumber =  next(self.get_number())
            self._number.append(mnumber)
        #mnumber = [next(self.get_number()) for _ in df.name]
        #self._number.extend(mnumber)
        self._title.extend(df.title.tolist())
        self._labels.extend(df.name.tolist())
        self._type.extend(df.type.tolist())
        try:
            self._d.extend(df.diameter.tolist())
        except TypeError:
            dtype = df['diameter'].apply(lambda x: x.convert('metre').value)
            self._d.extend(dtype.tolist())
        try:
            self._tw.extend(df.wall_thickness.tolist())
        except TypeError:
            dtype = df['wall_thickness'].apply(lambda x: x.convert('metre').value)
            self._tw.extend(dtype.tolist())
        #print('-->')
    #
    #
    #
    #
