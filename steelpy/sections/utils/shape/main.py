#
# Copyright (c) 2019 steelpy
#

# Python stdlib imports
from __future__ import annotations
#from dataclasses import dataclass
from collections.abc import Mapping
#import math
import re
#

# package imports
#
from steelpy.sections.utils.shape.utils import get_sect_properties, get_Isection
from steelpy.sections.utils.shape.tubular import TubularBasic
from steelpy.sections.utils.shape.solid import RectangleBasic, CircleBasic, Trapeziod
from steelpy.sections.utils.shape.ibeam import IbeamBasic
from steelpy.sections.utils.shape.box import BoxBasic
from steelpy.sections.utils.shape.channel import ChannelBasic
from steelpy.sections.utils.shape.tee import TeeBasic
from steelpy.sections.utils.shape.angle import AngleBasic
#
from steelpy.utils.dataframe.main import DBframework
#
#
#import numpy as np
#

#-------------------------------------------------
#
class ShapeBasic(Mapping):
    """ """
    #__slots__ = []
    
    def __init__(self):
        """
        """
        pass
    # -----------------------------------------------
    #
    def __len__(self):
        return len(self._labels)

    def __iter__(self):
        return iter(self._labels)

    def __contains__(self, value):
        return value in self._labels
    #
    #def __delitem__(self, shape_name: str | int) -> None:
    #    try:
    #        _nodes_empty = []
    #        _index = self._labels.index(shape_name)
    #        1/0
    #    except ValueError:
    #        raise KeyError(f'    *** warning section {shape_name} does not exist')
    #
    # -----------------------------------------------
    #
    #def properties(self):
    #    """
    #    """
    #    1/0
    #    return self._properties()    
    #
    #def get_number(self, start: int = 1):
    #    """
    #    """
    #    try:
    #        n = max(self._labels) + 1
    #    except ValueError:
    #        n = start
    #    #
    #    while True:
    #        yield n
    #        n += 1
    #
    
#
#
#
#-------------------------------------------------
#
#
class SectionMain(ShapeBasic):
    #__slots__ = ['_labels', '_number', '_title', '_type',
    #             '_tubular', '_solid', '_ibeam', '_box',
    #             '_channel', '_tee', '_angle', '_default']

    def __init__(self):
        """
        """
        super().__init__()
    #
    def __setitem__(self, shape_name: str | int,
                    properties: list[float] | dict[str, float] | str) -> None:
        """
        """
        try:
            self._labels.index(shape_name)
            raise Exception(f'Section {shape_name} already exist')
        except ValueError:
            #
            shape_type = properties[0]
            properties = get_sect_properties(properties[1:])
            #
            if re.match(r"\b(i((\_)?beam|section)?|w|m|s|hp|ub|uc|he|ipe|pg)\b",
                        shape_type, re.IGNORECASE):
                # [d, tw, bf, tf, bfb, tfb, r, title]
                properties =  get_Isection(properties)
                self._ibeam[shape_name] = properties

            elif re.match(r"\b(t(ee)?)\b", shape_type, re.IGNORECASE):
                self._tee[shape_name] = properties                  

            elif re.match(r"\b(tub(ular)?|pipe|chs)\b", shape_type, re.IGNORECASE):
                self._tubular[shape_name] = properties

            elif re.match(r"\b((solid|bar(\_)?)?rectangle|trapeziod|circular|round)\b",
                          shape_type, re.IGNORECASE):
                self._solid[shape_name] = [shape_type, *properties] 

            elif re.match(r"\b(b(ox)?|rhs|shs)\b", shape_type, re.IGNORECASE):
                self._box[shape_name] = properties 

            elif re.match(r"\b(c(hannel)?)\b", shape_type, re.IGNORECASE):
                self._channel[shape_name] = properties                 

            elif re.match(r"\b(l|angle)\b", shape_type, re.IGNORECASE):
                self._angle[shape_name] = properties               

            else:
                raise Exception(" section item {:} not recognized".format(shape_type))    
    #
    def __getitem__(self, shape_name: int|str):
        """
        node_name : node number
        """
        try:
            index = self._labels.index(shape_name)
            shape_type = self._type[index]
        except ValueError:
            raise KeyError(f'   *** Section {shape_name} does not exist')
        #
        if re.match(r"\b(tub(ular)?|pipe)\b", shape_type, re.IGNORECASE):
            return self._tubular[shape_name]

        elif re.match(r"\b((solid|bar(\_)?)?rectangle|trapeziod|circular|round)\b", shape_type, re.IGNORECASE):
            return self._solid[shape_name]
        
        elif re.match(r"\b(i((\_)?beam|section)?|w|m|s|hp|ub|uc|he|ipe|pg)\b", shape_type, re.IGNORECASE):
            return self._ibeam[shape_name]
        
        elif re.match(r"\b(b(ox)?|rhs|shs)\b", shape_type, re.IGNORECASE):
            return self._box[shape_name]
        
        elif re.match(r"\b(c(hannel)?)\b", shape_type, re.IGNORECASE):
            return self._channel[shape_name]
        
        elif re.match(r"\b(t(ee)?)\b", shape_type, re.IGNORECASE):
            return self._tee[shape_name]
        
        elif re.match(r"\b(l|angle)\b", shape_type, re.IGNORECASE):
            return self._angle[shape_name]
        
        else:
            raise IOError(f' Section type {shape_type} not recognised')
    # 
    #
    # -----------------------------------------------
    #
    @property
    def df(self):
        """ """
        db = DBframework()
        #header = ['number', 'name', 'title', 'type', 
        #          'diameter', 'wall_thickness',
        #          'height', 'web_thickness',
        #          'top_flange_width', 'top_flange_thickness',
        #          'bottom_flange_width', 'bottom_flange_thickness',
        #          'fillet_radius',
        #          'SA_inplane', 'SA_outplane',
        #          'shear_stress', 'build', 'compactness']
        #
        number = [x + 1 for x, item in enumerate(self._labels)]
        #
        datadf = {'number': number,
                  'name' : self._labels,
                  'title' : self._title,
                  'type': self._type, 
                  'diameter': [None for item in self._labels],
                  'wall_thickness': [None for item in self._labels],
                  'height': [None for item in self._labels],
                  'web_thickness': [None for item in self._labels],
                  'top_flange_width': [None for item in self._labels],
                  'top_flange_thickness': [None for item in self._labels],
                  'bottom_flange_width': [None for item in self._labels],
                  'bottom_flange_thickness': [None for item in self._labels],
                  'fillet_radius': [None for item in self._labels],
                  'SA_inplane': [1.0 for item in self._labels],
                  'SA_outplane': [1.0 for item in self._labels],
                  'shear_stress': ['maximum' for item in self._labels],
                  'build': ['welded' for item in self._labels],
                  'compactness': [None for item in self._labels]}        
        #
        #sec_df = db.DataFrame(data)
        #
        #for row in sec_df.itertuples():
        #    section = self.__getitem__(shape_name=row.name)
        #    data = section._data_df()
        #    sec_df[data.keys()].iloc[row.Index] =  data.values()
        #    #for key, item in data.items():
        #    #    sec_df[key].iloc[row.Index] = item
        #
        for idx, name in enumerate(self._labels):
            section = self.__getitem__(shape_name=name)
            data = section._data_df()
            for key, item in data.items():
                datadf[key][idx] = item
        #
        return db.DataFrame(datadf)
    #
    @df.setter
    def df(self, df):
        """ """
        #1 / 0
        #
        group = df.groupby("type", sort=False)
        for shape_type, section in group:
            if re.match(r"\b(i((\_)?beam|section)?|w|m|s|hp|ub|uc|he|ipe|pg)\b",
                        shape_type, re.IGNORECASE):
                # Bottom Flange
                try:
                    section['bottom_flange_width'] =  section.apply(lambda x: x['top_flange_width']
                                                                    if x['bottom_flange_width']== ""
                                                                    else x['bottom_flange_width'], axis=1)
                except KeyError:
                    section['bottom_flange_width'] = section['top_flange_width']
                #
                try:
                    section['bottom_flange_thickness'] =  section.apply(lambda x: x['top_flange_thickness']
                                                                    if x['bottom_flange_thickness']== ""
                                                                    else x['bottom_flange_thickness'], axis=1)
                except KeyError:
                    section['bottom_flange_thickness'] = section['top_flange_thickness']                
                #
                try:
                    section['fillet_radius'] =  section.apply(lambda x: float(0.0)
                                                                    if x['fillet_radius']== ""
                                                                    else x['fillet_radius'], axis=1)
                except KeyError:
                    section['fillet_radius'] = float(0.0)                
                #
                self._ibeam.df= section
            
            elif re.match(r"\b(t(ee)?)\b", shape_type, re.IGNORECASE):
                self._tee.df= section
            
            elif re.match(r"\b(tub(ular)?|pipe|chs)\b", shape_type, re.IGNORECASE):
                self._tubular.df= section
            
            elif re.match(r"\b((solid|bar(\_)?)?rectangle|trapeziod|circular|round)\b",
                           shape_type, re.IGNORECASE):
                self._solid.df= section
            
            elif re.match(r"\b(b(ox)?|rhs|shs)\b", shape_type, re.IGNORECASE):
                self._box.df= section
            
            elif re.match(r"\b(c(hannel)?)\b", shape_type, re.IGNORECASE):
                self._channel.df= section
            
            elif re.match(r"\b(l|angle)\b", shape_type, re.IGNORECASE):
                self._angle.df= section
            
            else:
                raise Exception(" section item {:} not recognized".format(shape_type))
    #     
#
#
#-------------------------------------------------
#
#
class SectionBasicXX(ShapeBasic):
    #__slots__ = ['_labels', '_number', '_title', '_type',
    #             '_tubular', '_solid', '_ibeam', '_box',
    #             '_channel', '_tee', '_angle', '_default']

    def __init__(self):
        """
        """
        super().__init__()
        self._default: str | None = None
        self._labels: list[str|int] = []
        self._number: list[int] = []
        self._title: list[str] = []
        self._type: list = []
    #
    def __setitem__(self, shape_name: str | int,
                    properties: list[float] | dict[str, float] | str) -> None:
        """
        """
        try:
            self._labels.index(shape_name)
            raise Exception(f'Section {shape_name} already exist')
        except ValueError:
            #
            shape_type = properties[0]
            properties = get_sect_properties(properties[1:])
            #
            self._labels.append(shape_name)
            self._type.append(shape_type)
            #mnumber = next(self.get_number())
            #self._number.append(mnumber)            
            #
            if re.match(r"\b(i((\_)?beam|section)?|w|m|s|hp|ub|uc|he|ipe|pg)\b",
                        shape_type, re.IGNORECASE):
                # [d, tw, bf, tf, bfb, tfb, r, title]
                properties =  get_Isection(properties)
                self._ibeam[shape_name] = properties

            elif re.match(r"\b(t(ee)?)\b", shape_type, re.IGNORECASE):
                self._tee[shape_name] = properties                  

            elif re.match(r"\b(tub(ular)?|pipe|chs)\b", shape_type, re.IGNORECASE):
                self._tubular[shape_name] = properties

            elif re.match(r"\b((solid|bar(\_)?)?rectangle|trapeziod|circular|round)\b",
                          shape_type, re.IGNORECASE):
                self._solid[shape_name] = [shape_type, *properties] 

            elif re.match(r"\b(b(ox)?|rhs|shs)\b", shape_type, re.IGNORECASE):
                self._box[shape_name] = properties 

            elif re.match(r"\b(c(hannel)?)\b", shape_type, re.IGNORECASE):
                self._channel[shape_name] = properties                 

            elif re.match(r"\b(l|angle)\b", shape_type, re.IGNORECASE):
                self._angle[shape_name] = properties               

            else:
                raise Exception(" section item {:} not recognized".format(shape_type))    
    #
    def __getitem__(self, shape_name: int):
        """
        node_name : node number
        """
        try:
            index = self._labels.index(shape_name)
            shape_type = self._type[index]
        except ValueError:
            raise KeyError(f'   *** Section {shape_name} does not exist')
        #
        if re.match(r"\b(tub(ular)?|pipe)\b", shape_type, re.IGNORECASE):
            return self._tubular[shape_name]

        elif re.match(r"\b((solid|bar(\_)?)?rectangle|trapeziod|circular|round)\b", shape_type, re.IGNORECASE):
            return self._solid[shape_name]
        
        elif re.match(r"\b(i((\_)?beam|section)?|w|m|s|hp|ub|uc|he|ipe|pg)\b", shape_type, re.IGNORECASE):
            return self._ibeam[shape_name]
        
        elif re.match(r"\b(b(ox)?|rhs|shs)\b", shape_type, re.IGNORECASE):
            return self._box[shape_name]
        
        elif re.match(r"\b(c(hannel)?)\b", shape_type, re.IGNORECASE):
            return self._channel[shape_name]
        
        elif re.match(r"\b(t(ee)?)\b", shape_type, re.IGNORECASE):
            return self._tee[shape_name]
        
        elif re.match(r"\b(l|angle)\b", shape_type, re.IGNORECASE):
            return self._angle[shape_name]
        
        else:
            raise IOError(f' Section type {shape_type} not recognised')
    #    
    #
    #
    #
#
#
#-------------------------------------------------
#
#
#
#@dataclass
#class ShapeGeometry:
#    """ """
#    name: str | int
#    #shape_type: str
#    #row: list
#    #
def ShapeGeometry(shape_type: str, geometry: list):
        """ """
        if re.match(r"\b(tub(ular)?|pipe)\b", shape_type, re.IGNORECASE):
            return TubularBasic(name=geometry[0], 
                                diameter=geometry[3], thickness=geometry[4]) 
    
        #elif re.match(r"\b((solid|bar(\_)?)?rectangle|trapeziod|circular|round)\b", shape_type, re.IGNORECASE):
        #    return self._solid[shape_name]
        elif re.match(r"\b((solid|bar(\_)?)?circular|round)\b", shape_type, re.IGNORECASE):
            d = geometry[3]
            return CircleBasic(name=geometry[0], d=d, type=shape_type)
    
        elif re.match(r"\b((solid|bar(\_)?)?rectangle)\b", shape_type, re.IGNORECASE):
            d = geometry[3]
            wb = geometry[7]
            return RectangleBasic(name=geometry[0], depth=d, width=wb,
                                  type=shape_type)
    
        elif re.match(r"\b((solid|bar(\_)?)?trapeziod)\b", shape_type, re.IGNORECASE):
            d = geometry[5]
            wb = geometry[7]
            wt = geometry[9]            
            c = abs(wt - wb) / 2.0
            return Trapeziod(name=geometry[0], depth=d, width=wb,
                             a=wt, c=c, type=shape_type)    
        
        elif re.match(r"\b(i((\_)?beam|section)?|w|m|s|hp|ub|uc|he|ipe|pg)\b", shape_type, re.IGNORECASE):
            return IbeamBasic(name=geometry[0], 
                              d=geometry[5], tw=geometry[6],
                              bft=geometry[7], tft=geometry[8],
                              bfb=geometry[9], tfb=geometry[10])
        
        elif re.match(r"\b(b(ox)?|rhs|shs)\b", shape_type, re.IGNORECASE):
            return BoxBasic(name=geometry[0], 
                    d=geometry[5], tw=geometry[6],
                    b=geometry[7], tb=geometry[8])
        
        elif re.match(r"\b(c(hannel)?)\b", shape_type, re.IGNORECASE):
            return ChannelBasic(name=geometry[0], 
                        d=geometry[5], tw=geometry[6],
                        b=geometry[7], tb=geometry[8])
        
        elif re.match(r"\b(t(ee)?)\b", shape_type, re.IGNORECASE):
            return TeeBasic(name=geometry[0], 
                            d=geometry[5], tw=geometry[6],
                            b=geometry[7], tb=geometry[8])
        
        elif re.match(r"\b(l|angle)\b", shape_type, re.IGNORECASE):
            return AngleBasic(name=geometry[0], 
                              d=geometry[5], tw=geometry[6],
                              b=geometry[7], r=0)
        
        else:
            raise IOError(f' Section type {shape_type} not recognised')
# 
