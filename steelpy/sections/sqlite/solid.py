# 
# Copyright (c) 2009 steelpy
#

# Python stdlib imports
from __future__ import annotations
#from array import array
from dataclasses import dataclass
#from typing import NamedTuple
import re
#
#
# package imports
#
from steelpy.sections.sqlite.utils import SectionMainSQL
from steelpy.sections.utils.operations import get_sect_properties
from steelpy.sections.utils.shape.solid import (CircleSolid,
                                                RectangleSolid,
                                                TrapezoidSolid)
#
#
# ----------------------------------------
#      Basic Solid Shapes
# ----------------------------------------
#
#
@dataclass
class SolidSectionSQL(SectionMainSQL):
    
    def __init__(self, component, db_file:str):
        """ """
        self.db_file = db_file
        super().__init__(component=component,
                         db_file=self.db_file)
    #
    #
    def __setitem__(self, shape_name: int|str, parameters: list) -> None:
        """
        parameters : [circular bar, d, FAvy, FAvz, shear_stress, build, compactness, title]
                     [rectangle bar, h, b, FAvy, FAvz, shear_stress, build, compactness, title],
                     [trapezoid bar, h, b, a, FAvy, FAvz, shear_stress, build, compactness, title]
        """
        try:
            self._labels.index(shape_name)
            raise Exception('element {:} already exist'.format(shape_name))
        except ValueError:
            shape_type = parameters.shape
            #
            if re.match(r"\b((solid(_|-|\s*)?)?circular|round((_|-|\s*)?bar)?)\b", shape_type, re.IGNORECASE):
                section = (shape_name,
                           "Circular Bar",        # shape type
                           parameters.d, None,   # diameter, wall_thickness
                           None, None,            # height, web_thickness
                           None, None,            # top_flange_width, top_flange_thickness
                           None, None,            # bottom_flange_width, bottom_flange_thickness
                           None,                  # root radius
                           *parameters[4:])       # FAvy, FAvz, shear_stress, build, compactness, title

            elif re.match(r"\b((solid(_|-|\s*)?)?square|rectangle((_|-|\s*)?bar)?)\b", shape_type, re.IGNORECASE):
                section = (shape_name,
                           "Rectangle Bar",           # shape type
                           None, None,                # diameter, wall_thickness
                           parameters.h, None,        # height, web_thickness
                           parameters.b, None,        # top_flange_width, top_flange_thickness
                           parameters.b, None,        # bottom_flange_width, bottom_flange_thickness
                           None,                      # root radius
                           *parameters[4:])           # FAvy, FAvz, shear_stress, build, compactness, title

            elif re.match(r"\b((solid(_|-|\s*)?)?trapezoid((_|-|\s*)?bar)?)\b", shape_type, re.IGNORECASE):
                section = (shape_name,
                           "Trapezoid Bar",           # shape type
                           None, None,                # diameter, wall_thickness
                           parameters.h, None,       # height, web_thickness
                           parameters.b, None,       # top_flange_width, top_flange_thickness
                           parameters.a, None,       # bottom_flange_width, bottom_flange_thickness
                           None,                      # root radius
                           *parameters[4:])           # FAvy, FAvz, shear_stress, build, compactness, title

            else:
                raise Exception(f" section type {shape_type} not recognized")
            #
            number = self.push_section(section)
    #
    def __getitem__(self, shape_name: str | int):
        """
        """
        try:
            index = self._labels.index(shape_name)
            #number = self._number[index]
        except ValueError:
            raise Exception(f" section name {shape_name} not found")


        shape_type = self._type[index]
        row = self.get_section(shape_name)
        
        if re.match(r"\b((solid(_|-|\s*)?)?circular|round((_|-|\s*)?bar)?)\b", shape_type, re.IGNORECASE):
            d = row[3]
            return CircleSolid(name=shape_name, d=d)

        elif re.match(r"\b((solid(_|-|\s*)?)?square|rectangle((_|-|\s*)?bar)?)\b", shape_type, re.IGNORECASE):
            d = row[5]
            wb = row[7]
            return RectangleSolid(name=shape_name, depth=d, width=wb)

        elif re.match(r"\b((solid(_|-|\s*)?)?trapezoid((_|-|\s*)?bar)?)\b", shape_type, re.IGNORECASE):
            d = row[5]
            wb = row[7]
            wt = row[9]            
            c = abs(wt - wb) / 2.0
            return TrapezoidSolid(name=shape_name, depth=d, width=wb, a=wt, c=c)

        else:
            raise Exception(f" section type {shape_type} not recognized")
    #    
#
#
@dataclass
class RectangleSQLite(SectionMainSQL):
    __slots__ = ['_properties', 'name', 'number', 'db_file']

    def __init__(self, name:str|int,
                 d: float, w: float,
                 db_file:str,
                 build:str = 'welded',
                 shear_stress:str = 'maximum',
                 FAvy:float = 1.0, FAvz:float = 1.0):
        """
        Parameters
        ----------
        d : Height
        w : Width
        """
        #RectangleBasic.__init__(self)
        self.name = name
        self._properties = None
        self.db_file = db_file
        compactness = None
        section = (self.name,
                   "Rectangle",   # shape type
                   None, None,    # diameter, wall_thickess
                   d, None,       # height, web_thickness
                   w, None,       # top_flange_width, top_flange_thickness
                   w, None,       # bottom_flange_width, bottom_flange_thickness
                   FAvy, FAvz,
                   shear_stress, build,
                   compactness,
                   None,)     # title
        # push data to sqlite table
        SectionSQLite.__init__(self, db_file=self.db_file, section=section)
    #
    #
    @property
    def d(self):
        return self.get_item(item="height")

    @d.setter
    def d(self, value:float):
        """ """
        value = get_sect_properties([value])
        self.update_item(item='height', value=value[0])
        self.push_property()
    #
    #
    @property
    def w(self):
        return self.get_item(item="top_flange_width")

    @w.setter
    def w(self, value:float):
        """ """
        value = get_sect_properties([value])
        self.update_item(item='top_flange_width', value=value[0])
        self.push_property()
    #
#
#
@dataclass
class CircleSQLite(SectionMainSQL):
    __slots__ = ['_properties', 'name', 'number', 'db_file']
    
    def __init__(self, name:str|int,
                 d:float|None, 
                 db_file:str,
                 build:str = 'welded', 
                 shear_stress:str = 'maximum',
                 FAvy:float = 1.0, FAvz:float = 1.0):
        """
        Parameters
        ----------
        d : diametre
        Shear Stress: MAXIMUM / AVERAGE
        """
        #CircleBasic.__init__(self)
        self.name = name
        self._properties = None
        self.db_file = db_file
        compactness = None
        section = (self.name, 
                   "Circular Bar",  # shape type
                   d, None,         # diameter, wall_thickess
                   None, None, # height, web_thickness
                   None, None, # top_flange_width, top_flange_thickness
                   None, None, # bottom_flange_width, bottom_flange_thickness
                   FAvy, FAvz,
                   shear_stress, build,
                   compactness,
                   None,)     # title
        # push data to sqlite table
        #SectionSQLite.__init__(self, db_file=self.db_file,
        #                       section=section)
        super().__init__(db_file=self.db_file, section=section)
    #
    #
    @property
    def d(self):
        """
        D: diameter
        """
        return self.get_item(item="diameter")

    @d.setter
    def d(self, diameter: float):
        """
        """
        diameter = get_sect_properties([diameter])
        self.update_item(item='diameter', value=diameter[0])
        self.push_property()
    #    
#
#

#
#
#