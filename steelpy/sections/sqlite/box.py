# 
# Copyright (c) 2009 steelpy
#
#
# Python stdlib imports
from __future__ import annotations
#from typing import Union
#
#
# package imports
#
from steelpy.sections.sqlite.utils import SectionMainSQL
from steelpy.sections.utils.shape.box import BoxBasic
#
#
#
#
# ----------------------------------------
#      Standard Section Profiles
# ----------------------------------------

class BoxSQL(SectionMainSQL):
    __slots__ = [ '_properties',
                  'name', 'number', 'db_file' ]

    def __init__(self, component: int, db_file: str):
        """ """
        self.db_file = db_file
        # push data to sqlite table
        super().__init__(component=component,
                         db_file=self.db_file)
    #
    def __setitem__(self, shape_name: int|str, parameters: list) -> None:
        """
        parameters = []
        """
        try:
            self._labels.index(shape_name)
            raise Exception('element {:} already exist'.format(shape_name))
        except ValueError:
            #self._labels.append(shape_name)
            #
            d = parameters[0]
            tw = parameters[1]
            b = parameters[2]
            tb = parameters[3]
            FAvy = 1
            FAvz = 1
            shear_stress:str = 'maximum'
            build:str = 'welded'            
            compactness = None
            section = (shape_name,
                       "Box",  # shape type
                       None, None, # diameter, wall_thickess
                       d, tw,      # height, web_thickness
                       b, tb,      # top_flange_width, top_flange_thickness
                       b, tb,      # bottom_flange_width, bottom_flange_thickness
                       None,       # root radius
                       FAvy, FAvz,
                       shear_stress, build,
                       compactness,
                       None,)     # title
            number = self.push_section(section)
            #self._number.append(number)
    #
    def __getitem__(self, shape_name: str | int):
        """
        """
        try:
            index = self._labels.index(shape_name)
            #number = self._number[index]
        except ValueError:
            raise Exception(f" section name {shape_name} not found")
        #
        row = self.get_section(shape_name)
        return BoxBasic(name=row[0], 
                        d=row[5], tw=row[6],
                        b=row[7], tb=row[8])        
    #
    @property
    def d(self):
        """
        D: diameter
        """
        return self.get_item ( item="height" )

    @d.setter
    def d(self, diameter: Union[ Units, float ]):
        """
        """
        diameter = get_sect_properties ( [ diameter ] )
        self.update_item ( item='height', value=diameter[ 0 ] )
        self.push_property ()

    #
    @property
    def tw(self):
        """
        """
        return self.get_item ( item="web_thickness" )

    @tw.setter
    def tw(self, thickness: Union[ Units, float ]):
        """
        """
        thickness = get_sect_properties ( [ thickness ] )
        self.update_item ( item='web_thickness', value=thickness[ 0 ] )
        self.push_property ()
    #
    #
    #
    @property
    def b(self):
        """
        D: diameter
        """
        return self.get_item ( item="top_flange_width" )

    @b.setter
    def b(self, diameter: Union[ Units, float ]):
        """
        """
        diameter = get_sect_properties ( [ diameter ] )
        self.update_item ( item='top_flange_width', value=diameter[ 0 ] )
        self.push_property ()

    #
    @property
    def tb(self):
        """
        """
        return self.get_item ( item="top_flange_thickness" )

    @tb.setter
    def tb(self, thickness: Union[ Units, float ]):
        """
        """
        thickness = get_sect_properties ( [ thickness ] )
        self.update_item ( item='top_flange_thickness', value=thickness[ 0 ] )
        self.push_property ()
    #
    #
    #
#
#