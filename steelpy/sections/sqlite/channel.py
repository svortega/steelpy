# 
# Copyright (c) 2019-2023 steelpy
#
#
# Python stdlib imports
from __future__ import annotations
#import math
#from typing import Union
#import sys
#

# package imports
#from ..inmemory.channel import ChannelBasic
from steelpy.sections.sqlite.utils import SectionItemSQL
#from steelpy.sections.process.operations import  get_sect_properties

#
#
#
#
# ----------------------------------------
#      Standard Section Profiles
# ----------------------------------------

#
class ChannelSQLite(SectionItemSQL):
    __slots__ = ['name', 'number', 'db_file', '_properties']

    def __init__(self, component: int, db_file: str):
        """ """
        self.db_file = db_file
        # push data to sqlite table
        super().__init__(component=component,
                         db_file=self.db_file)
    #
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
                       "Channel",  # shape type
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
        return ChannelBasic(name=row[0], 
                            d=row[5], tw=row[6],
                            b=row[7], tb=row[8])     
    #
    @property
    def d(self):
        """ """
        return self.get_item ( item="height" )

    @d.setter
    def d(self, diameter: Units|float ):
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
    def tw(self, thickness: Units|float):
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
    def b(self, diameter: Units|float):
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
    def tb(self, thickness: Units|float):
        """
        """
        thickness = get_sect_properties ( [ thickness ] )
        self.update_item ( item='top_flange_thickness', value=thickness[ 0 ] )
        self.push_property ()
    #
#
#
#