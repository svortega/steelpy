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
from steelpy.sections.utils.shape.channel import ChannelBasic
from steelpy.sections.sqlite.utils import SectionMainSQL
from steelpy.sections.utils.operations import get_sect_properties
#
#
#
#
# ----------------------------------------
#      Standard Section Profiles
# ----------------------------------------

#
class ChannelSQL(SectionMainSQL):
    __slots__ = ['name', 'number', 'db_file', 
                 '_properties', '_mesh_id']

    def __init__(self, mesh_id: int, db_file: str):
        """ """
        super().__init__(mesh_id=mesh_id, db_file=db_file)
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
            #
            section = (shape_name,
                       "Channel",                   # shape type
                       None, None,                  # diameter, wall_thickess
                       parameters.d, parameters.tw, # height, web_thickness
                       parameters.b, parameters.tb, # top_flange_width, top_flange_thickness
                       parameters.b, parameters.tb, # bottom_flange_width, bottom_flange_thickness
                       parameters.r,                # root radius
                       *parameters[6:])             # FAvy, FAvz, shear_stress, build, compactness, title
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