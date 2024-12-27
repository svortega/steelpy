# 
# Copyright (c) 2009 steelpy
#

# Python stdlib imports
from __future__ import annotations
#from typing import NamedTuple
#

# package imports
from steelpy.sections.utils.shape.angle import AngleBasic
from steelpy.sections.sqlite.utils import SectionMainSQL

#
#
#
# ----------------------------------------
#      Standard Section Profiles
# ----------------------------------------
#

class AngleSQL(SectionMainSQL):
    """ """
    __slots__ = ['_properties', '_mesh_id', 
                 'name', 'number', 'db_file']
    
    def __init__(self, mesh_id: int, db_file:str):
        """
        Parameters
        ----------
        d : diametre
        t : wall Thickness
        Shear Stress: MAXIMUM / AVERAGE
        """
        super().__init__(mesh_id=mesh_id, db_file=db_file)
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
                       "Angle",                     # shape type
                       None, None,                  # diameter, wall_thickess
                       parameters.d, parameters.tw, # height, web_thickness
                       parameters.b, parameters.tb, # top_flange_width, top_flange_thickness
                       None, None,                  # bottom_flange_width, bottom_flange_thickness
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
        return AngleBasic(name=row[0], 
                          d=row[5], tw=row[6],
                          b=row[7], r=0)
    #    