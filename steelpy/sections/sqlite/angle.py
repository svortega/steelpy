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
    __slots__ = ['_properties', # 'diameter', 'thickness', 
                 'name', 'number', 'db_file']
    
    def __init__(self, component: int, db_file:str):
                 #name:Union[str, int],
                 #d:Union[float,None], t:Union[float,None], 
                 #db_file:str,
                 #build:str = 'welded', 
                 #shear_stress:str = 'maximum',
                 #FAvy:float = 1.0, FAvz:float = 1.0):
        """
        Parameters
        ----------
        d : diametre
        t : wall Thickness
        Shear Stress: MAXIMUM / AVERAGE
        """
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
            bf = parameters[2]
            tf = parameters[3]            
            FAvy = 1
            FAvz = 1
            shear_stress:str = 'maximum'
            build:str = 'welded'            
            compactness = None
            section = (shape_name, 
                       "Angle",  # shape type
                       None, None,    # diameter, wall_thickess
                       d, tw,         # height, web_thickness
                       bf, tf,        # top_flange_width, top_flange_thickness
                       None, None,      # bottom_flange_width, bottom_flange_thickness
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
        return AngleBasic(name=row[0], 
                          d=row[5], tw=row[6],
                          b=row[7], r=0)
    #    