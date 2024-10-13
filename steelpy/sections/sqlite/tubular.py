# 
# Copyright (c) 2009 steelpy
#

# Python stdlib imports
from __future__ import annotations
#import re
#

# package imports
from steelpy.sections.utils.shape.tubular import TubularBasic
from steelpy.sections.sqlite.utils import SectionMainSQL
from steelpy.utils.sqlite.utils import create_connection
#
#from steelpy.utils.dataframe.main import DBframework
#
# ----------------------------------------
#      Standard Section Profiles
# ----------------------------------------
#
#
#
class TubularSQL(SectionMainSQL):
    """ """
    __slots__ = ['diameter', 'thickness', 'build', 'type',
                 'shear_stress', 'compactness', '_properties',
                 'FAvy', 'FAvz', 'name', 'number', 'db_file']
    
    def __init__(self, component:int, db_file:str):
        """
        Parameters
        ----------
        d : diametre
        t : wall Thickness
        Shear Stress: MAXIMUM / AVERAGE
        """
        self.db_file = db_file
        super().__init__(component=component,
                         db_file=self.db_file)
        #
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
            section = (shape_name, 
                       parameters.shape,  # shape type
                       parameters.d,      # diameter
                       parameters.tw,     # wall_thickness
                       None, None,        # height, web_thickness
                       None, None,        # top_flange_width, top_flange_thickness
                       None, None,        # bottom_flange_width, bottom_flange_thickness
                       None,              # root radius
                       *parameters[3:])   # FAvy, FAvz, shear_stress, build, compactness, title
            number = self.push_section(section)
    #
    def __getitem__(self, shape_name: str | int):
        """
        """
        try:
            index = self._labels.index(shape_name)
        except ValueError:
            raise Exception(f" section name {shape_name} not found")
        #
        conn = create_connection(self.db_file)
        with conn:        
            row = self._pull_section(conn, shape_name)        
        #
        return TubularBasic(name=row[0], 
                            diameter=row[3],
                            thickness=row[4])    
    #
    #