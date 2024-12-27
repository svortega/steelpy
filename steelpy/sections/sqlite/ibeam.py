# 
# Copyright (c) 2009 steelpy
#

# Python stdlib imports
from __future__ import annotations
#from collections import namedtuple
from typing import NamedTuple
#
#
# package imports
from steelpy.sections.utils.shape.ibeam import IbeamBasic
from steelpy.sections.sqlite.utils import SectionMainSQL
from steelpy.sections.utils.operations import get_sect_properties
from steelpy.utils.sqlite.utils import create_connection
#
#
#
#
class IbeamSQL(SectionMainSQL):
    __slots__ = ['_default', 'db_file', '_mesh_id']

    def __init__(self, mesh_id: int, db_file:str):
        """ 
        Parameters
        ----------
        d   : Section Height   
        tw  : Web thickness   
        bf  : flange base (top)
        tf  : flange thickness (top)
        bfb : Bottom flange base   
        tfb : Bottom flange thickness
        r   : root radius
        """
        super().__init__(mesh_id=mesh_id, db_file=db_file)
    #
    #
    #
    def __setitem__(self, shape_name: int|str,
                    parameters: list|NamedTuple) -> None:
        """
        parameters = []
        """
        try:
            self._labels.index(shape_name)
            raise Exception('element {:} already exist'.format(shape_name))
        except ValueError:
            section = (shape_name, 
                       parameters.shape, # shape type
                       None, None,       # diameter, wall_thickess
                       parameters.d,     # height,
                       parameters.tw,    # web_thickness
                       parameters.bf,    # top_flange_width,
                       parameters.tf,    # top_flange_thickness
                       parameters.bfb,   # bottom_flange_width,
                       parameters.tfb,   # bottom_flange_thickness
                       parameters.r,     # fillet_radius
                       *parameters[8:])  # FAvy, FAvz, shear_stress, build, compactness, title
            number = self.push_section(section)
            #self._number.append(number)
    #
    def __getitem__(self, shape_name: str | int):
        """
        """
        try:
            self._labels.index(shape_name)
            #number = self._number[index]
        except ValueError:
            raise Exception(f" section name {shape_name} not found")
        #
        conn = create_connection(self.db_file)
        with conn:        
            row = self._pull_section(conn, shape_name)
        #
        #1 / 0
        return IbeamBasic(name=row[0], 
                          d=row[5], tw=row[6],
                          bft=row[7], tft=row[8],
                          bfb=row[9], tfb=row[10],
                          root_radius=row[11])
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
    @property
    def tw(self):
        return self.get_item(item="web_thickness")

    @tw.setter
    def tw(self, value:float):
        """ """
        value = get_sect_properties([value])
        self.update_item(item='web_thickness', value=value[0])
        self.push_property()
    #
    #
    @property
    def bf(self):
        return self.get_item(item="top_flange_width")

    @bf.setter
    def bf(self, value:float):
        """ """
        value = get_sect_properties([value])
        self.update_item(item='top_flange_width', value=value[0])
        self.push_property()
    #
    @property
    def tf(self):
        return self.get_item(item="top_flange_thickness")

    @tf.setter
    def tf(self, value:float):
        """ """
        value = get_sect_properties([value])
        self.update_item(item='top_flange_thickness', value=value[0])
        self.push_property()
    #
    #
    @property
    def bfb(self):
        return self.get_item(item="bottom_flange_width")

    @bfb.setter
    def bfb(self, value:float):
        """ """
        value = get_sect_properties([value])
        self.update_item(item='bottom_flange_width', value=value[0])
        self.push_property()
    #
    @property
    def tfb(self):
        return self.get_item(item="bottom_flange_thickness")

    @tfb.setter
    def tfb(self, value:float):
        """ """
        value = get_sect_properties([value])
        self.update_item(item='bottom_flange_thickness', value=value[0])
        self.push_property()
    #    
#
