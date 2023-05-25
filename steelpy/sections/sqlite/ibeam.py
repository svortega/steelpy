# 
# Copyright (c) 2019-2023 steelpy
#

# Python stdlib imports
from __future__ import annotations
from collections import namedtuple
#from typing import Union, NamedTuple
#

# package imports
from ..inmemory.ibeam import IbeamBasic
from .operations import SectionSQLite
from steelpy.sections.process.operations import get_sect_properties
#
#

class IbeamSQLite(SectionSQLite):
    __slots__ = ['name', '_properties', 'db_file']

    def __init__(self, db_file:str):
                 #name:str|int,
                 #d: float, tw: float,
                 #bf: float, tf: float,
                 #db_file:str,
                 #bfb: float|None =None,
                 #tfb: float|None=None,
                 #build:str = 'welded',
                 #shear_stress:str = 'maximum',
                 #FAvy:float = 1.0, FAvz:float = 1.0,
                 #root_radius: float = 0):
        """ 
        Parameters
        ----------
        d   : Section Height   
        tw  : Web thickness   
        bf  : flange base (top)
        tf  : flange thickness (top)
        bfb : Bottom flange base   
        tfb : Bottom flange thickness
        """
        #IbeamBasic.__init__(d, tw, bft, tft, bfb, tfb)
        #self.name = name
        #self._properties = None
        self.db_file = db_file
        #
        #if not bfb: bfb = bf
        #if not tfb: tfb = tf
        #
        # push data to sqlite table
        #SectionSQLite.__init__(self, db_file=self.db_file, section=section)
        super().__init__(db_file=self.db_file)
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
            self._labels.append(shape_name)
            #
            d = parameters[0]
            tw = parameters[1]
            bf = parameters[2]
            tf = parameters[3]
            if not (bfb := parameters[4]):
                bfb = bf
            if not (tfb := parameters[5]):
                tfb = tf
            #if not (tfb := parameters[5]):
            #    tfb = tf
            FAvy = 1
            FAvz = 1
            shear_stress:str = 'maximum'
            build:str = 'welded'
            compactness = None
            section = (shape_name,
                       None,  # title
                       "I section",   # shape type
                       None, None,    # diameter, wall_thickess
                       d, tw,         # height, web_thickness
                       bf, tf,        # top_flange_width, top_flange_thickness
                       bfb, tfb,      # bottom_flange_width, bottom_flange_thickness
                       FAvy, FAvz,
                       shear_stress, build,
                       compactness,)
            number = self.push_section(section)
            self._number.append(number)
    #
    def __getitem__(self, shape_name: str | int):
        """
        """
        try:
            index = self._labels.index(shape_name)
            number = self._number[index]
        except ValueError:
            raise Exception(f" section name {shape_name} not found")
        #
        row = self.get_section(number)
        #
        return IbeamBasic(name=row[0], 
                          d=row[5], tw=row[6],
                          bft=row[7], tft=row[8],
                          bfb=row[9], tfb=row[10])
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
