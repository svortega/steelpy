# 
# Copyright (c) 2009 steelpy
#
#
# Python stdlib imports
from __future__ import annotations
#from collections.abc import Mapping
#import re
#
#
# package imports
from steelpy.utils.sqlite.utils import create_connection, create_table
#
from .angle import AngleSQL
from .tubular import TubularSQL
from .tee import TeeSQLite
from .channel import ChannelSQL
from .box import BoxSQL
from .ibeam import IbeamSQL
from .solid import SolidSectionSQL
#
from steelpy.sections.utils.operations import get_sect_df
from steelpy.sections.sqlite.utils import (SectionMainSQL,
                                           get_sections,
                                           get_section_df)

from steelpy.utils.dataframe.main import DBframework
#
#
#
class SectionSQL(SectionMainSQL):
    __slots__ = ['_tubular', '_solid', '_ibeam', '_box',
                 '_channel', '_tee', '_angle', '_default',
                 'db_file', '_mesh_id']
                 # '_labels', '_number', '_title', '_type',

    def __init__(self, db_file: str,
                 mesh_id:int, 
                 db_system: str = "sqlite"):
        """
        """
        super().__init__(mesh_id=mesh_id, db_file=db_file)
        #
        self._tubular = TubularSQL(mesh_id=mesh_id,
                                   db_file=db_file)
        
        self._solid = SolidSectionSQL(mesh_id=mesh_id,
                                      db_file=db_file)
        
        self._ibeam = IbeamSQL(mesh_id=mesh_id,
                               db_file=db_file)
        
        self._box = BoxSQL(mesh_id=mesh_id,
                           db_file=db_file)
        
        self._channel = ChannelSQL(mesh_id=mesh_id,
                                   db_file=db_file)
        #
        self._tee = TeeSQLite(mesh_id=mesh_id,
                              db_file=db_file)
        
        self._angle = AngleSQL(mesh_id=mesh_id,
                               db_file=db_file)
        #
        conn = create_connection(self.db_file)
        with conn:         
            self._new_table()
    #
    # -----------------------------------------------
    #
    def _new_table(self) -> None:
        """ """
        conn = create_connection(self.db_file)
        #
        table = "CREATE TABLE IF NOT EXISTS Section (\
                    number INTEGER PRIMARY KEY NOT NULL,\
                    name NOT NULL,\
                    mesh_id INTEGER NOT NULL REFERENCES Mesh(number),\
                    SA_inplane DECIMAL, \
                    SA_outplane DECIMAL,\
                    shear_stress TEXT, \
                    build TEXT, \
                    compactness TEXT, \
                    title TEXT);"
        create_table(conn, table)
        #
        # Section Geometry
        #
        table = "CREATE TABLE IF NOT EXISTS SectionGeometry (\
                    number INTEGER PRIMARY KEY NOT NULL,\
                    section_id INTEGER NOT NULL REFERENCES Section(number),\
                    type TEXT NOT NULL,\
                    diameter DECIMAL,\
                    wall_thickness DECIMAL,\
                    height DECIMAL,\
                    web_thickness DECIMAL,\
                    top_flange_width DECIMAL,\
                    top_flange_thickness DECIMAL,\
                    bottom_flange_width DECIMAL,\
                    bottom_flange_thickness DECIMAL,\
                    fillet_radius DECIMAL,\
                    web_orientation);"
        create_table(conn, table)
        #
        # Section Properties
        #
        table = "CREATE TABLE IF NOT EXISTS SectionProperty (\
                    number INTEGER PRIMARY KEY NOT NULL,\
                    section_id INTEGER NOT NULL REFERENCES Section(number),\
                    area DECIMAL,\
                    Zc DECIMAL,\
                    Yc DECIMAL,\
                    Iy DECIMAL,\
                    Zey DECIMAL,\
                    Zpy DECIMAL,\
                    ry DECIMAL,\
                    Iz DECIMAL,\
                    Zez DECIMAL,\
                    Zpz DECIMAL,\
                    rz DECIMAL,\
                    J DECIMAL,\
                    Cw DECIMAL);"
                    #alpha_sy DECIMAL,\
                    #alpha_sz DECIMAL);"
        
        create_table(conn, table)
    #
    # -----------------------------------------------
    #
    @property
    def df(self):
        """ raw data for dataframe"""
        conn = create_connection(self.db_file)
        with conn:        
            data = get_section_df(conn, mesh_id=self._mesh_id)
        db =  DBframework()
        header = ['name', 'type',
                  'diameter', 'wall_thickness',
                  'height', 'web_thickness',
                  'top_flange_width', 'top_flange_thickness',
                  'bottom_flange_width', 'bottom_flange_thickness',
                  'fillet_radius', 'web_orientation', 
                  'SA_inplane', 'SA_outplane',
                  'shear_stress', 'build', 'compactness', 'title',
                  'area[m^2]', 'Zc[m]', 'Yc[m]',
                  'Iy[m^4]', 'Zey[m^3]', 'Zpy[m^3]', 'ry[m]',
                  'Iz[m^4]', 'Zez[m^3]', 'Zpz[m^3]', 'rz[m]',
                  'J[m^4]', 'Cw[m^6]']
                  #'alpha_sy', 'alpha_sz']
        df = db.DataFrame(data=data, columns=header)
        return df

    @df.setter
    def df(self, df):
        """ """
        sectdf, propdf = get_sect_df(df)
        sectdf['mesh_id'] = self._mesh_id
        # Section
        header = ['name', 'mesh_id',
                  #'SA_inplane', 'SA_outplane',
                  'shear_stress', 'build',
                  'compactness', 'title']
        conn = create_connection(self.db_file)
        with conn:
            sectdf[header].to_sql('Section', conn,
                                  index_label=header,
                                  if_exists='append', index=False)
        # Section
        with conn:
            rows = get_sections(conn, mesh_id=self._mesh_id)
            sectid = {item[0]: item[1] for item in rows}
        #
        # SectionGeometry
        sectdf['section_id'] = sectdf['name'].apply(lambda x: sectid[x])
        header = ['section_id', 'type',
                  'diameter', 'wall_thickness',
                  'height', 'web_thickness',
                  'top_flange_width', 'top_flange_thickness',
                  'bottom_flange_width', 'bottom_flange_thickness',
                  'fillet_radius']
        conn = create_connection(self.db_file)
        with conn:
            sectdf[header].to_sql('SectionGeometry', conn,
                                  index_label=header,
                                  if_exists='append', index=False)
        #
        # SectionProperty
        propdf['section_id'] = propdf['name'].apply(lambda x: sectid[x])
        header = ['section_id',
                  'area', 'Zc', 'Yc',
                  'Iy', 'Zey', 'Zpy', 'ry', 
                  'Iz', 'Zez', 'Zpz', 'rz',
                  'J', 'Cw']
                  #'alpha_sy', 'alpha_sz']
        conn = create_connection(self.db_file)
        with conn:
            propdf[header].to_sql('SectionProperty', conn,
                                  index_label=header,
                                  if_exists='append', index=False)
        #

#


