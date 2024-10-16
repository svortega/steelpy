# 
# Copyright (c) 2009 steelpy
#
#
# Python stdlib imports
from __future__ import annotations
#from collections.abc import Mapping
import re
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
from steelpy.sections.sqlite.utils import SectionMainSQL
#from steelpy.utils.dataframe.main import DBframework
#
#
#
class SectionSQL(SectionMainSQL):
    __slots__ = ['_tubular', '_solid', '_ibeam', '_box',
                 '_channel', '_tee', '_angle', '_default',
                 'db_file', '_component']
                 # '_labels', '_number', '_title', '_type',

    def __init__(self, db_file: str,
                 component:int, 
                 db_system: str = "sqlite"):
        """
        """
        super().__init__(component=component,
                         db_file=db_file)
        #
        #self.db_file = db_file
        #self._component = component
        #
        self._tubular = TubularSQL(component=component,
                                   db_file=db_file)
        
        self._solid = SolidSectionSQL(component=component,
                                      db_file=db_file)
        
        self._ibeam = IbeamSQL(component=component,
                               db_file=db_file)
        
        self._box = BoxSQL(component=component,
                           db_file=db_file)
        
        self._channel = ChannelSQL(component=component,
                                   db_file=db_file)
        #
        self._tee = TeeSQLite(component=component,
                              db_file=db_file)
        
        self._angle = AngleSQL(component=component,
                               db_file=db_file)
        #
        conn = create_connection(self.db_file)
        with conn:         
            self._new_table()
    #
    #
    # def push_sections(self):
    #    """
    #    """
    #    conn = create_connection(self.bd_file)
    #    with conn:
    #        for key, section in self._sections.items():
    #            sect_number = self._push_section_table(conn, section)
    #            section.number = sect_number
    #            self._push_property_table(conn, section)
    #        conn.commit()
    #
    # def _push_property_table(self, conn, section):
    #    """ """
    #    project = (section.number, *section.properties)
    #    sql = 'INSERT INTO  SectionProperty(number, area, Zc, Yc,\
    #                                        Iy, Zey, Zpy, ry,\
    #                                        Iz, Zez, Zpz, rz,\
    #                                        J, Cw)\
    #                                        VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?)'
    #    cur = conn.cursor()
    #    cur.execute(sql, project)
    #
    # def _push_section_table(self, conn, section) -> int:
    #    """
    #    """
    #    project = section._get_section_table()
    #    sql = 'INSERT INTO  Section(name, title, type, diameter, wall_thickness,\
    #                                    height, web_thickness,\
    #                                    top_flange_width, top_flange_thickness,\
    #                                    bottom_flange_width, bottom_flange_thickness)\
    #                                    VALUES(?,?,?,?,?,?,?,?,?,?,?)'
    #    cur = conn.cursor()
    #    cur.execute(sql, project)
    #    return cur.lastrowid
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
                    fillet_radius DECIMAL);"
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
                    Cw DECIMAL,\
                    alpha_sy DECIMAL,\
                    alpha_sz DECIMAL);"
        
        create_table(conn, table)
    #
    #
    # -----------------------------------------------
    #
    @property
    def df(self):
        """ raw data for dataframe"""
        1/0
        #return self._sections.df

    @df.setter
    def df(self, df):
        """ """
        sectdf, propdf = get_sect_df(df)
        sectdf['mesh_id'] = self._component
        # Section
        header = ['name', 'mesh_id',
                  'SA_inplane', 'SA_outplane',
                  'shear_stress', 'build',
                  'compactness', 'title']
        conn = create_connection(self.db_file)
        with conn:
            sectdf[header].to_sql('Section', conn,
                                  index_label=header,
                                  if_exists='append', index=False)
        # Section
        with conn:
            rows = get_sections(conn, component=self._component)
            sectid = {item[1]: item[0] for item in rows}
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
                  'J', 'Cw',
                  'alpha_sy', 'alpha_sz']
        conn = create_connection(self.db_file)
        with conn:
            propdf[header].to_sql('SectionProperty', conn,
                                  index_label=header,
                                  if_exists='append', index=False)
        #
#
def get_sections(conn, component: int):
    """
    """
    query = (component, )
    table = "SELECT * from Section \
             WHERE mesh_id = ?;"
    cur = conn.cursor()
    cur.execute(table, query)
    row = cur.fetchall()
    return row
#
#


