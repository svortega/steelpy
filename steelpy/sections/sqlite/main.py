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
            self._create_table()
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
    def _create_table(self) -> None:
        """ """
        conn = create_connection(self.db_file)
        #
        table = "CREATE TABLE IF NOT EXISTS Section (\
                    number INTEGER PRIMARY KEY NOT NULL,\
                    name NOT NULL,\
                    component_id INTEGER NOT NULL REFERENCES Component(number),\
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

#
#
#
def get_shapeSQL(shape_type, geometry):
    """ """
    if re.match(r"\b(tub(ular)?|pipe)\b", shape_type, re.IGNORECASE):
        shape = TubularQSL()
        return TubularQSL(name=geometry[0], 
                            diameter=geometry[3], thickness=geometry[4]) 

    #elif re.match(r"\b((solid|bar(\_)?)?rectangle|trapeziod|circular|round)\b", shape_type, re.IGNORECASE):
    #    return self._solid[shape_name]
    elif re.match(r"\b((solid|bar(\_)?)?circular|round)\b", shape_type, re.IGNORECASE):
        d = geometry[3]
        return CircleBasic(name=geometry[0], d=d, type=shape_type)

    elif re.match(r"\b((solid|bar(\_)?)?square|rectangle)\b", shape_type, re.IGNORECASE):
        d = geometry[3]
        wb = geometry[7]
        return RectangleBasic(name=geometry[0], depth=d, width=wb,
                              type=shape_type)

    elif re.match(r"\b((solid|bar(\_)?)?trapeziod)\b", shape_type, re.IGNORECASE):
        d = geometry[5]
        wb = geometry[7]
        wt = geometry[9]            
        c = abs(wt - wb) / 2.0
        return Trapeziod(name=geometry[0], depth=d, width=wb,
                         a=wt, c=c, type=shape_type)    
    
    elif re.match(r"\b(i((\_)?beam|section)?|w|m|s|hp|ub|uc|he|ipe|pg)\b", shape_type, re.IGNORECASE):
        return IbeamBasic(name=geometry[0], 
                          d=geometry[5], tw=geometry[6],
                          bft=geometry[7], tft=geometry[8],
                          bfb=geometry[9], tfb=geometry[10])
    
    elif re.match(r"\b(b(ox)?|rhs|shs)\b", shape_type, re.IGNORECASE):
        return BoxBasic(name=geometry[0], 
                d=geometry[5], tw=geometry[6],
                b=geometry[7], tb=geometry[8])
    
    elif re.match(r"\b(c(hannel)?)\b", shape_type, re.IGNORECASE):
        return ChannelBasic(name=geometry[0], 
                    d=geometry[5], tw=geometry[6],
                    b=geometry[7], tb=geometry[8])
    
    elif re.match(r"\b(t(ee)?)\b", shape_type, re.IGNORECASE):
        return TeeBasic(name=geometry[0], 
                        d=geometry[5], tw=geometry[6],
                        b=geometry[7], tb=geometry[8])
    
    elif re.match(r"\b(l|angle)\b", shape_type, re.IGNORECASE):
        return AngleBasic(name=geometry[0], 
                          d=geometry[5], tw=geometry[6],
                          b=geometry[7], r=0)
    
    else:
        raise IOError(f' Section type {shape_type} not recognised')


