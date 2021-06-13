# 
# Copyright (c) 2019-2021 fem2ufo
# 

# Python stdlib imports
from collections.abc import Mapping
import re
from typing import NamedTuple, Dict, List, Tuple, Union


# package imports
from steelpy.f2uModel.sections.shapes.tubular import TubularInMemory, TubularSQLite
from steelpy.f2uModel.sections.shapes.tee import TeeSQLite, TeeInMemory
from steelpy.f2uModel.sections.shapes.channel import ChannelSQLite, ChannelInMemory
from steelpy.f2uModel.sections.shapes.box import BoxSQLite, BoxInMemory
from steelpy.f2uModel.sections.shapes.ibeam import IbeamInMemory, IbeamSQLite
from steelpy.f2uModel.sections.shapes.solid import CircleInMemory, CircleSQLite, RectangleSQLite, RectangleInMemory
from steelpy.f2uModel.results.sqlite.operation.process_sql import create_connection, create_table
from steelpy.f2uModel.sections.process.io_sections import PropertyOut, get_sect_properties



# ---------------------------------
#
class SectionSQL(Mapping):
    __slots__ = ['db_file', '_sections', '_default']

    def __init__(self, db_file: str,
                 db_system:str="sqlite"):
        """
        """
        self.db_file = db_file
        self._sections : Dict = {}
        self._default: Union[ str, None ] = None
        # create node table
        self._create_table()
    #
    def __setitem__(self, shape_name:Union[str, int], 
                    properties: Union[List[float], Dict[str, float], str]) -> None:
        """
        """
        try:
            self._sections[shape_name]
            raise IOError( '   Section {:} already exist'.format(shape_name))
        except KeyError:
            if isinstance (properties, str):
                shape_type = properties
                #properties = []
                properties = [None, None,
                              None, None,
                              None, None,
                              None, None,
                              None, None]
            else:
                shape_type = properties[0]
                properties = get_sect_properties(properties[1:])
            #
            if re.match(r"\b(i((\_)?beam|section)?|w|m|s|hp|ub|uc|he|ipe)\b"
                        ,shape_type, re.IGNORECASE):
                self._sections[shape_name] = IbeamSQLite(name=shape_name, db_file=self.db_file,
                                                         d=properties[0], tw=properties[1],
                                                         bf=properties[2], tf=properties[3],
                                                         bfb=properties[4], tfb=properties[5])
            
            elif re.match (r"\b(t(ee)?)\b", shape_type, re.IGNORECASE):
                self._sections[ shape_name ] = TeeSQLite(name=shape_name, db_file=self.db_file,
                                                         d=properties[0], tw=properties[1],
                                                         b=properties[2], tb=properties[3])
                
            elif re.match (r"\b(tub(ular)?|pipe|chs)\b", shape_type, re.IGNORECASE):
                self._sections[shape_name] = TubularSQLite(name=shape_name, db_file=self.db_file,
                                                           d=properties[0], t=properties[1])
            
            #elif re.match(r"\b((solid|bar(\_)?)?rectangle|trapeziod)\b", shape_type, re.IGNORECASE):
            elif re.match(r"\b((solid|bar(\_)?)?rectangle)\b", shape_type, re.IGNORECASE):
                self._sections[shape_name] = RectangleSQLite(name=shape_name, db_file=self.db_file,
                                                             d=properties[0], w=properties[1])
                #self._sections[shape_name] = Trapeziod(cls=self)
            
            elif re.match(r"\b((solid|bar(\_)?)?circular|round)\b", shape_type, re.IGNORECASE):
                self._sections[shape_name] = CircleSQLite(name=shape_name, db_file=self.db_file,
                                                           d=properties[0])

            elif re.match ( r"\b(b(ox)?|rhs|shs)\b", shape_type, re.IGNORECASE ):
                self._sections[ shape_name ] = BoxSQLite(name=shape_name, db_file=self.db_file,
                                                         d=properties[0], tw=properties[1],
                                                         b=properties[2], tb=properties[3])

            elif re.match ( r"\b(c(hannel)?)\b", shape_type, re.IGNORECASE ):
                self._sections[ shape_name ] = ChannelSQLite(name=shape_name, db_file=self.db_file,
                                                           d=properties[0], tw=properties[1],
                                                           b=properties[2], tb=properties[3])

            else:
                raise Exception(" section item {:} not recognized".format(shape_type))
            #
            #self._sections[shape_name].name = shape_name
            #
            #if properties:
            #    self._sections[shape_name].geometry(*properties)
            #    conn = create_connection(self.db_file)
            #    with conn:
            #        sect_number = self._push_section_table(conn, self._sections[shape_name])
            #        self._sections[shape_name].number = sect_number
            #        self._push_property_table(conn, self._sections[shape_name])
            #        conn.commit()
            #        #self._number.append(sect_number)
            #
            self._sections[shape_name].push_property()
            #
            #conn = create_connection(self.db_file)
            #with conn:
            #    self._push_property_table(conn, self._sections[shape_name])
            #print('-->')
    #
    def __getitem__(self, shape_name):
        """
        """
        if shape_name in self._sections:
            return self._sections[shape_name]
        else:
            raise Exception(" section name {:} not found".format(shape_name))
    #
    #def push_sections(self):
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
    #def _push_property_table(self, conn, section):
    #    """ """
    #    project = (section.number, *section.properties)
    #    sql = 'INSERT INTO  tb_SecProperties(number, area, Zc, Yc,\
    #                                        Iy, Zey, Zpy, ry,\
    #                                        Iz, Zez, Zpz, rz,\
    #                                        J, Cw)\
    #                                        VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?)'
    #    cur = conn.cursor()
    #    cur.execute(sql, project)
    #
    #def _push_section_table(self, conn, section) -> int:
    #    """
    #    """
    #    project = section._get_section_table()
    #    sql = 'INSERT INTO  tb_Sections(name, title, type, diameter, wall_thickess,\
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
        table_sections = "CREATE TABLE IF NOT EXISTS tb_Sections (\
                            number INTEGER PRIMARY KEY NOT NULL,\
                            name INTEGER NOT NULL,\
                            title TEXT,\
                            type TEXT NOT NULL,\
                            diameter DECIMAL,\
                            wall_thickess DECIMAL,\
                            height DECIMAL,\
                            web_thickness DECIMAL,\
                            top_flange_width DECIMAL,\
                            top_flange_thickness DECIMAL,\
                            bottom_flange_width DECIMAL,\
                            bottom_flange_thickness DECIMAL,\
                            SA_inplane DECIMAL, \
                            SA_outplane DECIMAL,\
                            shear_stress TEXT, \
                            build TEXT, \
                            compactness TEXT);"
        #
        table_properties = "CREATE TABLE IF NOT EXISTS tb_SecProperties (\
                            number INTEGER PRIMARY KEY NOT NULL,\
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
        conn = create_connection(self.db_file)
        create_table(conn, table_sections)
        create_table(conn, table_properties)
    #
    #
    def __len__(self):
        return len(self._sections)
    
    def __iter__(self):
        """
        """
        return iter(self._sections)    
    
    def __contains__(self, value):
        return value in self._sections
    #
   #def get_properties(self):
   #    """
   #    """
   #    summary = {}
   #    for key, item in self._sections.items():
   #        summary[key] = item._get_properties()
   #    return summary
#
#
def get_sections(conn, component_name):
    """
    """
    sections = {}
    cur = conn.cursor()
    cur.execute("SELECT tb_Sections.name, tb_Sections.number, tb_SecProperties.*\
                from tb_Sections, tb_SecProperties\
                WHERE  tb_SecProperties.number = tb_Sections.number;")
    rows = cur.fetchall()
    for row in rows:
        #data = [row[0], *row[4:]]
        #print(row)
        sections[row[0]] = PropertyOut._make(row[3:])
    #conn.close()
    #print("--->")
    return sections
#
def get_sectionSQL(conn, section_number:int):
    """
    """
    cur = conn.cursor()
    #cur.execute("SELECT tb_Sections.name, tb_Sections.number, tb_SecProperties.*\
    #            from tb_Sections, tb_SecProperties \
    #            WHERE tb_Sections.number = {:} \
    #            AND tb_SecProperties.number = tb_Sections.number;".format(section_number))
    #
    cur.execute("SELECT * from tb_SecProperties \
                WHERE tb_SecProperties.number = {:};".format(section_number))    
    row = cur.fetchone()
    sections = PropertyOut(*row[1:])
    #conn.close()
    #print("--->")
    return sections
#
# ---------------------------------
#
class SectionInMemory(Mapping):
    
    __slots__ = ['_sections', '_default']
    
    def __init__(self):
        """
        alpha : Thermal expansion ratio (in 1/K at 20C )
        """
        #  ----- Section -----
        self._sections : Dict = {}
        self._default :Union[str,None] = None
    #
    #
    def __setitem__(self, shape_name: Union[str, int],
                    properties: Union[List[float], Dict[str, float], str]) -> None:
        """
        """
        try:
            self._sections[shape_name]
            raise IOError('   error section {:} already exist'.format(shape_name))
        except KeyError:
            if isinstance (properties, str):
                shape_type = properties
                properties = [None, None,
                              None, None,
                              None, None,
                              None, None,
                              None, None]
            else:
                shape_type = properties[0]
                properties = get_sect_properties(properties[1:])            
            #
            if re.match(r"\b(i((\_)?beam|section)?|w|m|s|hp|ub|uc|he|ipe)\b"
                        ,shape_type, re.IGNORECASE):
                self._sections[shape_name] = IbeamInMemory(name=shape_name,
                                                           d=properties[0], tw=properties[1],
                                                           bf=properties[2], tf=properties[3],
                                                           bfb=properties[4], tfb=properties[5])
            
            elif re.match (r"\b(t(ee)?)\b", shape_type, re.IGNORECASE):
                self._sections[ shape_name ] = TeeInMemory(name=shape_name,
                                                           d=properties[0], tw=properties[1],
                                                           b=properties[2], tb=properties[3])
                
            elif re.match (r"\b(tub(ular)?|pipe)\b", shape_type, re.IGNORECASE):
                self._sections[shape_name] = TubularInMemory(name=shape_name,
                                                             d=properties[0], t=properties[1])
            
            #elif re.match(r"\b((solid|bar(\_)?)?rectangle|trapeziod)\b", shape_type, re.IGNORECASE):
            elif re.match(r"\b((solid|bar(\_)?)?rectangle)\b", shape_type, re.IGNORECASE):
                self._sections[shape_name] = RectangleInMemory(name=shape_name,
                                                               d=properties[0], w=properties[1])
                #self._sections[shape_name] = Trapeziod(cls=self)
            
            elif re.match(r"\b((solid|bar(\_)?)?circular|round)\b", shape_type, re.IGNORECASE):
                self._sections[shape_name] = CircleInMemory(name=shape_name, d=properties[0])

            elif re.match ( r"\b(b(ox)?|rhs|shs)\b", shape_type, re.IGNORECASE ):
                self._sections[ shape_name ] = BoxInMemory(name=shape_name,
                                                           d=properties[0], tw=properties[1],
                                                           b=properties[2], tb=properties[3])

            elif re.match ( r"\b(c(hannel)?)\b", shape_type, re.IGNORECASE ):
                self._sections[ shape_name ] = ChannelInMemory(name=shape_name,
                                                               d=properties[0], tw=properties[1],
                                                               b=properties[2], tb=properties[3])

            else:
                raise Exception(" section item {:} not recognized".format(shape_type))
            #print('-->')
            #self._sections[shape_name].name = shape_name
            #self._sections[shape_name].number = len(self._sections)
            #if properties[0]:
            #    self._sections[shape_name].geometry(*properties)
    
    def __getitem__(self, shape_name: Union[str, int]):
        """
        """
        if shape_name in self._sections:
            return self._sections[shape_name]
        else:
            raise Exception(" section name {:} not found".format(shape_name))   
    #
    def get_item_by_number(self, shape_name):
        """
        """
        _items = {_item.number:key for key, _item in self._sections.items()}
        try:
            _name = _items[shape_name]
            return self.__getitem__(_name)
        except KeyError:
            raise KeyError('Invalid section number')
    #
    #
    def __delitem__(self, shape_name: str) -> None:
        """
        """
        #_number = self._sections[shape_name].number
        del self._sections[shape_name]
    
    def __len__(self):
        return len(self._sections)
    
    def __iter__(self):
        """
        """
        return iter(self._sections)    
    
    def __contains__(self, value):
        return value in self._sections
    #
    #def get_properties(self):
    #    """
    #    """
    #    summary = {}
    #    for key, item in self._sections.items():
    #        summary[key] = item._get_properties()
    #    return summary
#
# ---------------------------------
#
class Sections(Mapping):
    
    __slots__ = ['_sections', '_default']
    
    def __init__(self, mesh_type:str,
                 db_file:Union[str,None]):
        """
        """
        if mesh_type != "inmemory":
            self._sections = SectionSQL(db_file=db_file,
                                        db_system=mesh_type)
        else:
            self._sections = SectionInMemory()
    #
    def __setitem__(self, shape_name:Union[str, int], 
                    shape_type:str) -> None:
        """
        """
        self._sections[shape_name] = shape_type
    #
    def __getitem__(self, shape_name:str):
        """
        """
        return self._sections[shape_name]
    #
    #
    def __delitem__(self, shape_name: str) -> None:
        """
        """
        del self._sections[shape_name]
    
    def __len__(self):
        return len(self._sections)
    
    def __iter__(self):
        """
        """
        return iter(self._sections)    
    
    def __contains__(self, value):
        return value in self._sections
    #
    #def __str__(self) -> str:
    #    """ """
    #    output = []
    #    output.append("\n")
    #    output.append("_______________________________________________________________________________________\n")
    #    output.append("\n")
    #    output.append("                              SECTION DERIVED PROPERTIES"+"\n")
    #    output.append("\n")
    #    output.append("Member ID      Area[mm^2]  I   [mm^4]  S   [mm^3]  Z   [mm^3]  ShapeFctor  r    [mm]\n")
    #    #output.append("Number        Awx  [mm^2]  Iyy [mm^4]  Syy [mm^3]  Zyy [mm^3]  SCeny [mm]  ry   [mm]"+"\n")
    #    output.append("               Mass[kg/m]  Ip  [mm^4]  J   [mm^4]"+"\n")
    #    output.append(".......................................................................................\n")
    #    output.append("\n")
    #    return output
    #
    def __str__(self, units:str="si") -> str:
        """ """
        unit_sec = " m"
        unit_mas = "kg/m"
        space = " "
        #
        output = "\n"
        output += "{:}\n".format(80*"_")
        output += "\n"
        output += f"{30*space}SECTION PROPERTIES REPORT{15*space}UNITS [{unit_sec}]\n"
        output += "\n"
        output += f"{48*space}Web{14*space}Flanges\n"
        output += (f"Member Name{5*space}Type{6*space}Diametre{3*space}Thickness"
                   f"{2*space}Height{5*space}Top Width{2*space}Bot Width\n")
        output += f"{48*space}Thickness{2*space}Top Thick{2*space}Bot Thick\n"
        output += f"{70*space}Fillet\n"
        output += "{:}\n".format(80*".")
        output += "\n"
        for name, section in self._sections.items():
            output += "{:<15s} ".format(str(name))
            output += section._dimension()
        #
        output += "\n"
        output += "{:}\n".format(80*"_")
        output += "\n"
        output += "                               SECTION DERIVED PROPERTIES\n"
        output += "{:}\n".format(80*"_")
        output += "\n"
        output += (f"{15*space}Area[{unit_sec}^2] Ixx [{unit_sec}^4] Iyy [{unit_sec}^4]"
                   f" Yp    [{unit_sec}] rx    [{unit_sec}] J   [{unit_sec}^4]\n")
        output += (f"{26*space}Sxx [{unit_sec}^3] Syy [{unit_sec}^3] SCeny [{unit_sec}]"
                   f" ry    [{unit_sec}] Cw  [{unit_sec}^6]\n")
        output += f"{26*space}Zxx [{unit_sec}^3] Zyy [{unit_sec}^3] SCenx [{unit_sec}] Mass[{unit_mas}]\n"
        output += "{:}\n".format(80*".")
        output += "\n"
        for name, section in self._sections.items():
            output += "{:<14s} ".format(str(name))
            output += section.properties.__str__()
        #print("-->")
        return output
    #
    #
    @property
    def default(self):
        """ """
        return self._sections._default
    
    @default.setter
    def default(self, shape_name):
        """ """
        try:
            self._sections[shape_name]
        except KeyError:
            raise IOError(f'section {shape_name} missing')
        self._sections._default = shape_name
    #
    def get_properties(self):
        """
        """
        summary = {}
        for key, item in self._sections.items():
            item.push_property()
            #item.properties
            #summary[key] = item._get_properties()
        #return summary
#    
# ---------------------------------
#