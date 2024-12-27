# 
# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
from collections.abc import Mapping
#import re
#
# package imports
#
from steelpy.sections.inmemory.main import SectionIM
from steelpy.sections.sqlite.main import SectionSQL
from steelpy.sections.utils.operations import get_section


# ---------------------------------
#
# ---------------------------------
#
class Section(Mapping):
    __slots__ = ['_sections', '_default', '_mesh_id']

    def __init__(self,
                 mesh_id:int,
                 db_file: str | None = None, 
                 mesh_type: str = 'inmemory'):
        """
        """
        self._mesh_id = mesh_id
        #
        if mesh_type != "inmemory":
            self._sections = SectionSQL(db_file=db_file,
                                        mesh_id=mesh_id, 
                                        db_system=mesh_type)
        else:
            self._sections = SectionIM()
    #
    # -----------------------------------------------
    #
    def __setitem__(self, shape_name: str | int,
                    properties: str|list|dict) -> None:
        """
        """
        properties = get_section(properties)
        self._sections[shape_name] = properties  
    
    def __getitem__(self, shape_name: str | int):
        """
        """
        try:
            return self._sections[shape_name]
        except KeyError:
            raise Exception(f" section name {shape_name} not found")    

    #
    # -----------------------------------------------
    #
    #
    def __delitem__(self, shape_name: str | int) -> None:
        del self._sections[shape_name]

    def __len__(self):
        return len(self._sections)

    def __iter__(self):
        return iter(self._sections)

    def __contains__(self, value):
        return value in self._sections

    #     
    def __str__(self, units: str = "si") -> str:
        """ """
        unit_sec = " m"
        unit_mas = "kg/m"
        space = " "
        #
        output = "\n"
        output += "{:}\n".format(80 * "_")
        output += "\n"
        output += f"{30 * space}SECTION PROPERTIES REPORT{15 * space}UNITS [{units}]\n"
        output += "\n"
        output += f"{48 * space}-- WEB --  ------ FLANGE ------\n"
        output += f"{37 * space}Wall{7 * space}{11 * space}-- TOP --{2 * space}-BOTTOM -\n"
        output += (f"Name{12 * space}Type{6 * space}Diameter{3 * space}Thickness"
                   f"{2 * space}Height{5 * space}Width{6 * space}Width\n")
        output += f"{48 * space}Thickness{2 * space}Thickness{2 * space}Thickness\n"
        output += f"{70 * space}Fillet\n"
        output += "{:}\n".format(80 * ".")
        output += "\n"
        for name, section in self._sections.items():
            output += "{:<15s} ".format(str(name))
            output += section._dimension()
        #
        output += "\n"
        output += "{:}\n".format(80 * "_")
        output += "\n"
        output += "                               SECTION DERIVED PROPERTIES\n"
        output += "{:}\n".format(80 * "_")
        output += "\n"
        output += (f"{15 * space}Area[{unit_sec}^2] Ixx [{unit_sec}^4] Iyy [{unit_sec}^4]"
                   f" Yp    [{unit_sec}] rx    [{unit_sec}] J   [{unit_sec}^4]\n")
        output += (f"{26 * space}Sxx [{unit_sec}^3] Syy [{unit_sec}^3] SCeny [{unit_sec}]"
                   f" ry    [{unit_sec}] Cw  [{unit_sec}^6]\n")
        output += f"{26 * space}Zxx [{unit_sec}^3] Zyy [{unit_sec}^3] SCenx [{unit_sec}] Mass[{unit_mas}]\n"
        output += "{:}\n".format(80 * ".")
        output += "\n"
        for name, section in self._sections.items():
            output += "{:<14s} ".format(str(name))
            prop = section.properties(poisson=0)
            output += prop.__str__()
        # print("-->")
        return output

    #
    # -----------------------------------------------
    #
    @property
    def default(self):
        """ """
        return self._sections._default

    @default.setter
    def default(self, shape_name: int|str):
        """ """
        try:
            self._sections[shape_name]
        except KeyError:
            raise IOError(f'section {shape_name} missing')
        self._sections._default = shape_name

    #
    #
    def get_item_by_number(self, shape_name: str | int):
        """
        """
        _items = {_item.number: key
                  for key, _item in self._sections.items()}
        try:
            _name = _items[shape_name]
            return self.__getitem__(_name)
        except KeyError:
            raise KeyError('Invalid section number')
    #     
    #
    # -----------------------------------------------
    #
    def tubular(self, values:None|list=None,
                df=None):
        """Tubular section"""
        1/0
        if values:
            if isinstance(values, list):
                #self._sections.tubular
                print('-->')
                1/0
            else:
                raise IOError('material input not valid')
        #
        # dataframe input
        try:
            df.columns
            df = df.drop_duplicates(['name'], keep='first')
            #
            # shape_type = properties[0]
            # properties = get_sect_properties(properties[1:])
            #
            self._sections._labels.extend(df.name.tolist())
            self._sections._title.extend(df.title.tolist())
            self._sections._number.extend([next(self._sections.get_number())
                                           for _ in df.name])
            self._sections._type.extend(df.type.tolist())
            #
            self._sections._tubular.df = df
        except AttributeError:
            pass

        #print("here")
        return self._sections._tubular
    #
    #
    # -----------------------------------------------
    #
    @property
    def df(self):
        """ raw data for dataframe"""
        return self._sections.df

    @df.setter
    def df(self, df):
        """ """
        self._sections.df = df
#    
# ---------------------------------
#

#
#
#
