# 
# Copyright (c) 2009 steelpy
#

# Python stdlib imports
from __future__ import annotations
from collections.abc import Mapping
import re
#

# package imports
from .utils.operations import get_mat_properties, get_mat_df
from .sqlite.isotropic import MaterialSQL
from .inmemory.isotropic import MaterialIM


#
#
#
class Material(Mapping):
    """This module stores material classes.
       -------------------------------
       Linear material:
       - Isotropic
       
       -----------------------------
    """
    __slots__ = ['_material']

    def __init__(self,
                 component:int,
                 db_file: str | None = None, 
                 mesh_type: str = 'inmemory') -> None:
        """
        """
        if mesh_type != "inmemory":
            self._material = MaterialSQL(db_file=db_file,
                                         component=component, 
                                         db_system=mesh_type)
        else:
            self._material = MaterialIM()

    #
    def __setitem__(self, material_name: str | int,
                    properties: tuple|list|dict[str, float]|str) -> None:
        """
        [elastic, Fy, Fu, E, G, Poisson, density, alpha, title(optional)]
        """
        # get material type
        properties = get_mat_properties(properties)
        self._material[material_name] = properties

    #
    def __getitem__(self, material_name: str):
        """
        """
        return self._material[material_name]

    #
    @property
    def default(self):
        """ """
        return self._material.default

    @default.setter
    def default(self, material_name):
        """ """
        try:
            self._material[material_name]
        except KeyError:
            raise IOError(f'material {material_name} missing')

        self._material.default = material_name

    #
    def __len__(self) -> float:
        return len(self._material)

    def __iter__(self):
        """
        """
        return iter(self._material)

    def __contains__(self, value) -> bool:
        return value in self._material

    #
    def __str__(self, units: str = "si") -> str:
        """ """
        stress = "N/mm2"
        density = "kg/m3"
        space = " "
        #
        output = "\n"
        output += "{:}\n".format(80 * "_")
        output += "\n"
        output += f"{33 * space}MATERIAL PROPERTIES\n"
        output += "\n"
        output += (f"Member ID      Fy [{stress}] Fu [{stress}] E  [{stress}] "\
                   f"G  [{stress}] Poisson    Rho[{density}]\n")
        output += "\n"
        output += "{:}\n".format(80 * ".")
        output += "\n"
        for name, mat in self._material.items():
            output += "{:<14s} ".format(str(name))
            output += mat.__str__()
            # output += "{:<14s} {:1.4E} {:1.4E} {:1.4E} {:1.4E} {:1.4E} {:1.4E}\n"\
            #    .format(str(name), mat.Fy.value, mat.Fu.value, mat.E.value,
            #            mat.G.value, mat.poisson, mat.density.value)
        #
        return output

    #
    #
    @property
    def elastic(self):
        """
        Linear elastic material
        """
        return self._material.elastic
    
    @property
    def linear(self):
        """
        Linear elastic material
        """
        return self._material.elastic
    #
    @property
    def df(self):
        """Data frame format"""
        return self._material._elastic.df
    #
    @df.setter
    def df(self, df):
        """ """
        try:
            columns = list(df.columns)
            df = get_mat_df(df)
            group = df.groupby("type")
            # Elastic type
            try:
                elastic = group.get_group("elastic")
                self._material._elastic.df = elastic
            except KeyError:
                # nonlin = group.get_group("plastic")
                raise IOError('Material type not valid')
        except AttributeError:
            raise IOError('Material type not valid')
        #print('--')
#
#
