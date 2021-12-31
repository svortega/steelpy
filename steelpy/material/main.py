# 
# Copyright (c) 2009-2022 steelpy
# 


# Python stdlib imports
from collections.abc import Mapping
from typing import NamedTuple, Dict, List, Iterable, Union
#

# package imports
import steelpy.material.operations as operations
from steelpy.material.matsql import MaterialSQL
from steelpy.material.matinmemory import MaterialInmemory
#
#
#
class Materials(Mapping):
    """This module stores material classes.
       -------------------------------
       Nonlinear materials:
       - Nonlinear elastic material
       - Bilinear elastoplastic material
         Yield criterion:
          - von Mises
          - Treca
          - Morh-Coulomb
          - Drucker-Prager
         Hardering rules:
          - Isotropic
          - Kinematic
          - Isotropic + Kinematic
       - Multilinear plastic material
       - Rigid-plastic material
       
       -----------------------------
    """
    __slots__ = ['_material', '_default']
    
    def __init__(self, mesh_type:str='inmemory',
                 db_file:Union[str,None]=None) -> None:
        """
        """
        if mesh_type != "inmemory":
            self._material = MaterialSQL(db_file=db_file,
                                         db_system=mesh_type)
        else:
            self._material = MaterialInmemory()
        #
        self._default = None
    #
    
    def __setitem__(self, material_name:Union[str, int], 
                    properties: Union[List[float], Dict[str, float], str]) -> None:
        """
        """
        # get material type
        if isinstance (properties, str):
            material_type = operations.find_material_type(properties) 
            properties = []
        else:
            material_type = operations.find_material_type(properties[0])
            properties = properties[1:]
        #
        # set material default plastic
        if 'curve' == material_type :
            raise Exception('--> Mat type No ready')
        elif 'elastic' == material_type :
            properties = operations.get_linear_mat_properties(properties)
        else:
            raise IOError(' material type {:} not recognised'
                          .format(material_type))            
        #
        self._material[material_name] = [material_type, *properties]
    #
    def __getitem__(self, material_name:str):
        """
        """
        return self._material[material_name]
    #
    @property
    def default(self):
        """ """
        return self._default
    
    @default.setter
    def default(self, material_name):
        """ """
        try:
            self._material[material_name]
        except KeyError:
            raise IOError(f'material {material_name} missing')
            
        self._default = material_name
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
    def __str__(self, units:str="si") -> str:
        """ """
        stress = "N/mm2"
        density = "kg/m3"
        space = " "
        #
        output = "\n"
        output += "{:}\n".format(80*"_")
        output += "\n"
        output += f"{33*space}MATERIAL PROPERTIES\n"
        output += "\n"
        output += (f"Member ID      Fy [{stress}] Fu [{stress}] E  [{stress}] "
                   f"G  [{stress}] Poisson    Rho[{density}]\n")
        output += "\n"
        output += "{:}\n".format(80*".")
        output += "\n"
        for name, mat in self._material.items():
            output += "{:<14s} ".format(str(name))
            output += mat.__str__()         
            #output += "{:<14s} {:1.4E} {:1.4E} {:1.4E} {:1.4E} {:1.4E} {:1.4E}\n"\
            #    .format(str(name), mat.Fy.value, mat.Fu.value, mat.E.value,
            #            mat.G.value, mat.poisson, mat.density.value)
        return output
#
#