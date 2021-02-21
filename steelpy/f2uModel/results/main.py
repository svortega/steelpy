#
# Copyright (c) 2009-2021 fem2ufo
#

# Python stdlib imports
from typing import NamedTuple, Dict, List, Iterable, Union

#
# package imports
from steelpy.f2uModel.results.sqlite.main import ResultSQL
from steelpy.f2uModel.results.inmemory.main import ResultInmemory

#
#
class Results:
    """
    """
    __slots__ = ['_results']
    
    def __init__(self, basic_loads,
                 mesh_type:str,
                 db_file:Union[str,None]):
        """
        """
        if mesh_type != "inmemory":
            self._results = ResultSQL(db_file=db_file,
                                      db_system=mesh_type)
        else:
            self._results = ResultInmemory(basic_loads=basic_loads)
    #
    #def __setitem__(self, material_name:Union[str, int],
    #                material_type:str) -> None:
    #    """
    #    """
    #    self._results[material_name] = material_type
    ##
    #def __getitem__(self, material_name:str):
    #    """
    #    """
    #    return self._results[material_name]
    #
    @property
    def node(self):
        """ """
        return NodeType(self)
    #
    #
    @property
    def element(self):
        """ """
        return BeamType(self)
#
#
class NodeType:
    __slots__ = ['_cls']
    
    def __init__(self, cls):
        """ """
        self._cls = cls
    #
    @property
    def deflection(self):
        """ """
        return self._cls._results.node_displacement
    
    @deflection.setter
    def deflection(self, value:List):
        """ """
        self._cls._results.node_displacement = value
    #
    def print_deflections(self):
        """ """
        self._cls._results.print_node_displacement()
    #
    def reaction(self):
        """ """
        pass    
#
#
class BeamType:
    __slots__ = ['_cls']
    
    def __init__(self, cls):
        """ """
        self._cls = cls
    #
    #def axial(self):
    #    """ """
    #    pass
    ##
    #def torsion(self):
    #    """ """
    #    pass
    ##
    #def shear(self):
    #    """ """
    #    pass
    ##
    #def bending(self):
    #    """ """
    #    pass
    ##
    @property
    def force(self):
        """ """
        return self._cls._results.element_force
    
    @force.setter
    def force(self, value):
        """ """
        self._cls._results.element_force = value

    def print_forces(self):
        """ """
        self._cls._results.print_element_forces()
#
#