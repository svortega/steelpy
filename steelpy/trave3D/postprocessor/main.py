#
# Copyright (c) 2009-2023 fem2ufo
#
from __future__ import annotations
# Python stdlib imports
import pickle
#from typing import NamedTuple

#
# package imports
# steelpy.trave3D.postprocessor
#from .operations import get_reactions
from .output import ResultInmemory

#
#
class Results:
    """
    """
    __slots__ = ['_results', '_mesh', '_m2D']
    
    def __init__(self, mesh, m2D: bool):
        """
        """
        self._results = ResultInmemory(m2D=m2D)
        self._mesh = mesh
        self._m2D = m2D
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
    # -----------------
    #
    #@property
    #def mesh(self):
    #    """
    #    FE model mesh
    #    """
    #    return self._mesh
    #
    #@property
    #def load(self):
    #    """
    #    FE model load
    #    """
    #    return self._mesh.load()
    #
    #@property
    #def K(self):
    #    """Global Stiffnes matrix"""
    #    file = open ( "stfmx.f2u", "rb" )
    #    jbc = pickle.load( file )
    #    stf = pickle.load( file )
    #    file.close()
    #    return stf
    #
        #
    # -----------------
    #
    def postprocess(self, df_ndisp, df_nforce, df_membf, df_reactions):
        """Postprocess"""
        #
        #self.process_deflection(df_ndisp=df_ndisp, df_nload=df_nload)
        self._results._displacement = df_ndisp
        self._results._node_force = df_nforce
        self._results._reaction = df_reactions
        self._results._beam_force = df_membf
    #
    #
    #
    #
    # -----------------
    #
    @property
    def node(self):
        """ """
        return NodeType(self)
    #
    @property
    def beam(self):
        """ """
        return BeamType(self)
    #
    # -----------------
    #
    def print(self):
        """ """
        self.node.print_deflections()
        self.node.print_reactions()
        self.beam.print_forces()
        self.beam.print_displacement()
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
        return self._cls._results._displacement
    
    #@deflection.setter
    #def deflection(self, df_ndisp):
    #    """ """
    #    self._cls._results._displacement = df_ndisp
    #    self.process_deflection(df_ndisp=df_ndisp)
    #
    #def process_deflection(self, df_ndisp):
    #    """ """
    #    # geometry
    #    mesh = self._cls._mesh
    #    elements = mesh.elements()
    #    df_nload = mesh._load._df_nodal
    #    boundaries = mesh.boundaries()
    #    # loading
    #    load =  self._cls.load
    #    basic_load = load.basic()
    #    # node force
    #    nforce = node_force(elements, basic_load,
    #                        df_ndisp, df_nload)
    #    self._cls._results._node_force = nforce
    #    # node reactions
    #    reactions = get_reactions(boundaries, nforce)
    #    self._cls._results._reaction = reactions
    #
    #
    def print_deflections(self):
        """ """
        self._cls._results.print_node_displacement()
    #
    @property
    def reaction(self):
        """ """
        #return self._cls._results.node_reaction
        #ndisp = self.deflection
        # -----------------------------------
        # geometry
        #mesh = self._cls._mesh
        #elements = mesh.elements()
        #nodes = mesh.nodes()
        #boundaries = mesh.boundaries()
        # -----------------------------------
        #nforce = self._cls._results._node_force
        #
        #reactions = get_reactions(nodes, boundaries, nforce)
        #self._cls._results._reaction = reactions
        return self._cls._results._reaction

    #
    def print_reactions(self):
        """ """
        self._cls._results.print_node_reactions()
#
#
class BeamType:
    __slots__ = ['_cls']
    
    def __init__(self, cls):
        """ """
        self._cls = cls
        #self._get_force()
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
        return self._cls._results.beam_force
    
    #@force.setter
    #def force(self, value):
    #    """ """
    #    self._cls._results.element_force = value
    #
    def print_displacement(self):
        """ """
        self._cls._results.print_element_disp()
    #
    #
    def print_forces(self):
        """ """
        self._cls._results.print_element_forces()
#
#