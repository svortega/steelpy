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
from .operations import (get_reactions, beam_end_force, beam_force,
                         updape_ndf, updape_memberdf, updape_memberdf2)
from .output import ResultInmemory

#
#
class Results:
    """
    """
    __slots__ = ['_results', '_mesh']
    
    def __init__(self, mesh):
        """
        """
        self._results = ResultInmemory()
        self._mesh = mesh
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
    @property
    def mesh(self):
        """
        FE model mesh
        """
        return self._mesh
    #
    @property
    def load(self):
        """
        FE model load
        """
        return self._mesh.load()
    #
    @property
    def K(self):
        """Global Stiffnes matrix"""
        file = open ( "stfmx.f2u", "rb" )
        jbc = pickle.load( file )
        stf = pickle.load( file )
        file.close()
        return stf
    #
        #
    # -----------------
    #
    def postprocess(self, df_ndisp):
        """Postprocess"""
                # geometry
        mesh = self._mesh
        elements = mesh.elements()
        boundaries = mesh.boundaries()
        #
        # loading
        load =  self._mesh.load()
        basic_load = load.basic()
        load_combination = load.combination()
        df_comb = load_combination.to_basic()
        #
        #
        # TODO: node load should be in SQL
        df_nload = load._df_nodal
        #
        # get beam end node forces
        df_nforce = beam_end_force(elements, df_ndisp)
        #
        # -----------------------------------
        # get beam force along lenght 
        df_membf = beam_force(elements, 
                              basic_load,
                              df_ndisp=df_ndisp, 
                              df_nforce=df_nforce)
        #
        #
        # ---------------------------------------
        # Update df to include load combinations
        # ---------------------------------------
        #
        df_ndisp = updape_ndf(dfnode=df_ndisp, dfcomb=df_comb,
                              values=['x', 'y', 'z', 'rx', 'ry', 'rz'])
        #
        df_nload = updape_ndf(dfnode=df_nload, dfcomb=df_comb,
                              values=['Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz'])
        #
        df_nforce = updape_memberdf(dfmemb=df_nforce, dfcomb=df_comb,
                                    values=['Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz'])
        #
        df_memb = updape_memberdf2(dfmemb=df_membf, dfcomb=df_comb,
                                    values=['F_Vx', 'F_Vy', 'F_Vz', 'F_Mx', 'F_My', 'F_Mz',
                                            'F_wx', 'F_wy', 'F_wz', 'F_phix', 'F_thetay', 'F_thetaz'])
        #
        # -------------------------------------
        # node reactions
        # -------------------------------------
        #
        df_reactions = get_reactions(boundaries, df_nforce)
        #
        #self.process_deflection(df_ndisp=df_ndisp, df_nload=df_nload)
        self._results._displacement = df_ndisp
        self._results._node_force = df_nforce
        self._results._reaction = df_reactions
        self._results._beam_force = df_memb
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
        #return self._cls._results.element_force
        return self._cls._results.beam_force
    
    #@force.setter
    #def force(self, value):
    #    """ """
    #    self._cls._results.element_force = value
    #
    #def _get_force(self):
    #    """ """
    #    ndisp = self._cls._results._displacement
    #    nforce = self._cls._results._node_force
    #    # -----------------------------------
    #    # geometry
    #    mesh = self._cls.mesh
    #    elements = mesh.elements()
    #    #nodes = mesh.nodes()
    #    #boundaries = mesh.boundaries()
    #    # loading
    #    load =  self._cls.load
    #    basic_load = load.basic()
    #    #nload = load._df_nodal
    #    #
    #    load_combination = load.combination()
    #    # -----------------------------------
    #    # beam 
    #    membf = beam_force(elements, 
    #                       basic_load,
    #                       load_combination=load_combination,
    #                       #node_load=nload,
    #                       node_disp=ndisp, 
    #                       node_force=nforce)
    #    self._cls._results._beam_force = membf
    #
    #
    def print_forces(self):
        """ """
        self._cls._results.print_element_forces()
#
#