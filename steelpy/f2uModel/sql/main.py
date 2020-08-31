#
# Copyright (c) 2009-2020 fem2ufo
#

# Python stdlib imports
import os
from itertools import chain
from collections import Counter
#
#
# package imports
from steelpy.f2uModel.sql.operation.process_sql import create_connection
from steelpy.f2uModel.sql.write_model.f2u_model import populate_component
from steelpy.f2uModel.sql.write_model.concept import populate_concepts
from steelpy.f2uModel.sql.write_model.materials import populate_materials
from steelpy.f2uModel.sql.write_model.sections import populate_sections
from steelpy.f2uModel.sql.write_model.nodes import populate_nodes, populate_boundaries
from steelpy.f2uModel.sql.write_model.elements import populate_elements
from steelpy.f2uModel.sql.write_model.basic_load import populate_basic_loading
from steelpy.f2uModel.sql.write_model.load_comb import populate_load_combination
from steelpy.f2uModel.sql.write_model.results import populate_displacements, populate_reacctions, populate_force
from steelpy.f2uModel.sql.write_model.postprocess import populate_stress
#
from steelpy.f2uModel.sql.read_model.nodes import get_nodes, get_boundaries
from steelpy.f2uModel.sql.read_model.materials import get_materials
from steelpy.f2uModel.sql.read_model.sections import get_sections
from steelpy.f2uModel.sql.read_model.elements import ElementSQL, get_connectivities
from steelpy.f2uModel.sql.read_model.basic_load import get_basic_load
from steelpy.f2uModel.sql.read_model.results import get_node_displacements, get_element_forces
#from steelpy.f2uModel.sql.read_model.results import 
#
#
# --------------------------
#
#
#
class f2uDB:
    
    __slots__ = ["component_name", "bd_file"]
    
    def __init__(self, component_name:str):
        """
        """
        self.component_name:str = str(component_name)
        self.bd_file = self.component_name + "_f2uDB.db"
        # remove file if exist
        try:
            os.remove(self.bd_file)
        except FileNotFoundError:
            pass
        #
        conn = create_connection(self.bd_file)
        with conn:
            populate_component(conn, component_name)
    #
    #
    def write_concept(self, concept):
        """write out concept data"""
        self.concepts = concept.beam
        #print('-->')

    #
    def write_geometry(self, geometry):
        """
        """
        # make connection with database
        #conn = create_connection(self.bd_file)
        # nodes
        self.nodes = geometry.nodes
        self.boundaries = geometry.boundaries
        # elements
        self.elements = geometry.elements
    #   
    def write_load(self, load):
        """ """
        #self._load = load
        # loading
        self.basic_load = load._basic
        self.load_combination = load._combination
    #
    def write_results(self, results):
        """ """
        #
        # self.displacements
        # self.reactions
        # self.element_forces
        # self.element_stress
        #
        print('-->')
    #
    #
    @property
    def concepts(self):
        """ """
        self._concepts = ConceptsSQL(self.bd_file)
        return self._concepts
    
    @concepts.setter
    def concepts(self, value):
        """ """
        conn = create_connection(self.bd_file)
        with conn:
            populate_concepts(conn, value)
    #
    #
    @property
    def nodes(self):
        """
        """
        conn = create_connection(self.bd_file)
        nodes = get_nodes(conn, self.component_name)
        conn.close()
        return nodes
    
    @nodes.setter
    def nodes(self, value):
        """
        """
        conn = create_connection(self.bd_file)
        with conn:
            populate_nodes(conn, value)
        #conn.commit()
        #conn.close()
        #print('--')
    #
    @property
    def boundaries(self):
        """
        """
        conn = create_connection(self.bd_file)
        boundaries = get_boundaries(conn, self.component_name)
        conn.close()
        return boundaries
    
    @boundaries.setter
    def boundaries(self, value):
        """
        """
        conn = create_connection(self.bd_file)
        with conn:
            populate_boundaries(conn, value)         
    #
    #
    @property
    def materials(self):
        """
        """
        conn = create_connection(self.bd_file)
        materials = get_materials(conn, self.component_name)
        conn.close()
        return materials
    
    @materials.setter
    def materials(self, value):
        """
        """
        conn = create_connection(self.bd_file)
        with conn:
            populate_materials(conn, value)        
    #
    #
    @property
    def sections(self):
        """
        """
        conn = create_connection(self.bd_file)
        sections = get_sections(conn, self.component_name)
        conn.close()
        return sections
    
    @sections.setter
    def sections(self, value):
        """
        """
        conn = create_connection(self.bd_file)
        with conn:
            populate_sections(conn, value)        
    #
    #
    @property
    def elements(self):
        """
        """
        #_elements = ElementSQL(self.bd_file)
        return ElementSQL(self.bd_file)
    
    @elements.setter
    def elements(self, value):
        """
        """
        conn = create_connection(self.bd_file)
        with conn:
            populate_elements(conn, value)
    #
    #def get_element(self, item_name):
    #    """
    #    """
    #    return get_element(self.bd_file, item_name)
    #
    @property
    def connectivities(self):
        """
        """
        conn = create_connection(self.bd_file)
        connectivities = get_connectivities(conn)
        conn.close()
        return connectivities
    #
    #
    @property
    def basic_load(self):
        """
        """
        conn = create_connection(self.bd_file)
        with conn:
            load = get_basic_load(conn, self.component_name)
        #conn.close()
        return load
    
    @basic_load.setter
    def basic_load(self, value):
        """
        """
        conn = create_connection(self.bd_file)
        with conn:
            populate_basic_loading(conn, value)
    #
    #
    @property
    def load_combination(self):
        """ """
        populate_load_combination
    
    @load_combination.setter
    def load_combination(self, value):
        """ """
        conn = create_connection(self.bd_file)
        with conn:
            populate_load_combination(conn, value)       
    #
    #
    @property
    def displacements(self):
        """
        """
        conn = create_connection(self.bd_file)
        disp = get_node_displacements(conn)
        conn.close()
        return disp
    #
    @displacements.setter
    def displacements(self, value):
        """
        """
        conn = create_connection(self.bd_file)
        with conn:
            populate_displacements(conn, self.nodes, value)
    #
    #
    @property
    def reactions(self):
        """
        """
        pass
        #conn = create_connection(self.bd_file)
        #disp = get_node_displacements(conn)
        #conn.close()
        #return disp
    
    @reactions.setter
    def reactions(self, value):
        """
        """
        conn = create_connection(self.bd_file)
        with conn:
            populate_reacctions(conn, self.nodes, value)  
    #
    @property
    def element_forces(self):
        """
        """    
        conn = create_connection(self.bd_file)
        forces = get_element_forces(conn)
        conn.close()
        return forces
    
    @element_forces.setter
    def element_forces(self, value):
        """
        """
        conn = create_connection(self.bd_file)
        with conn:
            populate_force(conn, value)       
        #print("-->")    
    #
    #
    @property
    def element_stress(self):
        """
        """
        pass
    
    @element_stress.setter
    def element_stress(self, value):
        """
        """
        conn = create_connection(self.bd_file)
        with conn:
            populate_stress(conn, value)       
        print("-->")  
    #
    @property
    def get_free_nodes(self):
        """
        find nodes not sharing elements
        """       
        connectivities = self.connectivities
        connectivities = [conn for conn in connectivities.values()]
        #column
        flat = list(chain.from_iterable(connectivities))
        return [k for k, v in Counter(flat).items() if v == 1]      