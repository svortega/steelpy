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
from steelpy.f2uModel.sql.dump_model.materials import populate_materials
from steelpy.f2uModel.sql.dump_model.sections import populate_sections
from steelpy.f2uModel.sql.dump_model.nodes import populate_nodes, populate_boundaries
from steelpy.f2uModel.sql.dump_model.elements import populate_elements
from steelpy.f2uModel.sql.dump_model.basic_load import populate_basic_loading
from steelpy.f2uModel.sql.dump_model.load_comb import populate_load_combination
from steelpy.f2uModel.sql.dump_model.results import populate_displacements, populate_reacctions, populate_force
from steelpy.f2uModel.sql.dump_model.postprocess import populate_stress
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
#
def populate_sql(component_name, component):
    """
    """
    #
    #model = component[component_name]
    model = component
    #
    try:
        #os.remove(r'C:\Users\svort\Dropbox\Python\fem2ufo\f2u\V055\examples\jacket\f2uDB.db')
        #os.remove(r'/home/chava/Dropbox/Python/steelpy/steelpy_03/f2uDB.db')
        os.remove('f2uDB.db')
    except FileNotFoundError:
        pass    
    #
    #conn = create_connection(r'C:\Users\svort\Dropbox\Python\fem2ufo\f2u\V055\examples\jacket\f2uDB.db')
    #conn = create_connection(r'/home/chava/Dropbox/Python/steelpy/steelpy_03/f2uDB.db')
    conn = create_connection('f2uDB.db')
    #
    #for _rel in model.mesh.releases:
    #    print(_rel)
    #
    #if conn is not None:
    #    create_table(conn, _table_nodes)
    #else:
    #    raise IOError("Error! cannot create the database connection.")
    #
    # loading
    #populate_loading(conn, model.load)    
    #
    #
    # materials
    populate_materials(conn, model.materials)
    # sections
    populate_sections(conn, model.sections)    
    # nodes
    populate_nodes(conn, model.nodes, model.boundaries)
    # elements
    populate_elements(conn, model.elements)
    #
    #
    #
    conn.commit()
    conn.close()
#
def populate_loading_sql(component_name, component):
    """
    """
    #
    #conn = create_connection(r'D:\svort\Dropbox\Python\fem2ufo\f2u\V055\examples\jacket\f2uDB.db')
    #conn = create_connection(r'/home/chava/Dropbox/Python/steelpy/steelpy_03/f2uDB.db')
    conn = create_connection('f2uDB.db')
    #
    # loading
    #for component_name, _component in model.load.items():
    populate_basic_loading(conn, component)
    #
    conn.commit()
    conn.close()    
#
def populate_component_sql(model):
    """
    """
    #
    try:
        #os.remove(r'D:\svort\Dropbox\Python\fem2ufo\f2u\V055\examples\jacket\f2uDB.db')
        os.remove(r'/home/chava/Dropbox/Python/steelpy/steelpy_03/f2uDB.db')
    except FileNotFoundError:
        pass
    #
    #conn = create_connection(r'D:\svort\Dropbox\Python\fem2ufo\f2u\V055\examples\jacket\f2uDB.db')
    conn = create_connection(r'/home/chava/Dropbox/Python/steelpy/steelpy_03/f2uDB.db')
    #
    
    #for _rel in model.mesh.releases:
    #    print(_rel)
    #
    for component_name, _component in model.components.items():
        # materials
        populate_materials(conn, _component.materials, component_name)
        # sections
        populate_sections(conn, _component.sections, component_name)   
        # nodes
        populate_nodes(conn, _component.mesh, component_name)
        # elements
        populate_members(conn, _component.mesh, component_name)
    #
    #
    #
    conn.commit()
    conn.close()    
#
def populate_soil_sql(model):
    """
    """
    #conn = create_connection(r'C:\Users\svort\Dropbox\Python\fem2ufo\f2u\V040\examples\jacket\f2uDB.db')
    conn = create_connection(r'/home/chava/Dropbox/Python/steelpy/steelpy_03/f2uDB.db')
    #
    for component_name, _component in model.foundation.items():
        # materials
        populate_materials(conn, _component.materials, component_name)
    #
    conn.commit()
    conn.close()    
#
# --------------------------
#
#
#
class f2uDB:
    
    __slots__ = ["component_name", "bd_file", 
                 "_model", "_elements", "_load"]
    
    def __init__(self, model, component_name:str="f2uDB"):
        """
        """
        self.component_name:str = str(component_name)
        self.bd_file = self.component_name + ".db"
        self._model = model.mesh
        self._load = model._load
        #self._conn = create_connection(self.bd_file)
    #
    def read_model(self):
        """
        """
        #conn = create_connection(self.bd_file)
        nodes = self.nodes
        boundaries = self.boundaries
        #
        materials = self.materials
        sections = self.sections
        elements = self.elements
        #
        basic_load = self.basic_load
        #
        #conn.close()
        print("--->")
    #
    def dump_model(self):
        """
        """
        # remove file if exist
        try:
            os.remove(self.bd_file)
        #except PermissionError:
        #    
        except FileNotFoundError:
            pass
        # make connection with database
        #conn = create_connection(self.bd_file)
        # nodes
        self.nodes = self._model.nodes
        self.boundaries = self._model.boundaries        
        # materials
        self.materials = self._model.materials
        # sections
        self.sections = self._model.sections 
        # elements
        self.elements = self._model.elements
        # loading
        self.basic_load = self._load._basic
        self.load_combination = self._load._combination
        #
        #self.element_stress = None
        #
        # close
        #conn.commit()
        #conn.close()        
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
        self._elements = ElementSQL(self.bd_file)
        return self._elements
    
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