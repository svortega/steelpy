# 
# Copyright (c) 2009-2023 fem2ufo
#
# Python stdlib imports
from __future__ import annotations
import pickle
#

# package imports
# steelpy.f2uModel
from ..load.main import Load
# steelpy.f2uModel.mesh
from .inmemory.sets import Groups
from .nodes import Nodes
from .boundaries import Boundaries
from .elements import Elements
#
from steelpy.process.dataframe.main import DBframework
from .process.Kmatrix.assemble import assemble_Kmatrix_np, assemble_banded_Kmatrix
#

#
class Mesh:
    """
    mesh[beam_name] = [number, element1, element2, elementn]
    """
    __slots__ = ['_nodes', '_elements', '_load', 'data_type',
                 '_eccentricities', '_boundaries', '_groups',
                 'db_file', '_df', 
                 '_Kmatrix']

    def __init__(self, materials, sections,
                 mesh_type:str="sqlite",
                 db_file:str|None = None):
        """
        """
        #mesh_type2 = 'sqlite'
        self.db_file = db_file
        self.data_type = mesh_type
        self._nodes = Nodes(mesh_type=mesh_type,
                            db_file=self.db_file)
        self._boundaries = Boundaries(mesh_type=mesh_type,
                                      db_file=self.db_file)
        #
        self._elements = Elements(nodes=self._nodes,
                                  materials=materials,
                                  sections=sections,
                                  mesh_type=mesh_type,
                                  db_file=self.db_file)
        # groups
        self._groups = Groups()
        #
        self._load = Load(nodes=self._nodes,
                          elements=self._elements,
                          mesh_type=mesh_type,
                          db_file=self.db_file)
        # Ops
        self._df = DBframework()
        self._Kmatrix = False
    #
    def nodes(self, values:None|list=None,
              df=None):
        """
        """
        supports = self._boundaries.supports()
        #if values:
        if isinstance(values, list):
            for value in values:
                self._nodes[value[0]] = value[1:4]
                try:
                    #if isinstance(value[4], str):
                    #    supports[value[0]] = value[4]
                    #else:
                    supports[value[0]] = value[4:]
                except IndexError:
                    pass
        #
        # dataframe input
        try:
            df.columns   
            self._nodes.df(df)
        except AttributeError:
            pass
        #
        return self._nodes

    #
    def elements(self, values:None|list=None,
                 df=None):
        """
        """
        if isinstance(values, list):
        #if values:
            for value in values:
                self._elements[value[0]] = value[1:]
        #
        #
        # dataframe input
        try:
            df.columns   
            self._elements.df(df)
        except AttributeError:
            pass
        return self._elements
    #
    def boundaries(self, values:None|list=None,
                   df=None):
        """
        """
        #if values:
        if isinstance(values, list):
            for item in values:
                self._boundaries[item[0]] = item[1:]
        #
        # dataframe input
        try:
            df.columns   
            self._boundaries.df(df)
        except AttributeError:
            pass 
        #
        return self._boundaries
    #
    def load(self, values:None|list=None,
             df=None):
        """
        """
        if isinstance(values, list):
            for item in values:
                item
        #
        # dataframe input
        try:
            df.columns
            #self._boundaries.df(df)
        except AttributeError:
            pass
        #
        return self._load
    #
    @property
    def groups(self):
        """
        """
        return self._groups
    #
    def renumbering(self):
        """
        """
        print("** Renumbering Nodes")
        self._nodes.renumbering(self._elements)
        #for node in single_nodes:
        #    boundary = self._boundaries.node[node]
        #    if not boundary:
        #        self._boundaries.node[ node ] = 'free'
        print("** End Renumbering Nodes")
    #
    # --------------------
    # Mesh Operations
    # --------------------
    #
    def K(self, solver: str|None = None, log: bool = False,
          drop: bool = True):
        """Returns the model's global stiffness matrix.
        
        Solver: numpy/banded/sparse
        """
        #
        jbc= self._nodes.jbc(supports=self._boundaries._nodes)
        neq = self._nodes.neq(supports=self._boundaries._nodes)
        #
        # Banded matrix
        if solver == 'banded':
            iband = self._elements.max_bandwidth(jbc=jbc)
            aa =  assemble_banded_Kmatrix(elements=self._elements , jbc=jbc,
                                          neq=neq, iband=iband)
        else:
            # numpy matrix
            aa = assemble_Kmatrix_np(elements=self._elements, jbc=jbc, neq=neq)            
        #
        if drop:
            with open("stfmx.f2u", "wb") as f:
                pickle.dump(jbc, f)
                pickle.dump(aa, f)
        else:
            return jbc, aa     
        #print('---')
        return None, None
        
