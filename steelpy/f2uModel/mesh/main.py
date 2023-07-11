# 
# Copyright (c) 2009-2023 fem2ufo
#
# Python stdlib imports
from __future__ import annotations
import pickle
import re
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
from steelpy.f2uModel.mesh.process.main import Kmatrix
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
    def nodes(self, values:None|list|tuple=None,
              df=None):
        """
        """
        supports = self._boundaries.supports()
        #if values:
        if isinstance(values, (list,tuple)):
            if isinstance(values[0], (list,tuple)):
                for value in values:
                    self._nodes[value[0]] = value[1:4]
                    try:
                        if isinstance(value[4], str):
                            supports[value[0]] = value[4]
                        else:
                            supports[value[0]] = value[4:]
                    except IndexError:
                        pass
            else:
                self._nodes[values[0]] = values[1:4]
                try:
                    if isinstance(values[4], str):
                        supports[values[0]] = values[4]
                    else:
                        supports[values[0]] = values[4:]
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
    def elements(self, values:None|list|tuple=None,
                 df=None):
        """
        """
        if isinstance(values, (list, tuple)):
            if isinstance(values[0], (list,tuple)):
                for value in values:
                    self._elements[value[0]] = value[1:]
            else:
                self._elements[values[0]] = values[1:]
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
    def boundaries(self, values:None|list|tuple=None,
                   df=None):
        """
        """
        if isinstance(values, (list, tuple)):
            if isinstance(values[0], (list,tuple)):
                for item in values:
                    self._boundaries[item[0]] = item[1:]
            else:
                self._boundaries[values[0]] = values[1:]
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
    def load(self, values:None|list|tuple=None,
             df=None):
        """
        """
        if isinstance(values, (list, tuple)):
            if isinstance(values[0], (list,tuple)):
                for item in values:
                    if re.match(r"\b(basic(\_)?(load)?)\b", item[0], re.IGNORECASE):
                        self._load.basic(item[1:])
                    elif re.match(r"\b(comb(ination)?(\_)?(load)?)\b", item[0], re.IGNORECASE):
                        self._load.combination(item[1:])
                    else:
                        raise IOError(f'load {item[0]}')
            else:
                if re.match(r"\b(basic(\_)?(load)?)\b", values[0], re.IGNORECASE):
                    self._load.basic(values[1:])
                elif re.match(r"\b(comb(ination)?(\_)?(load)?)\b", values[0], re.IGNORECASE):
                    self._load.combination(values[1:])
                else:
                    raise IOError(f'load {values[0]}')                
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
          drop: bool = True, m2D:bool = False):
        """Returns the model's global stiffness matrix.
        
        Solver: numpy/banded/sparse
        """
        # get data
        jbc = self._nodes.jbc(supports=self._boundaries._nodes)
        neq = self._nodes.neq(supports=self._boundaries._nodes)
        #
        if m2D:
            jbc = jbc[['x', 'y', 'rz']]
        #
        aa = Kmatrix(elements=self._elements, jbc=jbc, neq=neq,
                     solver=solver, m2D=m2D)
        #
        if drop:
            with open("stfmx.f2u", "wb") as f:
                pickle.dump(jbc, f)
                pickle.dump(aa, f)
        else:
            return jbc, aa     
        #print('---')
        return None, None
        
