# 
# Copyright (c) 2009-2023 fem2ufo
#
# Python stdlib imports
from __future__ import annotations
from dataclasses import dataclass
from typing import NamedTuple
#import pickle
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
from steelpy.f2uModel.mesh.process.main import Kmatrix, Gmatrix, Mmatrix
#
#
#
#@dataclass
class MeshPlane(NamedTuple):
    #__slots__ = ['plane2D', '_index']
    
    #def __init__(self, plane2D:bool):
    #    """ """
    plane2D:bool
    #
    # -----------------------------------------------
    #
    @property
    def ndof(self) -> int:
        """ node dgree of freedom"""
        return len(self.dof)
    #
    @property
    def dof(self) -> list[str]:
        """ node degree of freedom"""
        if self.plane2D:
            return ['x', 'y', 'rz']
        return ['x', 'y', 'z', 'rx', 'ry', 'rz']
    #
    @property
    def index(self) -> dict:
        """ """
        return {'x': 0, 'y': 1, 'z': 2, 'rx': 3, 'ry': 4, 'rz': 5} 
    #
    # -----------------------------------------------
    #
    def getidx(self) -> list[int]:
        """get indexes"""
        return [self.index[item] for item in self.dof]
    #
    @property
    def index_off(self) -> list[int]:
        """return index"""
        idx3D = set(self.index.values())
        ndof = len(idx3D)
        #dofactive = set(self._index[item] for item in self.dof)
        dofactive = set(self.getidx())
        indexes = list(idx3D - dofactive)
        indexes.extend([item + ndof for item in indexes])
        return list(reversed(indexes)) # [2, 3, 4, 8, 9, 10]
    #
    # -----------------------------------------------
    #
    @property
    def hforce(self) -> list[str]:
        """ return force header"""
        force = ['Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']
        if self.plane2D:
            idx = self.getidx()
            #force = ['Fx', 'Fy', 'Mz']
            force = [force[item] for item in idx]
        return force
    #
    @property
    def hdisp(self) -> list[str]:
        """ return displacement header"""
        disp = ['x', 'y', 'z', 'rx', 'ry', 'rz']
        if self.plane2D:
            idx = self.getidx()
            disp = [disp[item] for item in idx]
        return disp
    #
    # beam 
    #
    @property
    def bhforce(self) -> list[str]:
        """ beam intergration force"""
        bforce = ['F_Vx', 'F_Vy', 'F_Vz', 'F_Mx', 'F_My', 'F_Mz']
        if self.plane2D:
            idx = self.getidx()
            bforce = [bforce[item] for item in idx]
        #
        return ['node_end', *bforce]
    #
    @property
    def bhdisp(self) -> list[str]:
        """ beam intergration displacement"""
        bdisp = ['F_wx', 'F_wy', 'F_wz', 'F_phix', 'F_thetay', 'F_thetaz']
        if self.plane2D:
            idx = self.getidx()
            bdisp = [bdisp[item] for item in idx]
        #
        return ['node_end', *bdisp]
    #
    #
    @property
    def colrename(self) -> dict:
        """"""
        cols = {'x':'Fx', 'y':'Fy', 'z':'Fz',
                'rx':'Mx', 'ry':'My', 'rz':'Mz'}
        if self.plane2D:
            cols = {key: cols[key] for key in self.dof}
        return cols    
#
#
class Mesh:
    """
    mesh[beam_name] = [number, element1, element2, elementn]
    """
    __slots__ = ['_nodes', '_elements', '_load', 'data_type',
                 '_eccentricities', '_boundaries', '_groups',
                 'db_file', '_df', '_Kmatrix', '_plane']

    def __init__(self, materials, sections,
                 mesh_type:str="sqlite",
                 db_file:str|None = None):
        """
        """
        #self._plane = Plane3D()
        self._plane = MeshPlane(plane2D=False)
        #self._ndof = 6
        #
        self.db_file = db_file
        self.data_type = mesh_type
        self._nodes = Nodes(mesh_type=mesh_type,
                            plane=self._plane,
                            db_file=self.db_file)
        self._boundaries = Boundaries(mesh_type=mesh_type,
                                      db_file=self.db_file)
        #
        self._elements = Elements(nodes=self._nodes,
                                  materials=materials,
                                  sections=sections,
                                  plane=self._plane,
                                  mesh_type=mesh_type,
                                  db_file=self.db_file)
        # groups
        self._groups = Groups()
        #
        self._load = Load(nodes=self._nodes,
                          elements=self._elements,
                          plane=self._plane,
                          mesh_type=mesh_type,
                          db_file=self.db_file)
        # Ops
        self._df = DBframework()
        self._Kmatrix:bool = False
        #self._plane2D:bool = False
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
    def plane(self, plane2D: bool) -> None:
        """ """
        self._plane = MeshPlane(plane2D)
        self._nodes.plane = self._plane
        self._elements.plane = self._plane
        self._load.plane = self._plane
    #
    def jbc(self):
        """ """
        return self._nodes.jbc(supports=self._boundaries._nodes)
    #
    def neq(self):
        """number of equations"""
        return self._nodes.neq(supports=self._boundaries._nodes)
    #
    def K(self, solver: str|None = None): # drop: bool = True, 
        """Returns the model's global stiffness matrix.
        
        Solver: numpy/banded/sparse
        """
        # get data
        #jbc = self._nodes.jbc(supports=self._boundaries._nodes)
        jbc = self.jbc()
        neq = self.neq()
        #
        #if self._plane.m2D:
        #    jbc = jbc[self._plane.dof]
        #
        Ka = Kmatrix(elements=self._elements, jbc=jbc, neq=neq,
                     ndof=self._plane.ndof, solver=solver)
        #
        #if drop:
        #    with open("stfmx.f2u", "wb") as f:
        #        pickle.dump(jbc, f)
        #        pickle.dump(aa, f)
        #else:
        return Ka     
        #print('---')
        #return None, None
    #
    def Kg(self, solver: str|None = None):
        """ Element global geometric stiffness matrix"""
        jbc = self.jbc()
        neq = self.neq()
        Kg = Gmatrix(elements=self._elements, jbc=jbc, neq=neq,
                     plane=self._plane.ndof, solver=solver)
        return Kg
    #
    def M(self, solver: str|None = None):
        """ Element global geometric stiffness matrix"""
        jbc = self.jbc()
        neq = self.neq()
        Ma = Mmatrix(elements=self._elements, jbc=jbc, neq=neq,
                     plane=self._plane.ndof, solver=solver)
        return Ma
    #
    #