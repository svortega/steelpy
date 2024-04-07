# 
# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
#from dataclasses import dataclass
from datetime import datetime as dt
from typing import NamedTuple
#import pickle
#import re
import os
#

# package imports
from steelpy.ufo.process.main import ufoBasicModel, ModelClassBasic
from steelpy.ufo.load.main import MeshLoad
from steelpy.ufo.mesh.sqlite.nodes import NodeSQL
from steelpy.ufo.mesh.sqlite.elements import ElementsSQL
from steelpy.ufo.mesh.sqlite.boundary import BoundarySQL
from steelpy.ufo.mesh.process.main import Ke_matrix, Kg_matrix, Km_matrix, Kt_matrix
from steelpy.ufo.mesh.elements.sets import Groups
from steelpy.ufo.plot.main import PlotMesh
#
from steelpy.sections.main import Section
from steelpy.material.main import Material
#
from steelpy.utils.dataframe.main import DBframework
from steelpy.utils.io_module.inout import check_input_file
from steelpy.utils.sqlite.utils import create_connection, create_table
#

#
#
class ConceptMesh(ModelClassBasic):
    """ Mesh Model Class"""
    __slots__ = ['_name', 'db_file', '_item']
    
    def __init__(self, name:str|None = None,
                 sql_file:str|None = None):
        """
        """
        super().__init__()
        self._name = name
        self._item:dict = {}
    #
    def __setitem__(self, name: int|str, title: int|str) -> None:
        """
        """
        try:
            self._labels.index(name)
            raise Exception(f'Item {name} already exist')
        except ValueError:
            self._labels.append(name)
            self._item[name] = Mesh(name=self._name,
                                    component=name)
            #            
            self._item[name]._set_type(component=name,
                                       comp_type='concept',
                                       title=title)
    #
    def __getitem__(self, name: int|str):
        """ """
        try:
            self._labels.index(name)
            return self._item[name]
        except ValueError:
            raise IndexError(f'Item {name} not valid')
    #
    #
    #
#
#
class Mesh(ufoBasicModel):
    """
    mesh[beam_name] = [number, element1, element2, elementn]
    """
    __slots__ = ['_name', 'db_file', '_plane', '_component', 
                 '_df', 'data_type',
                 '_nodes', '_elements', '_materials', '_sections',
                 '_load', '_boundaries', '_eccentricities', '_groups',
                 '_Kmatrix', '_build', '_solution']

    def __init__(self, name:str|int|None = None,
                 component:str|int|None = None, 
                 sql_file:str|None = None):
        """
        """
        mesh_type:str="sqlite"
        #super().__init__()
        #
        self._build = True
        if sql_file:
            #print('--')
            sql_file = check_input_file(file=sql_file,
                                        file_extension="db")
            self.db_file = sql_file
            self._build = False
            # fixme: name
            #self._name = sql_file.split('.')[0]
        else:
            self.db_file = self._get_file(name=name)
            self.data_type = mesh_type
            self._name = name
            # FIXME : how to handle component? 
            if not component:
                component = name
            #
            conn = create_connection(self.db_file)
            with conn:
                self._create_table(conn)
                comp_no = self._push_data(conn, component)
        #
        #self._component = component
        #
        # --------------------------------------------------
        #
        #self._plane = Plane3D()
        self._plane = MeshPlane(plane2D=False)
        #self._ndof = 6        
        #
        #self._materials = materials
        #self._sections = sections
        #
        # --------------------------------------------------
        #
        self._materials = Material(component=comp_no,
                                    mesh_type=mesh_type, 
                                    db_file=self.db_file)
    
        self._sections = Section(component=comp_no,
                                  mesh_type=mesh_type, 
                                  db_file=self.db_file)       
        #
        # --------------------------------------------------
        #
        self._nodes = NodeSQL(db_system=mesh_type,
                              plane=self._plane,
                              component=comp_no, 
                              db_file=self.db_file)
        #
        self._boundaries = BoundarySQL(db_system=mesh_type,
                                       component=comp_no,
                                       db_file=self.db_file)
        #
        self._elements = ElementsSQL(plane=self._plane,
                                     db_system=mesh_type,
                                     component=comp_no,
                                     db_file=self.db_file)
        #
        # --------------------------------------------------
        # groups
        self._groups = Groups()
        #
        #self._load:Load|None = None
        #nodes=self._nodes,
        #elements=self._elements,
        #boundaries=self._boundaries,         
        self._load = MeshLoad(plane=self._plane,
                              mesh_type=mesh_type,
                              component=comp_no,
                              db_file=self.db_file)
        #
        # --------------------------------------------------
        # Ops
        self._df = DBframework()
        self._Kmatrix:bool = False
        #self._plane2D:bool = False
        # mesh
        #self._plot = PlotMesh(mesh=self)
        #
        # Solution
        #self._solution =  UnSQL(load=self._load,
        #                        db_file=self.db_file,)
        #
    #
    #
    def _get_file(self, name: str):
        """ """
        #BASE_DIR = os.path.dirname(os.path.abspath(__file__))
        filename = name + ".db"
        path = os.path.abspath(filename)
        #self.db_file = path
        #directory = os.path.dirname(path)
        #
        #self.db_file:str = component + "_f2u.db"
        #if mesh_type != "ttt": #"inmemory":
        try: # remove file if exist
            os.remove(path)
        except FileNotFoundError:
            pass
        #
        return path
    #
    def check_file(self, name: str):
        """ """
        pass
    #
    # --------------------
    # Mesh items
    # -------------------- 
    #
    def node(self, values:None|list|tuple=None,
              df=None):
        """
        """
        supports = self._boundaries.support()
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
            self._nodes.df = df
        except AttributeError:
            pass
        #
        return self._nodes

    #
    #
    def boundary(self, values:None|list|tuple=None,
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
            self._boundaries.df = df
        except AttributeError:
            pass 
        #
        return self._boundaries
    #
    #def groups(self):
    #    """
    #    """
    #    return self._groups
    #
    #@property
    #def sections(self):
    #    """ """
    #    1 / 0
    #    return self._sections
    #
    #@property
    #def materials(self):
    #    """ """
    #    1 / 0
    #    return self._materials
    #    
    #
    # ------------------
    # SQL ops
    # ------------------
    #
    def _create_table(self, conn) -> None:
        """ """
        # conn = create_connection(self.db_file)
        table = "CREATE TABLE IF NOT EXISTS Component (\
                    number INTEGER PRIMARY KEY NOT NULL,\
                    name NOT NULL,\
                    type TEXT NOT NULL, \
                    units TEXT NOT NULL,\
                    superelement DECIMAL NOT NULL , \
                    date TEXT NOT NULL,\
                    title TEXT);"
        #
        create_table(conn, table)
        #
    #
    #
    def _push_data(self, conn,
                   component: str|int|None,
                   title: str|None = None):
        """ """
        #
        table = 'INSERT INTO Component(name, type, units,\
                                       superelement, date, title)\
                            VALUES(?,?,?,?,?,?)'
        #
        date = dt.now().strftime('%Y-%m-%d')
        data = (component, 'mesh', 'si', 0, date, title)
        # push
        cur = conn.cursor()
        out = cur.execute(table, data)
        return out.lastrowid
    #
    #
    #
    def _set_type(self, component: str|int,
                  comp_type: str, title: str|None):
        """ """
        
        time=dt.now().strftime('%Y-%m-%d')
        #item = 'concept'
        #
        query = (time, comp_type, title, component)
        table = f"UPDATE Component \
                 SET date = ?, \
                     type = ?, \
                     title = ? \
                 WHERE name = ?;"
        #
        conn = create_connection(self.db_file)
        with conn:          
            cur = conn.cursor()
            comp = cur.execute(table, query)
        #
        if not comp:
            raise IOError(f' component {component} not valid')
    #
    # --------------------
    # Mesh Operations
    # --------------------
    #
    def renumbering(self):
        """
        """
        print("** Renumbering Node")
        self._nodes.renumbering(self._elements)
        #for node in single_nodes:
        #    boundary = self._boundaries.node[node]
        #    if not boundary:
        #        self._boundaries.node[ node ] = 'free'
        print("** End Renumbering Node")
    #
    def build(self):
        """ """
        # sections 
        # self._sections.get_properties()
        # hydro
        #self._load._hydro.solve()
        #
        # FIXME : remove beam nodal load instead --> calculare node load
        if self._build :
            #for key, item in self._mesh.items():
            # FIXME: Wave
            self._load._hydro.process()
            # TODO : remove second _load for simplification
            self._load._basic.FER(elements= self._elements)        
        #
    #
    # --------------------
    # Matrix Operations
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
        #supports=self._boundaries._nodes
        return self._nodes.jbc()
    #
    #def neq(self):
    #    """number of equations"""
    #    return self._nodes.neq(supports=self._boundaries._nodes)
    #
    def Ke(self, sparse: bool = True):
           #condensed: bool = True):
        """Returns the model's global stiffness matrix.
        
        Solver: numpy/banded/sparse
        """
        # get data
        Ka = Ke_matrix(elements=self._elements,
                       nodes=self._nodes,
                       ndof=self._plane.ndof,
                       #condensed=condensed,
                       sparse=sparse)
        #
        return Ka
    #
    def Kg(self, D, sparse: bool = True):
        """
        Element global geometric stiffness matrix
        
        D : 
        """
        #D = D.set_index('node_name', inplace=True)
        Kg = Kg_matrix(elements=self._elements,
                       nodes=self._nodes,
                       D=D,
                       ndof=self._plane.ndof,
                       sparse=sparse)
        return Kg
    #
    def Kt(self, D, sparse: bool = True):
        """ Element Tangent Stiffness matrix"""
        Kt = Kt_matrix(elements=self._elements,
                       nodes=self._nodes,
                       D=D,
                       ndof=self._plane.ndof,
                       sparse=sparse)
        return Kt
    #
    def Km(self, sparse: bool = True):
        """ Element global mass matrix"""
        #
        Ma = Km_matrix(elements=self._elements,
                       nodes=self._nodes,
                       ndof=self._plane.ndof,
                       sparse=sparse)
        return Ma
    #
    # --------------------
    # Plotting
    # --------------------
    #
    def plot(self, figsize:tuple = (10, 10)):
        """ """
        #print('--')
        #self._plot.figsize = figsize
        #return self._plot
        return PlotMesh(cls=self, figsize=figsize)
    #
    #
    # --------------------
    # Tools
    # --------------------
    #
    #
    def to_excel(self, name: str|None = None):
        """dump mesh in excel file"""
        #
        if not name:
            name = self.db_file.split('.')
            name = f"{name[0]}.xlsx"
        #
        #
        nodes = self._nodes.df
        nodes.set_index('name', inplace=True)
        bound = self._boundaries.df
        bound.set_index('name', inplace=True)
        #
        nodes = nodes.join(bound)
        nodes.reset_index(inplace=True)
        #
        # Loading section
        #
        bload = self._load._case._basic._beams.df
        bload.rename(columns={'L0': 'a', 'L1': 'b', 'load_name': 'Load',
                              'element_name': 'Beam', 'load_type': 'Type',
                              'load_comment': 'Title',},
                     inplace=True)
        bload = bload[['Load', 'Beam', 'Type', 'a', 'b',
                       'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz',
                       'qx0', 'qy0', 'qz0', 'qx1', 'qy1', 'qz1',
                       'Title']]
        #
        nload = self._load._case._basic._nodes.df
        nload.rename(columns={'load_name': 'Load',
                              'node_name': 'Node',
                              'element_name': 'Beam',
                              'load_type': 'Type',
                              'load_comment': 'Title',},
                     inplace=True)
        # remove element's node load
        nload = nload[nload['Beam'].isnull()]
        nload = nload[['Load', 'Node', 'Type',
                       'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz',
                       'x', 'y', 'z', 'rx', 'ry', 'rz', 'Title']]
        #
        #
        with self._df.ExcelWriter(name) as writer:
            mat = self._materials.df
            mat.to_excel(writer, sheet_name='Material', index=False)
            #
            sect = self._sections.df
            sect.to_excel(writer, sheet_name='Section', index=False)
            #
            nodes.to_excel(writer, sheet_name='Node', index=False)
            #
            memb = self._elements.df
            memb.to_excel(writer, sheet_name='Element', index=False)
            #
            #nload = self._load._basic._nodes.df
            nload.to_excel(writer, sheet_name='Node_Load', index=False)
            #
            #bline = self._load._basic._beams.line.df
            #bline.to_excel(writer, sheet_name='Beam_Line_Load', index=False)
            #
            #bpoint = self._load._basic._beams.point.df
            #bpoint.to_excel(writer, sheet_name='Beam_Point_Load', index=False)
            #
            bload.to_excel(writer, sheet_name='Beam_Load', index=False)
        #
        #
        print(f'--- end writing excel {name}')
        #1 / 0
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