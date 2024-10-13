# 
# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
#from dataclasses import dataclass
from datetime import datetime as dt
from typing import NamedTuple
import time
#import os
#

# package imports
from steelpy.ufo.process.main import ModelClassBasic # ufoBasicModel, 
from steelpy.ufo.load.main import MeshLoad
from steelpy.ufo.mesh.sqlite.main import MeshSQL
from steelpy.ufo.mesh.process.main import Ke_matrix, Kg_matrix, Km_matrix, Kt_matrix
from steelpy.ufo.plot.main import PlotMesh
#
from steelpy.ufo.process.sets import Groups
from steelpy.ufo.process.node import node_renumbering
#
from steelpy.sections.main import Section
from steelpy.material.main import Material
#
from steelpy.utils.dataframe.main import DBframework
from steelpy.utils.sqlite.utils import create_connection, create_table
from steelpy.utils.sqlite.main import ClassBasicSQL
#
#
#
#
class ConceptMesh(ModelClassBasic):
    """ Mesh Model Class"""
    __slots__ = ['_component', 'db_file', '_item']
    
    def __init__(self, component:str,
                 sql_file:str|None = None):
        """
        """
        super().__init__()
        self._component = component
        self._item:dict = {}
        self._mesh = Mesh(component, sql_file=sql_file)
    #
    def __setitem__(self, name: int|str, title: int|str) -> None:
        """
        """
        try:
            self._labels.index(name)
            raise Exception(f'Item {name} already exist')
        except ValueError:
            self._labels.append(name)
            self._item[name] = MeshItem(name=self._name,
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
#
class Mesh(ClassBasicSQL):
    __slots__ = ['_name', 'db_file', '_item', 
                 '_plane', '_component']
    
    def __init__(self, component:int, sql_file:str):
        """
        """
        super().__init__(db_file=sql_file)
        self._component = component
        self._item:dict = {}
        #print('--')
    #
    #
    # --------------------------------------------
    #
    #
    @property
    def _labels(self):
        """ """
        query = (self._component, )
        table = 'SELECT name FROM Mesh \
                 WHERE component_id = ? \
                 ORDER BY number ASC ;'
        #
        conn = create_connection(self.db_file)
        with conn:        
            cur = conn.cursor()
            cur.execute(table, query)
            items = cur.fetchall()
        return [item[0] for item in items]
    #
    def __setitem__(self, name: int|str, title: str) -> None:
        """
        """
        try:
            self._labels.index(name)
            raise Exception(f' mesh {name} already exist')
        except ValueError:        
            #
            conn = create_connection(self.db_file)
            with conn:
                number = self._push_data(conn, name=name, title=title)
        #
        self._item[name] = MeshItem(name=name, component=number,
                                    sql_file=self.db_file)
        #1 / 0
    #
    def __getitem__(self, name: str|int) -> tuple:
        """
        node_number : node number
        """
        try:
            self._labels.index(name)
            return self._item[name]
        except ValueError:
            raise IndexError(f' mesh : {name} not valid')
    #
    #
    # --------------------------------------------
    # SQL ops
    #
    def _new_table(self, conn) -> None:
        """ """
        table = "CREATE TABLE IF NOT EXISTS Mesh (\
                    number INTEGER PRIMARY KEY NOT NULL,\
                    name NOT NULL,\
                    component_id INTEGER NOT NULL REFERENCES Component(number), \
                    plane TEXT NOT NULL, \
                    units TEXT NOT NULL,\
                    superelement DECIMAL NOT NULL , \
                    date TEXT NOT NULL,\
                    title TEXT);"
        create_table(conn, table)
    #
    def _push_data(self, conn,
                   name: str|int,
                   title: str, 
                   plane: str = '3D'):
        """ """
        table = 'INSERT INTO Mesh(name, component_id, plane, units,\
                                  superelement, date, title)\
                            VALUES(?,?,?,?,?,?,?)'
        #
        date = dt.now().strftime('%Y-%m-%d')
        data = (name, self._component, plane, 'si', 0, date, title)
        # push
        cur = conn.cursor()
        out = cur.execute(table, data)
        return out.lastrowid
    #
    def _pull_data(self, conn, name: int|str, item: str = '*'):
        """ """
        project = (name, self._component)
        table = f'SELECT Mesh.{item}, \
                    Component.name \
                FROM Mesh, Component \
                WHERE Mesh.name = ? \
                AND Mesh.component_id = ?'
        cur = conn.cursor()
        cur.execute(table, project)
        record = cur.fetchone()
        return record
    #
    def _set_type(self, component: str|int,
                  comp_type: str, title: str|None):
        """ """
        
        time=dt.now().strftime('%Y-%m-%d')
        #item = 'concept'
        #
        query = (time, comp_type, title, component)
        table = f"UPDATE Mesh \
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
            raise IOError(f' mesh {component} not valid')
    #
    def build(self):
        """ """
        #print('{:}'.format(52 * '-'))
        for key, item in self._item.items():
            #print(f' building mesh : {key}')
            #print('{:}'.format(52 * '-'))
            item.build()
        #print('{:}'.format(52 * '-'))
#
#
class MeshItem(MeshSQL):
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
        super().__init__(name=name,
                         component=component,
                         sql_file=sql_file)
        #
        #
        # --------------------------------------------------
        # TODO: should plane be out of mesh?
        self._plane = MeshPlane(plane2D=False)     
        #
        # --------------------------------------------------
        #
        self._materials = Material(component=self._component,
                                    mesh_type=self.data_type, 
                                    db_file=self.db_file)
    
        self._sections = Section(component=self._component,
                                  mesh_type=self.data_type, 
                                  db_file=self.db_file)       
        #
        # --------------------------------------------------
        # groups
        self._groups = Groups()
        #     
        self._load = MeshLoad(mesh_type=self.data_type,
                              component=self._component,
                              db_file=self.db_file)
        #
        # --------------------------------------------------
        # Ops
        self._df = DBframework()
        #
    #
    #
    # --------------------
    # Mesh items
    # -------------------- 
    #
    def node(self, values:None|list|tuple|dict = None,
              df=None):
        """
        """
        if values:
            if isinstance(values, dict):
                mname = values['name']
                if isinstance(mname, (list | tuple)):
                    db = DBframework()
                    dfnew = db.DataFrame(values)
                    self._nodes.df = dfnew
                else:
                    mname = values.pop('name')
                    self._nodes[mname] = values
            
            elif isinstance(values, (list,tuple)):
                if isinstance(values[0], (list,tuple)):
                    for value in values:
                        self._nodes[value[0]] = value[1:]
                elif isinstance(values[0], dict):
                    for value in values:
                        self._nodes[value['name']] = value
                else:
                    self._nodes[values[0]] = values[1:]
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
    #
    # --------------------
    # Mesh Operations
    # --------------------
    #
    def renumbering(self):
        """
        """
        print('{:}'.format(52 * '-'))
        start_time = time.time()
        #print("-- Renumbering Node")
        new_number = node_renumbering(nodes=self._nodes,
                                      elements=self._elements)
        self._nodes.renumbering(new_number)
        #for node in single_nodes:
        #    boundary = self._boundaries.node[node]
        #    if not boundary:
        #        self._boundaries.node[ node ] = 'free'
        #print("-- End Renumbering Node")
        uptime = time.time() - start_time
        print(f"** Renumbering Node: {uptime:1.4e} sec")
    #
    def build(self):
        """ """
        print('{:}'.format(52 * '-'))
        print(f' building mesh : {self._name}')
        print('{:}'.format(52 * '-'))
        start_time = time.time()
        # mesh checks
        orphan_nodes = self._nodes.orphan(self._elements)
        if orphan_nodes:
            raise IOError(f' nodes {orphan_nodes} orphan')
        #
        # FIXME : remove beam nodal load instead --> calculate node load
        if self._build :
            self._load._hydro.process()
            # TODO : remove second _load for simplification
            self._load._basic.FER(elements=self._elements)
        #
        uptime = time.time() - start_time
        print(f"** Mesh Building: {uptime:1.4e} sec")
        print('{:}'.format(52 * '-'))
    #
    # --------------------
    # Matrix Operations
    # --------------------
    #
    def plane(self, plane2D: bool) -> None:
        """ """
        self._plane = MeshPlane(plane2D)
        #self._nodes.plane = self._plane
        #self._elements.plane = self._plane
        #self._load.plane = self._plane
    #
    def jbc(self):
        """ """
        #supports=self._boundaries._nodes
        return self._nodes.jbc(dof=self._plane.dof)
    #
    #def neq(self):
    #    """number of equations"""
    #    return self._nodes.neq(supports=self._boundaries._nodes)
    #
    def Ke(self, sparse:bool = True):
           #condensed: bool = True):
        """Returns the model's global stiffness matrix.
        
        Solver: numpy/banded/sparse
        """
        # get data
        Ka = Ke_matrix(elements=self._elements,
                       nodes=self._nodes,
                       plane2D=self._plane.plane2D, 
                       #dof=self._plane.dof,
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
                       plane2D=self._plane.plane2D, 
                       #ndof=self._plane.ndof,
                       sparse=sparse)
        return Kg
    #
    def Kt(self, D, sparse: bool = True):
        """ Element Tangent Stiffness matrix"""
        Kt = Kt_matrix(elements=self._elements,
                       nodes=self._nodes,
                       D=D,
                       plane2D=self._plane.plane2D, 
                       #ndof=self._plane.ndof,
                       sparse=sparse)
        return Kt
    #
    def Km(self, sparse: bool = True):
        """ Element global mass matrix"""
        #
        Ma = Km_matrix(elements=self._elements,
                       nodes=self._nodes,
                       plane2D=self._plane.plane2D, 
                       #ndof=self._plane.ndof,
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