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
from steelpy.ufo.load.main import MeshLoad
from steelpy.ufo.mesh.sqlite.main import MeshSQL
from steelpy.ufo.mesh.process.main import (Ke_matrix, Km_matrix,
                                           Kg_matrix, 
                                           Kt_matrix, Kt_matrix_R)
from steelpy.ufo.plot.main import PlotMesh
#
from steelpy.ufo.utils.main import ModelClassBasic # ufoBasicModel, 
from steelpy.ufo.utils.sets import Groups
from steelpy.ufo.utils.node import node_renumbering
#
from steelpy.sections.main import Section
from steelpy.material.main import Material
#
from steelpy.utils.dataframe.main import DBframework
from steelpy.utils.sqlite.utils import create_connection, create_table
from steelpy.utils.sqlite.main import ClassBasicSQL
#
import numpy as np
#
#
#
class ConceptMesh(ModelClassBasic):
    """ Mesh Model Class"""
    __slots__ = ['_component', 'db_file', '_item', '_labels', '_mesh']
    
    def __init__(self, component:str,
                 sql_file:str|None = None):
        """
        """
        super().__init__(component)
        #
        self._item:dict = {}
        self._labels:list = []
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
                 '_plane', '_component', '_build']
    
    def __init__(self, component:int|str, sql_file:str,
                 build: bool):
        """
        """
        super().__init__(db_file=sql_file)
        self._component = component
        self._build = build
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
                mesh_id = self._push_data(conn, name=name, title=title)
        #
        self._item[name] = MeshItem(name=name,
                                    mesh_id=mesh_id,
                                    sql_file=self.db_file,
                                    build=self._build)
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
    __slots__ = ['_name', 'db_file', '_plane', '_mesh_id',
                 '_df', 'data_type', 
                 '_nodes', '_elements', '_materials', '_sections',
                 '_load', '_boundaries', '_eccentricities', '_groups',
                 '_Kmatrix', '_build', '_solution']

    def __init__(self, name:str|int, mesh_id:int,
                 sql_file:str, build: bool):
        """
        """
        super().__init__(name=name,
                         mesh_id=mesh_id,
                         sql_file=sql_file)
        #
        #self._mesh_id = mesh_id
        self._build = build
        #
        # --------------------------------------------------
        # TODO: should plane be out of mesh?
        self._plane = MeshPlane(plane2D=False)     
        #
        # --------------------------------------------------
        #
        self._materials = Material(mesh_id=self._id,
                                    mesh_type=self.data_type, 
                                    db_file=self.db_file)
    
        self._sections = Section(mesh_id=self._id,
                                 mesh_type=self.data_type,
                                 db_file=self.db_file)
        #
        # --------------------------------------------------
        # groups
        self._groups = Groups()
        #     
        self._load = MeshLoad(mesh_type=self.data_type,
                              mesh_id=self._id,
                              name=name,
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
    #
    def plane(self, plane2D: bool) -> None:
        """ """
        self._plane = MeshPlane(plane2D)
    #
    #
    # --------------------
    # Matrix Operations
    # --------------------
    #
    #
    def jbc(self, plane2D: bool|None = None):
        """ """
        if plane2D is None:
            plane2D=self._plane.plane2D
        #print(self._plane.hdisp)
        jbc = self._nodes.jbc(plane2D=plane2D)
        #if self._plane.plane2D:
        #    jbc.drop(columns=['z', 'rx', 'ry'], inplace=True)
        return jbc
    #
    #
    def _D_partition(self, jbc:DBframework):
        """ """
        
        jbc = jbc.stack()
        # List of the indices for degrees of freedom with unknown displacements.
        D1 = [i for i, item in enumerate(jbc)
              if item != 0]
        # List of the indices for degrees of freedom with known displacements.
        D2 = [i for i, item in enumerate(jbc)
              if item == 0]
        #
        return D1, D2    
    #
    # ------------------------------------------
    #
    def Ke(self): # sparse:bool = True
           #condensed: bool = True):
        """
        sparse : True/False
        
        Returns:
            Elements global stiffness matrix.
        """
        Ka = Ke_matrix(elements=self._elements,
                       nodes=self._nodes,
                       plane2D=self._plane.plane2D)
                       #sparse=sparse)
        return Ka
    #
    def Km(self):
        """
        sparse : True/False
        
        Returns:
            Elements global mass matrix
        """
        Ma = Km_matrix(elements=self._elements,
                       nodes=self._nodes,
                       plane2D=self._plane.plane2D)
                       #sparse=sparse)
        return Ma
    #
    # ------------------------------------------
    #
    def KgXX(self, Dn):
        """
        D : Nodal displacement 
        sparse : True/False
        
        Returns:
            Element Global Geometric Stiffness Matrix
        
        """
        head_disp = ['node_name', *self._plane.hdisp]
        #Dn.set_index('node_name', inplace=True)
        Kg = Kg_matrix(elements=self._elements,
                       nodes=self._nodes,
                       Un=Dn[head_disp].set_index('node_name'),
                       plane2D=self._plane.plane2D)
                       #sparse=sparse)
        return Kg
    #
    def Kg(self, Fb, Rb=None):
        """
        D : Nodal displacement 
        sparse : True/False
        
        Returns:
            Element Global Geometric Stiffness Matrix
        """
        plane2D = self._plane.plane2D
        # ---------------------------------------------
        # Beam section
        beam_temp = []
        Rbeam = {}
        for key, beam in self._elements._beams.items():
            Fb_item = Fb.loc[key]
            #
            #if not Rb:
            #    R = [np.eye(3), np.eye(3)]
            #else:
            #    R = Rb[key]
            beam_mod = beam.deformed(Fb=Fb_item['Fb_local'],
                                     R=Fb_item['R'],
                                     L=Fb_item['L'], 
                                     uv=Fb_item['uv'], 
                                     plane2D=plane2D)
            Rbeam[key] = beam_mod.R
            beam_temp.append(beam_mod)
        #
        Kg = Kg_matrix(elements=beam_temp,
                       nodes=self._nodes,
                       #Un=Fb['Fb_local'],
                       plane2D=self._plane.plane2D)
                       #sparse=sparse)
        return Kg, Rbeam      
    #
    def Kt(self, Dn):
        """
        D : Nodal displacement 
        sparse : True/False
        
        Returns:
            Element Global Tangent Stiffness matrix
        """
        head_disp = ['node_name', *self._plane.hdisp]
        
        Kt = Kt_matrix(elements=self._elements,
                       nodes=self._nodes,
                       Un=Dn[head_disp].set_index('node_name'),
                       plane2D=self._plane.plane2D)
                       #sparse=sparse)
        return Kt
    #
    def Kt_R(self, Fb, Rb=None):
        """
        D : Nodal displacement 
        sparse : True/False
        
        Returns:
            Element Global Tangent Stiffness matrix
        """
        #head_disp = ['node_name', *self._plane.hdisp]
        #Dn = Dn[head_disp].set_index('node_name')
        plane2D = self._plane.plane2D
        # ---------------------------------------------
        # Beam section
        beam_temp = []
        Rbeam = {}
        for key, beam in self._elements._beams.items():
            #nodes = beam.connectivity
            #Un_global = np.concatenate((Dn.loc[nodes[0]],
            #                            Dn.loc[nodes[1]]), axis=None)
            #Un_global = Dn.loc[nodes].stack(future_stack=True)
            #Fb_local = Fb.loc[key].to_numpy()
            Fb_item = Fb.loc[key]
            #
            #if not Rb:
            #    R = [np.eye(3), np.eye(3)]
            #else:
            #    R = Rb[key]
            
            beam_mod = beam.deformed(Fb=Fb_item['Fb_local'],
                                     R=Fb_item['R'],
                                     L=Fb_item['L'], 
                                     uv=Fb_item['uv'], 
                                     plane2D=plane2D)
            Rbeam[key] = beam_mod.R
            beam_temp.append(beam_mod)
            
        
        
        Kt = Kt_matrix_R(elements=beam_temp,
                         nodes=self._nodes,
                         Fb=Fb['Fb_local'],
                         plane2D=self._plane.plane2D)
                         #sparse=sparse)
        #
        return Kt, Rbeam
    #
    # ------------------------------------------
    #
    def Dn(self) -> DBframework.DataFrame|None:
        """
        Returns:
            Nodal global displacement dataframe
        """
        #df = DBframework()
        #col_disp = ['x', 'y', 'z', 'rx', 'ry', 'rz'] # self._plane.hdisp
        #
        #
        #col_grp = ['load_name', 'load_id',
        #           'load_level', 'load_title']
        #basic =  self._load._basic._Dnt()
        #comb =  self._load._combination._Dnt()
        #
        Dn = self._boundaries.df
        Dn = Dn.loc[Dn['type'].isin(['constrained',
                                     'displacement'])]
        #1 / 0
        if not Dn.empty:
            Dn.drop(columns=['title'], inplace=True)
            #
            #bnode[col_disp] = bnode[col_disp].apply(lambda x: 3 if x == 0 else x)
            #bnode[col_disp] = bnode[col_disp].apply(lambda x: 0 if x == 1 else x)
            #bnode.loc[bnode[col_disp] == 0, col_disp] = 3
            #
            # TODO : 1 to zero only on reatain condition
            #bnode[col_disp] = bnode[col_disp].mask(bnode[col_disp] == 0, 3)
            #bnode[col_disp] = bnode[col_disp].mask(bnode[col_disp] == 1, 0)
            #
            # Basic Load
            #if not basic.empty:
            #    grouped = basic.groupby(col_grp) 
            #    new = []
            #    for key, item in grouped:
            #        new_rows = bnode.copy()
            #        new_rows[col_grp] = key
            #        item = df.concat([new_rows, item], ignore_index=True)
            #        new.append(item)
            #    basic = df.concat(new, ignore_index=True)
            ##
            ## Combination
            #if not comb.empty:
            #    grouped = comb.groupby(col_grp)
            #    new = []
            #    for key, item in grouped:
            #        new_rows = bnode.copy()
            #        new_rows[col_grp] = key
            #        item = df.concat([new_rows, item], ignore_index=True)
            #        new.append(item)
            #    comb = df.concat(new, ignore_index=True)
        #
        #Dn = df.concat([self._load._basic._Dnt(),
        #                self._load._combination._Dnt()],
        #               ignore_index=True)
        #Dn = df.concat([basic, comb], ignore_index=True)        
        #if Dn.empty:
        #    return Dn
        #
        #1 / 0
        #
        #jbc = self.jbc()
        #dfjbc = jbc[jbc.any(axis=1)]
        #
        #Di = Dn.loc[Dn['node_name'].isin(dfjbc.index)]
        #Di = Di.loc[Di[col_disp].any(axis=1)]
        #
        #col_grp = ['load_name', 'load_id',
        #           'load_level', 'load_title',
        #           'system', 'mesh_name',
        #           'node_name', 'node_index',
        #           *col_disp]
        #Dn = Dn[col_grp]
        return Dn

    #
    def Fn(self) -> DBframework.DataFrame:
        """
        Returns:
            Nodal global force dataframe
        """
        df = DBframework()
        basic = self._load._basic._Fnt()
        comb = self._load._combination._Fnt()
        Fn = df.concat([basic, comb], ignore_index=True)
        #else:  # basic
        #Fn = self._load._basic.Fn()
        #Fn = df.concat([self._load._basic._Fnt(),
        #                self._load._combination._Fnt()],
        #               ignore_index=True)
        if Fn.empty:
            return Fn
        #
        col_force = ['Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz'] #self._plane.hforce
        #jbc = self.jbc()
        #dfjbc = jbc[jbc.any(axis=1)]
        #
        #Fi = Fn.loc[Fn['node_name'].isin(dfjbc.index)]
        #Fi = Fi.loc[Fi[col_force].any(axis=1)]
        #
        #if Fi.empty:
        #    Fi = Fn
        #
        columns = ['load_name', 'load_id',
                   'load_level', 'load_title',
                   'system', 'mesh_name',
                   'node_name', 'node_index',
                   *col_force]
        #Fi = Fn[columns]
        #return Fi
        return Fn[columns]

    #
    # ------------------------------------------
    #
    def _matrix_partitionX(self, Km):
        """
        Partitions a matrix into sub-matrices based on degree of freedom boundary conditions
        
        Km: The un-partitioned matrix (or vector) to be partitioned.

        Return:
            m1 :
            m2 :
        """
        #jbcc = jbcflat.values
        jbc = self.jbc().stack()
        #
        # List of the indices for degrees of freedom with unknown displacements.
        D1 = [i for i, item in enumerate(jbc)
              if item != 0]
        # List of the indices for degrees of freedom with known displacements.
        D2 = [i for i, item in enumerate(jbc)
              if item == 0]
        #
        # 1D vectors
        if Km.shape[1] == 1:
            # Partition the vector into 2 subvectors
            m1 = Km[D1, :]
            m2 = Km[D2, :]
            return m1, m2
        #
        # 2D matrices
        else:
            # Partition the matrix into 4 submatrices
            m11 = Km[D1, :][:, D1]
            m12 = Km[D1, :][:, D2]
            # m21 = Km[D2, :][:, D1]
            # m22 = Km[D2, :][:, D2]
            return m11, m12  # , m21, m22
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
                              'comment': 'Comment',},
                     inplace=True)
        bload = bload[['Load', 'Beam', 'Type', 'a', 'b',
                       'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz',
                       'qx0', 'qy0', 'qz0', 'qx1', 'qy1', 'qz1',
                       'Comment']]
        #
        nload = self._load._case._basic._nodes.df
        nload.rename(columns={'load_name': 'Load',
                              'node_name': 'Node',
                              'element_name': 'Beam',
                              'load_type': 'Type',
                              'comment': 'Comment',},
                     inplace=True)
        # remove element's node load
        nload = nload[nload['Beam'].isnull()]
        nload = nload[['Load', 'Node', 'Type',
                       'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz',
                       'x', 'y', 'z', 'rx', 'ry', 'rz', 'Comment']]
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
    #
    #
#
# TODO : MeshPlane should be eventually removed
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
    @property
    def colrename_f2u(self) -> dict:
        """"""
        cols = {'Fx':'x', 'Fy':'y', 'Fz':'z',
                'Mx':'rx', 'My':'ry', 'Mz':'rz'}
        if self.plane2D:
            cols = {'Fx':'x', 'Fy':'y', 'Mz':'rz'}
        return cols    
#
#