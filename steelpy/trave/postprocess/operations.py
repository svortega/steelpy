# 
# Copyright (c) 2009 fem2ufo
#
# Python stdlib imports
from __future__ import annotations
from collections import namedtuple
from dataclasses import dataclass
import time
from typing import NamedTuple
from datetime import datetime as dt
#
#
# package imports
from steelpy.trave.beam.main import BeamBasic
from steelpy.trave.postprocess.beam import Beam
from steelpy.trave.postprocess.nodes import Nodes
#
from steelpy.utils.math.operations import linstep
from steelpy.utils.dataframe.main import DBframework
from steelpy.utils.sqlite.utils import create_connection, create_table
#
#from steelpy.trave.beam.roark.chapter10.table103_point import BTOpenSupports
#
import numpy as np
#
#
#
#
# -----------------------------------------------------------
# 
# -----------------------------------------------------------
#
#
@dataclass
class Results:
    __slots__ = ['_mesh', '_node', '_beam', 'shell', '_db_file']
    
    def __init__(self, mesh, db_file: str) -> None:
        """
        """
        self._mesh = mesh
        self._db_file = db_file
        self._node = Nodes(mesh=mesh, db_file=self._db_file)
        self._beam = Beam(mesh=mesh, db_file=self._db_file)
    #
    #def __call__(self, noderes, nodereac, beamres):
    #    """ """
    #    plane = self._mesh._plane
    #    self.node:tuple = Node(results=noderes, 
    #                           reac=nodereac,
    #                           plane=plane)
    #    #
    #    self.beam:tuple = Beam(results=beamres, plane=plane)        
    #
    def __str__(self) -> str:
        """ """
        output = "\n"
        output += self._node.__str__()
        output += self._beam.__str__()
        return output
    #
    #
    #def reactions(self):
    #    """ """
    #    pass
    #
    def nodes(self):
        """node results"""
        return self._node
    #
    def beam(self):
        """node results"""
        return self._beam
#
#
# -----------------------------------------------------------
#  Elements
# -----------------------------------------------------------
#
#
#
@dataclass
class MainProcess:
    __slots__ = ['_plane', '_mesh', '_db_file', '_name',
                 '_results']
    # '_load',  '_ndisp'
    
    def __init__(self, mesh, name: str, db_file: str) -> None:
        """
        """
        self._plane = mesh._plane
        self._mesh = mesh
        self._db_file = db_file
        self._name = name
        #
        #
        conn = create_connection(self._db_file)
        with conn:
            self._create_table(conn)        
        #      
        #
        self._results = Results(mesh=mesh, 
                                db_file=self._db_file)

    #
    # -----------------------------------------------------------
    #
    #    
    #
    # -----------------------------------------------------------
    # SQL ops
    #
    def _main_table(self, conn) -> None:
        """ """
        #
        table = "CREATE TABLE IF NOT EXISTS Result (\
                    number INTEGER PRIMARY KEY NOT NULL,\
                    name NOT NULL,\
                    type TEXT,\
                    mesh_name TEXT,\
                    units TEXT NOT NULL,\
                    plane TEXT NOT NULL,\
                    date TEXT);"
        #
        create_table(conn, table)
        #
        #
        table = 'INSERT INTO Component(name, type, mesh_name, \
                                     units, plane, date)\
                                     VALUES(?,?,?,?,?,?)'
        #
        time=dt.now().strftime('%Y-%m-%d')
        plane = '3D'
        if self._mesh._plane.plane2D:
            plane = '2D'
        data = (self._name, None, self._mesh._name, # name, type, mesh_name
                'si', plane, time)                  # units, plane, date
        # push
        cur = conn.cursor()
        cur.execute(table, data)         
    #
    #
    def _create_table(self, conn) -> None:
        """ """
        #       
        #
        # ---------------------------------------------------
        # Node
        # ---------------------------------------------------
        #
        # Node end forces
        table_nodes = "CREATE TABLE IF NOT EXISTS ResultNodeForce (\
                        number INTEGER PRIMARY KEY NOT NULL,\
                        load_name NOT NULL,\
                        component_name NOT NULL, \
                        load_level TEXT NOT NULL,\
                        load_system TEXT NOT NULL,\
                        node_name DECIMAL NOT NULL,\
                        Fx DECIMAL,\
                        Fy DECIMAL,\
                        Fz DECIMAL,\
                        Mx DECIMAL,\
                        My DECIMAL,\
                        Mz DECIMAL);"
        #
        create_table(conn, table_nodes)
        #
        #
        # Node Reactions
        table_nodes = "CREATE TABLE IF NOT EXISTS ResultNodeReaction (\
                        number INTEGER PRIMARY KEY NOT NULL,\
                        load_name NOT NULL,\
                        component_name NOT NULL, \
                        load_level TEXT NOT NULL,\
                        load_system TEXT NOT NULL,\
                        node_name DECIMAL NOT NULL,\
                        Fx DECIMAL,\
                        Fy DECIMAL,\
                        Fz DECIMAL,\
                        Mx DECIMAL,\
                        My DECIMAL,\
                        Mz DECIMAL);"
        #
        create_table(conn, table_nodes)        
        #
        #
        # ---------------------------------------------------
        #
        # Member global end forces
        table_nodes = "CREATE TABLE IF NOT EXISTS ResultMemberEndForce (\
                        number INTEGER PRIMARY KEY NOT NULL,\
                        load_name NOT NULL,\
                        component_name NOT NULL, \
                        load_level TEXT NOT NULL,\
                        load_system TEXT NOT NULL,\
                        element_name NOT NULL,\
                        node_name DECIMAL NOT NULL,\
                        Fx DECIMAL,\
                        Fy DECIMAL,\
                        Fz DECIMAL,\
                        Mx DECIMAL,\
                        My DECIMAL,\
                        Mz DECIMAL);"
        #
        create_table(conn, table_nodes)
        #        
        #
        # ---------------------------------------------------
        # Beam
        # ---------------------------------------------------
        #        
        # Beam local deflection
        table_nodes = "CREATE TABLE IF NOT EXISTS ResultBeamU (\
                        number INTEGER PRIMARY KEY NOT NULL,\
                        load_name NOT NULL,\
                        component_name NOT NULL, \
                        load_level TEXT NOT NULL,\
                        load_system TEXT NOT NULL,\
                        element_name NOT NULL,\
                        node_end DECIMAL NOT NULL,\
                        x DECIMAL,\
                        y DECIMAL,\
                        z DECIMAL,\
                        rx DECIMAL,\
                        ry DECIMAL,\
                        rz DECIMAL);"
        #
        create_table(conn, table_nodes)
        #
        # Beam Local internal forces
        table_nodes = ("CREATE TABLE IF NOT EXISTS ResultBeamForce (\
                        number INTEGER PRIMARY KEY NOT NULL,\
                        load_name NOT NULL,\
                        component_name NOT NULL, \
                        load_level TEXT NOT NULL,\
                        load_system TEXT NOT NULL,\
                        element_name NOT NULL,\
                        node_end DECIMAL NOT NULL,\
                        Fx DECIMAL,\
                        Fy DECIMAL,\
                        Fz DECIMAL,\
                        Mx DECIMAL,\
                        My DECIMAL,\
                        Mz DECIMAL,\
                        psi DECIMAL,\
                        B DECIMAL,\
                        Tw DECIMAL);")
        #
        create_table(conn, table_nodes)
        #
        #
        # Beam local stress
        table_nodes = "CREATE TABLE IF NOT EXISTS ResultBeamStress (\
                        number INTEGER PRIMARY KEY NOT NULL,\
                        load_name NOT NULL,\
                        component_name NOT NULL, \
                        load_level TEXT NOT NULL,\
                        load_system TEXT NOT NULL,\
                        element_name NOT NULL,\
                        node_end DECIMAL NOT NULL,\
                        stress_points DECIMAL NOT NULL,\
                        y DECIMAL,\
                        z DECIMAL,\
                        tau_x DECIMAL,\
                        tau_y DECIMAL,\
                        tau_z DECIMAL,\
                        sigma_x DECIMAL,\
                        sigma_y DECIMAL,\
                        sigma_z DECIMAL);"
        #
        create_table(conn, table_nodes)        
        #
        # ---------------------------------------------------
        # Shell
        # ---------------------------------------------------
        #         
    #
    # -----------------------------------------------------------
    # Operations Node, Element forces and displacement
    # -----------------------------------------------------------
    #
    #
    def solve(self, Un, steps: int,
              Pdelta: bool):
        """ """
        elements = self._mesh._elements
        #
        Un = Un.groupby(['load_level'])
        #
        if Pdelta:
            Un = Un.get_group(('combination', ))
            #print('pdelta')
            cload = self._mesh._load._combination
            load_func = cload.function(steps=steps,
                                       Pdelta=Pdelta)
            #df_comb = cload.to_basic()
            benl_df = cload.ENL()
        else:
            Un = Un.get_group(('basic', ))
            # Solve node and element foce and displacements
            #df_nload = basic_load._nodes.df
            bload = self._mesh._load._basic
            load_func = bload.function(steps=steps,
                                       Pdelta=Pdelta)
            benl_df = bload.ENL()
            #df_nload = self._mesh._load._basic.Fn()       
        #
        #
        df_beamf, df_nforce = self.solve_forces(elements, Un,
                                                load_func, benl_df,
                                                steps=steps, Pdelta=Pdelta)
        #
        return df_beamf, df_nforce
    #
    #
    def solve_forces(self, elements, df_ndisp,
                     load_func, benl_df,
                     steps: int, Pdelta: bool):
        """
        Beam's internal forces along its lenght
        """
        start_time = time.time()
        #
        ndgrp = df_ndisp.groupby(['load_name', 'component_name',
                                  'load_level',  'load_system'])
        #
        benlgrp = benl_df.groupby(['load_name', 'component_name',
                                  'load_level', 'load_system'])
        #
        blgrp = load_func.groupby(['load_name', 'component_name',
                                    'load_level'])
        #
        # Dummy Bending [V, M, theta, w]
        Fblank = [0, 0, 0, 0]
        # Dummy Torsion [T, Phi, Psi, B, Tw] - [T, theta, theta1, theta2, theta3]
        Fblank_t = [0, 0, 0, 0, 0]
        dummyf = np.array([0]*self._plane.ndof)
        #
        ndof = self._plane.ndof
        hdisp = ['node_name', *self._plane.hdisp]
        hforce = [*self._plane.hforce]
        #rforce = ['x', 'y', 'z', 'rx', 'ry', 'rz']
        tforce = ['Psi', 'B', 'Tw'] # torsion part
        #hload = ['node_name', *self._plane.hload]
        hload = ['axial', 'torsion', 'VM_inplane', 'VM_outplane']
        #
        dftemp: list = []
        mif_data: list = []
        for key, noded in ndgrp:
            # TODO: select basic only
            #if key[2] != 'basic':
            #    continue
            ndisp = noded[hdisp]
            ndisp.set_index('node_name', inplace=True)
            #
            for mname, element in elements.items():
                nodes = element.connectivity
                Tb = element.T
                #
                material = element.material
                section = element.section.properties(poisson=material.poisson)
                beam = BeamBasic(L=element.L, area=section.area, 
                                 Iy=section.Iy, Iz=section.Iz,
                                 J=section.J, Cw=section.Cw, 
                                 E=material.E, G=material.G,
                                 Asy=section.Asy,
                                 Asz=section.Asz,
                                 Pdelta=Pdelta)
                #
                #nodes = element.connectivity
                # ---------------------------------------------
                # displacement
                nd_global = np.concatenate((ndisp.loc[nodes[0]],
                                            ndisp.loc[nodes[1]]), axis=None)
                #
                # ---------------------------------------------
                # convert global end-node disp in beam's local system
                # [x,y,z,rx,ry,rz]
                nd_local = Tb @ nd_global
                #
                # --------------------------------------------
                # set beam to general response expresions --> R0
                # [V, M, theta, w]
                #
                #
                # convert beam end-node disp to force [F = Kd] in global system
                #FUan = element.Ke_local @ nd_local
                kt_local = element.Kt_local(nd_global)
                FUtn = kt_local @ nd_local
                #
                #TODO: confirm change reactions sign
                eq = NodeGenRespEq(nd_local, self._plane)
                #
                # ---------------------------------------------
                # Beam load (udl/point)
                try:
                    #
                    # [load_name, member_load_title, load_type, load_system, 
                    # beam_number, x, Fx, Fy, Fz]
                    #
                    mbload = blgrp.get_group(key[:3])
                    mbload = mbload.groupby(['element_name'])                    
                    mnload = mbload.get_group((mname, ))
                    mnload = mnload.groupby(['node_end'],
                                            as_index=False)[hload].sum()
                    #
                    benl = benlgrp.get_group(key)
                    benl = benl.groupby(['element_name'])                    
                    bnload = benl.get_group((mname, ))
                    #
                    # Displacement
                    #rload = bnload.groupby(['node_name'],
                    #                        as_index=False)[rforce].sum()
                    #rload.set_index('node_name', inplace=True)                     
                    #
                    # Torsion
                    mtload = bnload.groupby(['node_name'],
                                            as_index=False)[tforce].sum()
                    mtload.set_index('node_name', inplace=True)                    
                    #
                    # get beam end loads
                    #
                    nf_local = bnload.groupby(['node_name'],
                                               as_index=False)[hforce].sum()
                    nf_local.set_index('node_name', inplace=True)                 
                    #
                    FUen = np.concatenate((nf_local.loc[nodes[0]],
                                           nf_local.loc[nodes[1]]), axis=None)
                    #
                    FUen = Tb @ FUen
                    #
                    Pu = FUen - FUtn  
                    #
                    # get reactions with correct sign for postprocessing
                    #
                    #bload0 = nf_local.loc[nodes[0]]
                    #tload0 = mtload.loc[nodes[0]]
                    R0 = eq.R0(bload=Pu, #nf_local.loc[nodes[0]],
                               tload=mtload.loc[nodes[0]])
                    #
                    lbforce = [['local', mname, bstep.node_end,
                                *beam.response(x=bstep.node_end,
                                               R0=[R0.x, R0.t, R0.y, R0.z],
                                               Fx=[*bstep[2:]])]
                               for bstep in mnload.itertuples()]
                    #print('-->')
                except KeyError : # (KeyError, AttributeError, TypeError): # No load on beam
                    # --------------------------------------------
                    # [Fx,Fy,Fz,Mx,My,Mz]
                    Pu = -1 * FUtn
                    R0 = eq.R0(bload=Pu,
                               tload=Tload(0, 0, 0))
                    #
                    # [beam_number, load_title, x, Fx, Fy, Fz, Mx, My, Mz]
                    #
                    Lsteps = linstep(d=element.section.geometry.d,
                                     L=element.L, steps=steps)                    
                    #Lsteps = linspace(start=0, stop=element.L, num=steps+1, endpoint=True)
                    #
                    lbforce = [['local', mname,  xstep,
                                *beam.response(x=xstep, R0=[R0.x, R0.t, R0.y, R0.z],
                                               Fx=[Fblank, Fblank_t, Fblank, Fblank])]
                               for xstep in Lsteps]
                #
                #
                # ---------------------------------------------
                # convert beam end-node disp to force [F = Kd] in global system
                #
                gnforce = element.Ke_global @ nd_global
                #
                # ---------------------------------------------
                #
                dftemp.append([*key, mname, nodes[0], *gnforce[:ndof]])
                dftemp.append([*key, mname, nodes[1], *gnforce[ndof:]])                
                #
                # ---------------------------------------------
                # Member total force in local system
                #
                # Axial   [FP, blank, blank, Fu]
                # Torsion [T, Phi, Psi, B, Tw] - [T, theta, theta1, theta2, theta3]
                # Bending [V, M, theta, w]
                #
                mif_data.extend([[*key[:3], # load_name, component_name, load_level,
                                  *lbf[:3], # load_system, element_name, node_end
                                  *lbf[3],  # Axial
                                  *lbf[4],  # Torsion
                                  *lbf[5],  # Bending in plane
                                  *lbf[6]]  # Bending out plane
                                 for lbf in lbforce])
                #
                # ---------------------------------------------
                #
        #
        #
        # ---------------------------------------------
        #        
        df_nforce = self.dfmend(dftemp=dftemp)
        # Member internal forces
        df_mif = self.df_mif(mif_data)
        #
        #conn = create_connection(self._db_file)
        #with conn:        
        #    df_membf.to_sql('Qf', conn,
        #                 index_label=['load_name', 'load_level', 'system',
        #                              'element_name', 'node_end',
        #                              'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz',
        #                              'x', 'y', 'z', 'rx', 'ry', 'rz'], 
        #                 if_exists='append', index=False)
        #
        uptime = time.time() - start_time
        print(f"** Beam Force Calculation Time: {uptime:1.4e} sec")
        #print('ok')
        return df_mif, df_nforce
    #    
    #
    def df_mif(self, member_load):
        """ member internal forces"""
        # --------------------------------------------
        # Seting member df
        # --------------------------------------------
        #
        #1/0
        header = ['load_name', 'component_name', 'load_level',
                  #'load_id','load_title',
                  'load_system',
                  'element_name', 'node_end',
                  'Fx', 'blank1', 'blank2', 'x',  # axial
                  'Mx', 'rx', 'Psi',  'B',  'Tw', # torsion
                  'Fy', 'Mz', 'rz', 'y',          # bending in plane
                  'Fz', 'My', 'ry', 'z']          # bending out plane
        #
        db = DBframework()
        df_membf = db.DataFrame(data=member_load,
                                columns=header, index=None)
        # reorder columns
        df_membf.drop(columns=['blank1', 'blank2'])
        df_membf = df_membf[['load_name', 'component_name',
                             #'load_id',
                             'load_level',
                             #'load_title', 
                             'load_system',
                             'element_name', 'node_end',
                             'Fx', 'Fy', 'Fz',
                             'Mx', 'My', 'Mz',
                             'x', 'y', 'z',
                             'rx', 'ry', 'rz',
                             'Psi', 'B', 'Tw']]
        #
        #df_membf.rename({'F_Vx': 'Fx','F_Vy': 'Fy', 'F_Vz': 'Fz',
        #                 'F_Mx': 'Mx', 'F_My': 'My', 'F_Mz': 'Mz',
        #                 'F_wx': 'x', 'F_wy': 'y', 'F_wz': 'z',
        #                 'F_thetax': 'rx', 'F_thetay': 'ry', 'F_thetaz': 'rz'},
        #                axis=1, inplace=True)
        #
        return df_membf
    #
    #
    def dfmend(self, dftemp):
        """ """
        hforce = ['node_name', *self._plane.hforce]
        # get df 
        db = DBframework()
        header: list[str] = ['load_name', 'component_name',
                             #'load_id',
                             'load_level', 
                             # 'load_title',
                             'load_system',
                             'element_name',
                             *hforce]
        df_nforce = db.DataFrame(data=dftemp, columns=header, index=None)
        #return df_nforce[['load_name', 'component_name',
        #                  #'load_id','load_title',
        #                  'load_level',
        #                  'load_system', 'element_name',
        #                  *hforce]]
        return df_nforce
    #
    #
    def get_reactions(self, nforce):
        """ get nodal reactions """
        #
        nforce = (nforce.groupby(['load_name', 'component_name',
                                  #'load_id',
                                  'load_level',
                                  #'load_title',
                                  'load_system',
                                  'node_name'],
                                  as_index=False)
                   [self._plane.hforce].sum())
        #
        supports = self._mesh._boundaries.support()
        nname = list(supports)
        #
        nrsupp = nforce.loc[nforce['node_name'].isin(nname)]
        #
        #nrsupp = nforce.loc[nforce['node_name'].isin(nname)]
        #nrsupp2 =  nreacs.loc[nreacs['node_name'].isin(nname),
        #                      ['Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']].sum()
        # 
        #nreaction = (nrsupp.groupby(['load_name', 'load_id', 'load_type',
        #                              'load_system', 'load_title', 'node_name'])
        #             [self._plane.hforce].sum())
        # reset groupby to create new df
        #nreaction = nreaction.reset_index(names=['load_name', 'load_id', 'load_type',
        #                                         'load_system', 'load_title', 'node_name'])
        #
        return nforce, nrsupp    
    #
    #
    # -----------------------------------------------------------
    # Operations Element stress
    #
    #
    def solve_stress(self, bforce_df):
        """get elements stress"""
        #print("** Beam stress calculation")
        start_time = time.time()
        elements = self._mesh._elements
        members = bforce_df.groupby(['element_name',
                                     'load_name', 'component_name',
                                     'load_level', 'load_system'])
        #
        for key, item in members:
            member = elements[key[0].tolist()]
            section = member.section
            material = member.material
            stressdf = section.stress(df=item,
                                      E=material.E,
                                      G=material.G,
                                      poisson=material.poisson)
            conn = create_connection(self._db_file)
            with conn:            
                stressdf.to_sql('ResultBeamStress', conn,
                                if_exists='append', index=False)
        #
        uptime = time.time() - start_time
        print(f"** Beam stress calculation Time: {uptime:1.4e} sec")
    #
    # -----------------------------------------------------------
    # 
    #
    def results(self, Un, Pdelta: bool,
                beam_steps:int = 10):
        """
        beam_steps : Integration points beam element (10 default)
        """
        print("** Postprocessing")
        #      
        #
        df_beamf, df_Qn = self.solve(Un, beam_steps,
                                     Pdelta=Pdelta)
        #
        #
        if not Pdelta:
            # combination
            load_comb = self._mesh._load.combination()
            df_comb = load_comb.to_basic() 
            #
            df_beamf = update_memberdf(dfmemb=df_beamf,
                                       dfcomb=df_comb,
                                       item='node_end', 
                                       values=['Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz',
                                               'x', 'y', 'z', 'rx', 'ry', 'rz',
                                               'Psi', 'B', 'Tw'])
            #
            # End node force combination
            df_Qn = update_memberdf(dfmemb=df_Qn,
                                    dfcomb=df_comb,
                                    item='node_name',
                                    values=self._plane.hforce)
        #
        #
        #
        df_Qbeam = df_beamf[['load_name', 'component_name', 
                             'load_level', 'load_system',
                             'element_name', 'node_end',
                             'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz',
                             #'rx', 'ry', 'rz',
                             'Psi', 'B', 'Tw']]
        #
        df_Dbeam = df_beamf[['load_name', 'component_name', 
                             'load_level', 'load_system',
                             'element_name', 'node_end',
                             'x', 'y', 'z', 'rx', 'ry', 'rz']]
        # Push to sql
        conn = create_connection(self._db_file)
        with conn:        
            df_Qbeam.to_sql('ResultBeamForce', conn,
                            if_exists='append', index=False)
            #
            df_Dbeam.to_sql('ResultBeamU', conn, 
                            if_exists='append', index=False)
        #        
        #
        #
        #header = ['load_name', 'load_level', 'system',
        #          'element_name', 'node_name']
        #header.extend(self._plane.hforce)
        #conn = create_connection(self._db_file)
        #with conn:
        #    df_Qn.to_sql('Qm', conn,
        #                   index_label=header, 
        #                   if_exists='append', index=False)
        #
        #
        # node reacctions
        df_nforce, df_nreac = self.get_reactions(df_Qn)
        #
        # Push to sql
        #
        header = ['load_name', 'component_name', 
                  'load_level', 'load_system',
                  'element_name', 'node_name']
        header.extend(self._plane.hforce)        
        #
        conn = create_connection(self._db_file)
        with conn:
            df_Qn.to_sql('ResultMemberEndForce', conn,
                         index_label=header, 
                         if_exists='append', index=False)
            #
            header.remove('element_name')
            df_nforce.to_sql('ResultNodeForce', conn,
                             index_label=header, 
                             if_exists='append', index=False)
            #
            df_nreac.to_sql('ResultNodeReaction', conn,
                            index_label=header, 
                            if_exists='append', index=False)
        #
        #
        # displacments
        #df_ndisp = update_ndf(dfnode=Un, dfcomb=df_comb,
        #                      values=self._plane.hdisp)
        #
        #header = ['load_name',
        #          #'load_id',
        #          'load_level',
        #          'load_system',
        #          #'load_title',
        #          'node_name']
        #
        #
        #df_nres = Un.set_index(header).join(df_nforce.set_index(header))
        #df_nres = df_nres.reset_index()
        #
        #1 / 0
        #
        #self._results(noderes=df_nres, 
        #              nodereac=df_nreac,
        #              beamres=df_beam)
        #
        #
        #
        #
        # solve element stress
        #
        self.solve_stress(bforce_df=df_Qbeam)        
        #
        #1 / 0
        return self._results
    #    
    #
    # -----------------------------------------------------------
    #
    #
    #
    #
#
#
# -----------------------------------------------------------
# Nodes
# -----------------------------------------------------------
#
#
@dataclass
class NodeGenRespEq:
    """ Local system"""
    __slots__ = ['nd_local', 'plane']
    
    def __init__(self, nd_local: list,
                 plane: bool) -> None:
        """ """
        self.nd_local = nd_local
        self.plane = plane
    #
    def R0(self, bload: list, tload: tuple|list) -> tuple:
        """
        Axial   [FP, blank, blank, Fu]
        Bending [V, M, theta, w]
        Torsion [T, Phi, Psi, B, Tw]
        """
        # select node disp end 0
        nd0 = self.nd_local[:self.plane.ndof]
        bload = bload[:self.plane.ndof]
        #
        # Axial
        axial = bload[0]
        if np.isclose(a=axial, b=0.0, atol=0.01):
            axial = 0
        #
        R0x = [1 * axial,            # Fx
               0, 0,                 # blank, blank
               -1 * nd0[0]]          # dx
        #
        if self.plane.plane2D:
            R0t, R0y, R0z = self.D2(nd0, bload)
        else:
            R0t, R0y, R0z = self.D3(nd0, bload, tload)
        #
        # Inverting sign to calculate beam load along length
        #
        # In plane
        R0y[0] *= -1   # Vy
        #R0y[1] *= -1  # Mz
        #R0y[2] *= -1  # rz
        R0y[3] *= -1   # dy
        #
        # Out plane
        R0z[0] *= -1    # Vz
        # R0z[1] *= -1  # My
        # R0z[2] *= -1  # ry
        R0z[3] *= -1    # dz
        #
        return Req(R0x, R0t, R0y, R0z)
    #
    def D3(self, nd0:list, bload:list,
           tload:tuple):
        """
        3D plane
        nd0 = [dx,dy,dz,rx,ry,rz]
        bload : [Fx,Fy,Fz,mx,my,mz]
        """
        # In plane
        R0y = [bload[1],   # Vy
               bload[5],   # Mz
               nd0[5],     # rz
               nd0[1]]     # dy
        # Out plane
        R0z = [bload[2],  # Vz
               bload[4],  # My
               nd0[4],    # ry
               nd0[2]]    # dz
        #
        # Torsion
        T0 = 1 * bload[3]
        rx = 1 * nd0[3]
        #
        # Thin walled sections (Ibeam, Channel & Z)
        thetas = [1 * tload.Psi,
                  1 * tload.B,
                  1 * tload.Tw]
        #
        R0t = [T0, rx, *thetas]
        #
        return R0t, R0y, R0z
    #
    def D2(self, nd0:list, bload:list):
        """
        2D plane
        nd0 : [dx, dy, rz]
        bload : [Fx,Fy,mz]
        """
        # In plane
        R0y = [bload[1],   # Vy
               bload[2],   # Mz
               nd0[2],     # rz
               nd0[1]]     # dy
        # Out plane
        R0z = [0, 0, 0, 0]  # [Vz, My, ry, dz]
        R0t = [0, 0, 0, 0, 0]  # [T, Phi, Psi, B, Tw]
        return R0t, R0y, R0z
#
#
#
Req = namedtuple('Reactions', ['x', 't', 'y', 'z'])
#Pload = namedtuple('PointLoad', ['Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz'])
Tload = namedtuple('TorsionLoad', ['Psi', 'B', 'Tw'])
#
# -----------------------------------------------------------
# Combination process
# -----------------------------------------------------------
#
#
def update_memberdf(dfmemb, dfcomb,
                    item: str, 
                     values:list[str]):
    """
    Update node displacements to include lcomb
    """
    db = DBframework()
    # group basic load by name
    grp = dfmemb.groupby('load_name')
    combgrp = dfcomb.groupby('load_name')
    for key, combfactors in combgrp:
        for row in combfactors.itertuples():
            comb = grp.get_group(row.basic_load).copy()
            comb.loc[:, values] *= row.factor
            comb['load_level'] = 'combination'
            comb['load_name'] = row.load_name
            comb['load_id'] = row.load_id
            comb['load_title'] = row.load_title
            #
            try:
                dftemp = db.concat([dftemp, comb], ignore_index=True)
            except UnboundLocalError:
                dftemp = comb
    #
    try:
        dftemp = dftemp.groupby(['load_name', 'component_name',
                                 'load_level', 'load_system',
                                 #'load_title',  'load_id',
                                 'element_name', item],
                                  as_index=False)[values].sum()
        #test
        dfmemb = db.concat([dfmemb, dftemp], ignore_index=True)
    except UnboundLocalError:
        pass
    #
    return dfmemb
#
#
#
#

