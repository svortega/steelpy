# 
# Copyright (c) 2009 fem2ufo
#
from __future__ import annotations
# Python stdlib imports
from dataclasses import dataclass
import time
from typing import NamedTuple
#
#
# package imports
from steelpy.trave.beam.main import BeamBasic
from steelpy.trave.postprocessor.beam import Beam
from steelpy.trave.postprocessor.nodes import Nodes
#
from steelpy.utils.math.operations import linspace
from steelpy.utils.dataframe.main import DBframework
from steelpy.utils.sqlite.utils import create_connection, create_table
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
    __slots__ = ['_plane', '_mesh', '_results', '_db_file']
    # '_load',  '_ndisp'
    
    def __init__(self, mesh, db_file: str) -> None:
        """
        """
        self._plane = mesh._plane
        self._mesh = mesh
        self._db_file = db_file
        #
        self._results = Results(mesh=mesh, 
                                db_file=self._db_file)
        #
        conn = create_connection(self._db_file)
        with conn:
            self._create_table(conn)        
    #
    # -----------------------------------------------------------
    #
    #
    #def input(self, load, ndisp):
    #    """ """
    #    self._load = load
    #    self._ndisp = ndisp
    #
    #
    # -----------------------------------------------------------
    # SQL ops
    #
    #
    def _create_table(self, conn) -> None:
        """ """
        #
        # ---------------------------------------------------
        # Node
        # ---------------------------------------------------
        #
        # Node end forces
        table_nodes = "CREATE TABLE IF NOT EXISTS tb_Fn (\
                        number INTEGER PRIMARY KEY NOT NULL,\
                        load_name NOT NULL,\
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
        table_nodes = "CREATE TABLE IF NOT EXISTS tb_Rn (\
                        number INTEGER PRIMARY KEY NOT NULL,\
                        load_name NOT NULL,\
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
        table_nodes = "CREATE TABLE IF NOT EXISTS tb_Fm (\
                        number INTEGER PRIMARY KEY NOT NULL,\
                        load_name NOT NULL,\
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
        table_nodes = "CREATE TABLE IF NOT EXISTS tb_Dbeam (\
                        number INTEGER PRIMARY KEY NOT NULL,\
                        load_name NOT NULL,\
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
        table_nodes = "CREATE TABLE IF NOT EXISTS tb_Qbeam (\
                        number INTEGER PRIMARY KEY NOT NULL,\
                        load_name NOT NULL,\
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
                        B DECIMAL,\
                        Tw DECIMAL);"
        #
        create_table(conn, table_nodes)
        #
        #
        # Beam local stress
        table_nodes = "CREATE TABLE IF NOT EXISTS tb_Sbeam (\
                        number INTEGER PRIMARY KEY NOT NULL,\
                        load_name NOT NULL,\
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
    def solve(self, Un, steps: int):
        """ """
        elements = self._mesh._elements
        #
        #loadCase = self._mesh._load.case()
        basic_load = self._mesh._load._basic
        #df_nload = basic_load.node_df()
        df_nload = basic_load._nodes.df
        bload_func = basic_load.process(elements=elements,
                                        steps=steps) 
        #
        # Solve node and element foce and displacements
        # 
        df_beamf, df_nforce = self.solve_forces(elements, bload_func,
                                                Un, df_nload, steps)
        #
        return df_beamf, df_nforce
    #
    #
    def solve_forces(self, elements, bload_func,
                     df_ndisp, df_nload, steps):
        """
        Beam's internal forces along its lenght
        """
        #
        ndgrp = df_ndisp.groupby(['load_name',  'load_level',
                                  #'load_number', 'load_title',
                                  'load_system'])
        #nfgrp = df_nforce.groupby(['load_name',  'load_type', 'load_number', 'load_title'])
        nlgrp = df_nload.groupby(['load_name', 'load_level',
                                  #'load_number', 'load_title',
                                  'load_system'])
        blgrp = bload_func.groupby(['load_name',  'load_level'])
        #
        # Dummy Bending [V, M, theta, w]
        Fblank = [0, 0, 0, 0]
        # Dummy Torsion [T, B, Psi, Phi, Tw]
        Fblank_t = [0, 0, 0, 0, 0]
        dummyf = np.array([0]*self._plane.ndof)
        #
        ndof = self._plane.ndof
        hdisp = ['node_name', *self._plane.hdisp]
        hforce = ['node_name', *self._plane.hforce]
        #hload = ['node_name', *self._plane.hload]
        hload = ['axial', 'torsion', 'VM_inplane', 'VM_outplane']
        #
        dftemp: list = []
        member_load: list = []
        for key, noded in ndgrp:
            # TODO: select basic only
            if key[1] != 'basic':
                continue
            ndisp = noded[hdisp]
            ndisp.set_index('node_name', inplace=True)        
            #
            # check if basic load
            try:
                mbload = blgrp.get_group(key[:2])
                mbload = mbload.groupby(['element_name'])
                #
                nlval = nlgrp.get_group(key)
                nlval = nlval.groupby(['element_name'])
            except KeyError:
                mbload = {}
            #
            for mname, element in elements.items():
                nodes = element.connectivity
                Tlg = element.T
                #
                material = element.material
                section = element.section.properties()
                beam = BeamBasic(L=element.L, area=section.area, 
                                 Iy=section.Iy, Iz=section.Iz,
                                 J=section.J, Cw=section.Cw, 
                                 E=material.E, G=material.G)         
                #
                #nodes = element.connectivity
                # ---------------------------------------------
                # displacement
                gndisp = np.concatenate((ndisp.loc[nodes[0]],
                                         ndisp.loc[nodes[1]]), axis=None)
                #
                # ---------------------------------------------
                # convert beam end-node disp to force [F = Kd] in global system
                #
                gnforce = element.K @ gndisp
                #
                # ---------------------------------------------
                #
                dftemp.append([*key, mname, nodes[0], *gnforce[:ndof]])
                dftemp.append([*key, mname, nodes[1], *gnforce[ndof:]])                
                #             
                # ---------------------------------------------
                # convert global end-node disp in beam's local system
                #
                lndisp = Tlg @ gndisp
                #
                # --------------------------------------------
                # convert beam end-node disp to force [F = Kd] in global system
                #
                lnforce = element.k @ lndisp
                #
                # --------------------------------------------
                # set beam to general response expresions --> R0
                # [V, M, theta, w]
                #
                #TODO: confirm change reactions sign
                eq = NodeGenRespEq(lnforce, lndisp,
                                   ndof, self._plane.plane2D)
                #
                # ---------------------------------------------
                #
                try:
                    # Beam load (udl/point)
                    # [load_name, member_load_title, load_type, load_system, 
                    # beam_number, x, Fx, Fy, Fz]
                    mnload = mbload.get_group(mname)
                    mnload = mnload.groupby(['node_end'], as_index=False)[hload].sum()
                    #
                    nodeloads = nlval.get_group(mname)
                    nodeloads = nodeloads.groupby(['node_name'], as_index=False)[hforce].sum()
                    nodeloads = nodeloads.set_index('node_name')                 
                    #
                    blitem = np.concatenate((nodeloads.loc[nodes[0]],
                                             nodeloads.loc[nodes[1]]), axis=None)                
                    #
                    blitem = Tlg @ blitem
                    #
                    R0 = eq.R0(bload=blitem)
                    #R0 = eq.R0()
                    lbforce = [['local', mname, bstep.node_end,
                                *beam.response(x=bstep.node_end,
                                               R0=[R0.x, R0.t, R0.y, R0.z],
                                               Fx=[*bstep[2:]])]
                               for bstep in mnload.itertuples()]
                    #print('-->')
                except (KeyError, AttributeError, TypeError):
                    # No load on beam
                    R0 = eq.R0()
                    # [beam_number, load_title, x, Fx, Fy, Fz, Mx, My, Mz]
                    Lsteps = linspace(start=0, stop=element.L, num=steps+1, endpoint=True)
                    lbforce = [['local', mname,  xstep,
                                *beam.response(x=xstep, R0=[R0.x, R0.t, R0.y, R0.z],
                                               Fx=[Fblank, Fblank_t, Fblank, Fblank])]
                               for xstep in Lsteps]
                #
                # ---------------------------------------------
                # Member total force in local system
                #
                # Axial   [FP, blank, blank, Fu]
                # Torsion [T, B, Psi, Phi, Tw]
                # Bending [V, M, theta, w]
                #
                member_load.extend([[*key[:2], # load_name, load_type,
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
        #conn = create_connection(self._db_file)
        #with conn:
        #    self.push_Qm(conn, data=dftemp)
        #    self.push_Qf(conn, data=member_load)
        #
        # ---------------------------------------------
        #        
        df_nforce = self.dfmend(dftemp=dftemp)
        df_membf = self.df_mint(member_load)
        #
        #conn = create_connection(self._db_file)
        #with conn:        
        #    df_membf.to_sql('tb_Qf', conn,
        #                 index_label=['load_name', 'load_level', 'system',
        #                              'element_name', 'node_end',
        #                              'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz',
        #                              'x', 'y', 'z', 'rx', 'ry', 'rz'], 
        #                 if_exists='append', index=False)
        #
        #print('ok')
        return df_membf, df_nforce
    #    
    #
    def df_mint(self, member_load):
        """ """
        # --------------------------------------------
        # Seting member df
        # --------------------------------------------
        #
        header = ['load_name', 'load_level',
                  #'load_number','load_title',
                  'load_system',
                  'element_name', 'node_end',
                  'F_Vx', 'blank1', 'blank2', 'F_wx',          # axial
                  'F_Mx', 'F_B', 'F_psi', 'F_phix', 'F_Tw',    # torsion
                  'F_Vy', 'F_Mz', 'F_thetaz', 'F_wy',          # bending in plane
                  'F_Vz', 'F_My', 'F_thetay', 'F_wz']          # bending out plane
        #
        db = DBframework()
        df_membf = db.DataFrame(data=member_load, columns=header, index=None)
        # reorder columns
        df_membf = df_membf[['load_name',
                             #'load_number',
                             'load_level',
                             #'load_title', 
                             'load_system',
                             'element_name', 'node_end',
                             'F_Vx', 'F_Vy', 'F_Vz',
                             'F_Mx', 'F_My', 'F_Mz',
                             'F_wx', 'F_wy', 'F_wz',
                             'F_phix', 'F_thetay', 'F_thetaz',
                             'F_B', 'F_Tw']]
        #
        df_membf.rename({'F_Vx': 'Fx','F_Vy': 'Fy', 'F_Vz': 'Fz',
                         'F_Mx': 'Mx', 'F_My': 'My', 'F_Mz': 'Mz',
                         'F_wx': 'x', 'F_wy': 'y', 'F_wz': 'z',
                         'F_phix': 'rx', 'F_thetay': 'ry', 'F_thetaz': 'rz',
                         'F_B': 'B','F_Tw': 'Tw'},
                        axis=1, inplace=True)
        #
        return df_membf
    #
    #
    def dfmend(self, dftemp):
        """ """
        hforce = ['node_name', *self._plane.hforce]
        # get df 
        db = DBframework()
        header: list[str] = ['load_name', 'load_level',
                             #'load_number', 'load_title',
                             'load_system', 'element_name',
                             *hforce]
        df_nforce = db.DataFrame(data=dftemp, columns=header, index=None)
        return df_nforce[['load_name',
                          #'load_number','load_title',
                          'load_level',
                          'load_system', 'element_name',
                          *hforce]]    
    #
    #
    def get_reactions(self, nforce):
        """ get nodal reactions """
        #
        nforce = (nforce.groupby(['load_name',
                                  #'load_number',
                                  'load_level',
                                  'load_system',
                                  #'load_title',
                                  'node_name'],
                                  as_index=False)
                   [self._plane.hforce].sum())
        #
        supports = self._mesh._boundaries.supports()
        nname = list(supports)
        #
        nrsupp = nforce.loc[nforce['node_name'].isin(nname)]
        #
        #nrsupp = nforce.loc[nforce['node_name'].isin(nname)]
        #nrsupp2 =  nreacs.loc[nreacs['node_name'].isin(nname),
        #                      ['Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']].sum()
        # 
        #nreaction = (nrsupp.groupby(['load_name', 'load_number', 'load_type',
        #                              'load_system', 'load_title', 'node_name'])
        #             [self._plane.hforce].sum())
        # reset groupby to create new df
        #nreaction = nreaction.reset_index(names=['load_name', 'load_number', 'load_type',
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
        print("** Beam stress calculation")
        start_time = time.time()
        #
        elements = self._mesh._elements
        #
        members = bforce_df.groupby(['element_name',
                                     'load_name', 'load_level', 'load_system'])
        #
        for key, item in members:
            member = elements[key[0]]
            section = member.section
            stressdf = section.stress(df=item)
            conn = create_connection(self._db_file)
            with conn:            
                stressdf.to_sql('tb_Sbeam', conn,
                                if_exists='append', index=False)
        #
        uptime = time.time() - start_time
        print(f"** Time: {uptime:1.4e} sec")
    #
    # -----------------------------------------------------------
    # 
    #
    def results(self, Un, beam_steps:int = 10):
        """
        beam_steps : Integration points beam element (10 default)
        """
        print("** Postprocessing")
        start_time = time.time()
        #
        df_beamf, df_Qn = self.solve(Un, beam_steps)
        #
        #
        # combination
        load_comb = self._mesh._load.combination()
        df_comb = load_comb.to_basic() 
        #
        df_beam = update_memberdf(dfmemb=df_beamf,
                                  dfcomb=df_comb,
                                  values=['Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz',
                                          'x', 'y', 'z', 'rx', 'ry', 'rz',
                                          'B', 'Tw'])
        #
        df_Qbeam = df_beam[['load_name', 'load_level', 'load_system',
                            'element_name', 'node_end',
                            'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz',
                            'B', 'Tw']]
        #
        df_Dbeam = df_beam[['load_name', 'load_level', 'load_system',
                            'element_name', 'node_end',
                            'x', 'y', 'z', 'rx', 'ry', 'rz']]
        # Push to sql
        conn = create_connection(self._db_file)
        with conn:        
            df_Qbeam.to_sql('tb_Qbeam', conn,
                            if_exists='append', index=False)
            #
            df_Dbeam.to_sql('tb_Dbeam', conn, 
                            if_exists='append', index=False)
        #
        #
        # End node force combination
        df_Qn = update_memberdf(dfmemb=df_Qn,
                                dfcomb=df_comb,
                                values=self._plane.hforce)
        #
        #
        #header = ['load_name', 'load_level', 'system',
        #          'element_name', 'node_name']
        #header.extend(self._plane.hforce)
        #conn = create_connection(self._db_file)
        #with conn:
        #    df_Qn.to_sql('tb_Qm', conn,
        #                   index_label=header, 
        #                   if_exists='append', index=False)
        #
        #
        # node reacctions
        df_nforce, df_nreac = self.get_reactions(df_Qn)
        #
        # Push to sql
        #
        header = ['load_name', 'load_level', 'load_system',
                  'element_name', 'node_name']
        header.extend(self._plane.hforce)        
        #
        conn = create_connection(self._db_file)
        with conn:
            df_Qn.to_sql('tb_Fm', conn,
                         index_label=header, 
                         if_exists='append', index=False)
            #
            header.remove('element_name')
            df_nforce.to_sql('tb_Fn', conn,
                             index_label=header, 
                             if_exists='append', index=False)
            #
            df_nreac.to_sql('tb_Rn', conn,
                            index_label=header, 
                            if_exists='append', index=False)
        #
        #
        # displacments
        #df_ndisp = update_ndf(dfnode=Un, dfcomb=df_comb,
        #                      values=self._plane.hdisp)
        #
        #header = ['load_name',
        #          #'load_number',
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
        uptime = time.time() - start_time
        print(f"** Finish Time: {uptime:1.4e} sec")
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
    def push_QfXX(self, conn, data) -> None:
        """ push element forces to sql"""
        #head = 'load_name, load_level, system, element_name, node_end,'
        #cur = conn.cursor()
        #table = 'INSERT INTO tb_Qf( ' + head + ','.join(self._plane.hforce) + ','
        #table += ','.join(self._plane.hdisp) + ')'
        #vals = ['?' for item in self._plane.hforce]
        #vals.extend(['?' for item in self._plane.hdisp])
        #table += ' VALUES(?,?,?,?,?,' + ','.join(vals) + ')'
        table = 'INSERT INTO tb_Qf( load_name, load_level, system, element_name, node_end,\
                Fx, Fy, Fz, Mx, My, Mz, x, y, z, rx, ry, rz) \
                VALUES(?,?,?,?,?, ?,?,?,?,?,?, ?,?,?,?,?,?)'
        # push
        cur = conn.cursor()
        cur.executemany(table, data)
        print('ok')
    #
    #
    def push_QmXX(self, conn, data):
        """push element end forces to sql"""
        head = 'load_name, load_level, system, element_name, node_name,'
        table = 'INSERT INTO tb_Qm( ' + head + ','.join(self._plane.hforce) + ')'
        vals = ['?' for item in self._plane.hforce]
        table += ' VALUES(?,?,?,?,?,' + ','.join(vals) + ')'
        # push
        cur = conn.cursor()
        cur.executemany(table, data)
    #
    # -----------------------------------------------------------
    #
    #def end_force(self, df_ndisp, steps:int = 10):
    #    """Node force global system"""
    #    #
    #    if self._node_force:
    #        df_nforce = self._node_force
    #    else:
    #        df_membf, df_nforce = self.__call__(df_ndisp, steps)
    #        
    #    
    #    #beamdf = self.node_end_force(elements=self._elements,
    #    #                             df_ndisp=df_ndisp)
    #    # TODO : include additional element types
    #    #
    #    return df_nforce
    #
    #
    def beam_end_forceX(self, df_ndisp):
        """Node force global system"""
        elements = self._mesh._elements
        beamdf = self.node_end_force(elements=elements,
                                     df_ndisp=df_ndisp)
        # TODO : include additional element types
        #
        return beamdf    
    #
    def node_end_forceX(self, elements, df_ndisp):
        """
        Beam's end nodes force global system
        """
        dummyf = np.array([0]*self._plane.ndof)
        dispgrp = df_ndisp.groupby(['load_name', 'load_level',
                                    'load_number', 'load_title',
                                    'load_system'])
        #
        hdisp = ['node_name', *self._plane.hdisp]
        ndof = self._plane.ndof
        #
        dftemp = []
        for key, item in dispgrp:
            df1 = item[hdisp]
            df1.set_index('node_name', inplace=True)
            # start element big loop
            for mname, element in elements.items():
                nodes = element.connectivity
                # ---------------------------------------------
                # get beam end-node displacement in global system
                gndisp = np.concatenate((df1.loc[nodes[0]],
                                         df1.loc[nodes[1]]), axis=None)
                #
                # ---------------------------------------------
                # convert beam end-node disp to force [F = Kd] in global system
                #gnforce = mtxmul(element.K, gndisp)
                gnforce = element.K @ gndisp
                # ---------------------------------------------
                #
                dftemp.append([*key, mname, nodes[0], *gnforce[:ndof]])
                dftemp.append([*key, mname, nodes[1], *gnforce[ndof:]])
        #
        df_nforce = self.dfmend(dftemp=dftemp)
        return df_nforce
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
    def __init__(self, force:list, disp: list,
                 ndof: int, m2D: bool) -> None:
        self.force = force
        self.disp = disp
        self.ndof = ndof
        self.m2D = m2D
    #
    def R0(self, bload: list|None = None) -> tuple:
        """
        Axial   [FP, blank, blank, Fu]
        Bending [V, M, theta, w]
        Torsion [T, B, Psi, Phi, Tw]
        """
        bload0 = [0] * self.ndof
        try:
            if bload:
                bload0 = bload[:self.ndof]
        except (TypeError, ValueError):
            bload0 = bload[:self.ndof]
        #
        lnforce0 = self.force[:self.ndof]
        lndisp0 = self.disp[:self.ndof]
        if self.m2D:
            R0x = [1 * (lnforce0[0] + bload0[0]),
                   0, 0,
                   1 * lndisp0[0]] # Fx,0,0,Fu
            #
            R0y = [-1 * (lnforce0[1] - bload0[1]), # Vy
                   1 * (lnforce0[2] - bload0[2]), # Mz
                   -1 * lndisp0[2],  # thetaz
                   1 * lndisp0[1]]  # wy
            #
            R0z = [0, 0, 0, 0] # 0,0,0,0
            R0t = [0, 0, 0, 0, 0] # 0,0,0,0
        else:
            R0x = [lnforce0[0] + bload0[0],
                   0, 0, 1 * lndisp0[0]]
            #
            R0y = [-1 * (lnforce0[1] - bload0[1]),
                   1 * lnforce0[5] - bload0[5],
                   -1 * lndisp0[5],
                   1 * lndisp0[1]]
            #
            R0z = [-1 * (lnforce0[2] - bload0[2]),
                   1 * lnforce0[4] - bload0[4],
                   -1 * lndisp0[4],
                   1 * lndisp0[2]]
            #
            R0t = [1 * (lnforce0[3] + bload0[3]),
                   0, 0, 1 * lndisp0[3], 0]
        return Req(R0x, R0t, R0y, R0z)
#
#
class Req(NamedTuple):
    """ """
    x: list
    t: list
    y: list
    z: list
#
#
#
# -----------------------------------------------------------
# Combination process
# -----------------------------------------------------------
#
#
def update_memberdf(dfmemb, dfcomb, 
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
            comb['load_number'] = row.load_number
            comb['load_title'] = row.load_title
            #
            try:
                dftemp = db.concat([dftemp, comb], ignore_index=True)
            except UnboundLocalError:
                dftemp = comb
    #
    try:
        dftemp = dftemp.groupby(['load_name', 'load_number','load_level',
                                 'load_title', 'load_system',
                                 'element_name' ,'node_name'],
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
#def update_ndf(dfnode, dfcomb, 
#               values:list[str]): #['x', 'y', 'z', 'rx', 'ry', 'rz']
#    """
#    Update node displacements to include lcomb
#    """
#    db = DBframework()
#    # group basic load by name
#    ndgrp = dfnode.groupby('load_name')
#    # get combinations with basic loads 
#    #
#    combgrp = dfcomb.groupby('load_name')
#    for key, combfactors in combgrp:
#        for row in combfactors.itertuples():
#            comb = ndgrp.get_group(row.basic_load).copy()
#            comb.loc[:, values] *= row.factor
#            comb['load_level'] = 'combination'
#            comb['load_name'] = row.load_name
#            comb['load_number'] = row.load_number
#            comb['load_title'] = row.load_title
#            #
#            try:
#                dftemp = db.concat([dftemp, comb], ignore_index=True)
#            except UnboundLocalError:
#                dftemp = comb
#    #
#    try:
#        #check = dftemp.groupby(['node_name', 'c']).sum().reset_index()
#        dftemp = dftemp.groupby(['load_name', 'load_number','load_level',
#                                 'load_title', 'load_system','node_name'],
#                                  as_index=False)[values].sum()
#        #test
#        dfnode = db.concat([dfnode, dftemp], ignore_index=True)
#    except UnboundLocalError:
#        pass
#    #
#    return dfnode #, memb_comb
#
#

