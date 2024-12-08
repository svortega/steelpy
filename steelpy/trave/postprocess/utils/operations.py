# 
# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
from collections import namedtuple
from dataclasses import dataclass
import time
#from typing import NamedTuple
#from datetime import datetime as dt
#
#
# package imports
from steelpy.trave.beam.main import BeamBasic
#
from steelpy.utils.math.operations import linstep
from steelpy.utils.dataframe.main import DBframework
#from steelpy.utils.sqlite.utils import create_connection, create_table
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
#
# -----------------------------------------------------------
#  Elements
# -----------------------------------------------------------
#
#
#
@dataclass
class MainPostProcess:
    __slots__ = ['_plane', '_mesh', '_result_name'] #'_name',
    
    def __init__(self, mesh, result_name:int|str,) -> None: #, name: str
        """
        """
        self._mesh = mesh
        self._plane = mesh._plane
        self._result_name = result_name
    #
    # -----------------------------------------------------------
    # Operations Node, Element forces and displacement
    # -----------------------------------------------------------
    #
    def solve(self, Un, steps: int,
              Pdelta: bool):
        """
        Un : node global displacement result
        steps: Beam load locations (?)
        Pdelta: 2nd order flag
        """
        elements = self._mesh._elements
        Un_grp = Un.df.groupby(['load_level'])
        #
        if Pdelta:
            Un_set = Un_grp.get_group(('combination', ))
            Fb_local = self.solve_beam_end_disp(elements, Un_set, Pdelta)
            #
            comb_load = self._mesh._load._combination
            beam_enl = comb_load.FER_ENL()
            load_func = comb_load.function(steps=steps,
                                           Fb_local=Fb_local)
        else:
            Un_set = Un_grp.get_group(('basic', ))
            # Solve node and element foce and displacements
            basic_load = self._mesh._load._basic
            beam_enl = basic_load.FER_ENL()
            load_func = basic_load.function(steps=steps,
                                            Pa=0.0, factor=1.0)    
        #
        #
        df_beamf, df_nforce = self.solve_beam_forces(elements, Un_set,
                                                     load_function=load_func,
                                                     benl=beam_enl,
                                                     steps=steps, Pdelta=Pdelta)
        #
        return df_beamf, df_nforce
    #
    # -----------------------------------------------------------
    #
    def solve_beam_end_disp(self, elements,
                            Un:DBframework.DataFrame,
                            Pdelta: bool)-> DBframework.DataFrame:
        """
        Convert beam end-node disp to force [F = Kd] in global system

        Return:
        Fb_local : ?
        Fb_global : ?
        """
        # Setup 
        ndof:int = 6
        plane2D:bool = False
        Kglobal:str = 'Ke'
        if Pdelta:
            Kglobal:str = 'Kt'
        #
        Un_grp = Un.groupby(['load_name', 'mesh_name', 'load_level'])
        head_disp = ['node_name', 'x', 'y', 'z', 'rx', 'ry', 'rz']
        #
        Fb_temp: list = []
        for key, Un_item in Un_grp:
            Un_set = Un_item[head_disp].set_index('node_name')
            #
            for e_name, element in elements.items():
                nodes = element.connectivity
                #
                # ---------------------------------------------
                # beam end-nodes global displacement {U}
                nd_global = np.concatenate((Un_set.loc[nodes[0]],
                                            Un_set.loc[nodes[1]]),
                                           axis=None)
                #
                # --------------------------------------------
                # get beam [K] global
                Kb = getattr(element, Kglobal)(plane2D, nd_global)
                #
                # ---------------------------------------------
                # convert beam end-node disp to force ({F} = [K]{U}) in global system
                Fb_global = Kb @ nd_global
                #
                # ---------------------------------------------
                # Calculate beam end-nodes force in local system
                Tb = element.T3D()
                Fb_local = Tb @ Fb_global
                #
                # ---------------------------------------------
                # Build list with results
                Fb_temp.append([*key, 'local', e_name, nodes[0], *Fb_local[:ndof]])
                Fb_temp.append([*key, 'local', e_name, nodes[1], *Fb_local[ndof:]])
        #
        # ---------------------------------------------
        Fb_local = self._mf2df(Fb=Fb_temp)
        return Fb_local
    #
    def solve_beam_forces(self, elements,
                          Un:DBframework.DataFrame,
                          load_function:DBframework.DataFrame,
                          benl:DBframework.DataFrame,
                          steps: int, Pdelta: bool)-> (DBframework.DataFrame, DBframework.DataFrame):
        """
        Beam's internal forces along its lenght
        """
        #1/0
        start_time = time.time()
        #
        Klocal = 'Ke_local'
        Kglobal = 'Ke'
        if Pdelta:
            Klocal = 'Kt_local'
            Kglobal = 'Kt'
        #
        Un_grp = Un.groupby(['load_name', 'mesh_name',
                             'load_level',  'system'])
        #
        benlgrp = benl.groupby(['load_name', 'mesh_name',
                                'load_level', 'load_system'])
        #
        lf_grp = load_function.groupby(['load_name', 'mesh_name',
                                       'load_level'])
        #
        #
        ndof = self._plane.ndof
        hdisp = ['node_name', *self._plane.hdisp]
        #
        plane2d = self._plane.plane2D
        dftemp: list = []
        mif_data: list = []
        for key, node_item in Un_grp:
            # get node displacement variables
            Un_disp = node_item[hdisp].set_index('node_name')
            #
            for mname, element in elements.items():
                nodes = element.connectivity
                Tb = element.T(plane2d)
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
                # ---------------------------------------------
                # Select beam end-nodes displacement
                #
                ben_dglobal = np.concatenate((Un_disp.loc[nodes[0]],
                                              Un_disp.loc[nodes[1]]), axis=None)
                #
                # ---------------------------------------------
                # Convert global end-node disp in beam's local system
                # 
                ben_dlocal = Tb @ ben_dglobal
                #
                # --------------------------------------------
                # Convert beam end-node disp to force [F = Kd] in global system
                #
                k_local = getattr(element, Klocal)(plane2d, ben_dglobal)
                FUtn = k_local @ ben_dlocal
                #
                #TODO: confirm change reactions sign
                eq = NodeGenRespEq(ben_dlocal, self._plane)
                #
                # ---------------------------------------------
                # Beam load (udl/point)
                try:
                    #
                    mbload = lf_grp.get_group(key[:3])
                    benl = benlgrp.get_group(key)
                    lbforce = self._beam_load(mname, nodes,
                                              beam, Tb, FUtn,
                                              mbload, benl, eq)
                    ##
                    ## [load_name, member_load_title, load_type, load_system, 
                    ## beam_number, x, Fx, Fy, Fz]
                    ##
                    #mbload = lf_grp.get_group(key[:3])
                    #mbload = mbload.groupby(['element_name'])                    
                    #mnload = mbload.get_group((mname, ))
                    #mnload = mnload.groupby(['length'],
                    #                        as_index=False)[hload].sum()
                    ##
                    #benl = benlgrp.get_group(key)
                    #benl = benl.groupby(['element_name'])                    
                    #bnload = benl.get_group((mname, ))
                    ##
                    ## Torsion
                    #mtload = bnload.groupby(['node_name'],
                    #                        as_index=False)[tforce].sum()
                    #mtload.set_index('node_name', inplace=True)                    
                    ##
                    ## Get beam end loads
                    ##
                    #nf_local = bnload.groupby(['node_name'],
                    #                           as_index=False)[hforce].sum()
                    #nf_local.set_index('node_name', inplace=True)                 
                    ##
                    #FUen = np.concatenate((nf_local.loc[nodes[0]],
                    #                       nf_local.loc[nodes[1]]), axis=None)
                    ##
                    #Pu = Tb @ FUen - FUtn
                    ##
                    ## get reactions with correct sign for postprocessing
                    ##
                    #R0 = eq.R0(bload=Pu,
                    #           tload=mtload.loc[nodes[0]])
                    ##
                    #lbforce = [['local', mname, bstep.length,
                    #            *beam.response(x=bstep.length,
                    #                           R0=[R0.x, R0.t, R0.y, R0.z],
                    #                           Fx=[*bstep[2:]])]
                    #           for bstep in mnload.itertuples()]
                
                except KeyError :
                    # --------------------------------------------
                    #
                    lbforce = self._beam_no_load(mname,
                                                 beam, FUtn, eq, 
                                                 element.section.geometry,
                                                 steps)
                    ##
                    ## [Fx,Fy,Fz,Mx,My,Mz]
                    #Pu = -1 * FUtn
                    #R0 = eq.R0(bload=Pu,
                    #           tload=Tload(0, 0, 0))
                    ##
                    ## [beam_number, load_title, x, Fx, Fy, Fz, Mx, My, Mz]
                    ##
                    #Lsteps = linstep(d=element.section.geometry.d,
                    #                 L=element.L, steps=steps)
                    ##
                    #lbforce = [['local', mname,  xstep,
                    #            *beam.response(x=xstep, R0=[R0.x, R0.t, R0.y, R0.z],
                    #                           Fx=[Fblank, Fblank_t, Fblank, Fblank])]
                    #           for xstep in Lsteps]
                #
                # ---------------------------------------------
                # convert beam end-node disp to force [F = Kd] in global system
                #
                k_global = getattr(element, Kglobal)(plane2d, ben_dglobal)
                gnforce = k_global @ ben_dglobal
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
        df_nforce = self._mf2df(Fb=dftemp,
                                head_force=self._plane.hforce)
        # Member internal forces
        df_mif = self._mif2df(mif_data)
        #
        # ---------------------------------------------
        uptime = time.time() - start_time
        print(f"** Beam Force Calculation: {uptime:1.4e} sec")
        return df_mif, df_nforce
    #
    # -----------------------------------------------------------
    # Operations
    #
    def _beam_load(self, mname: int|str, nodes: list,
                   beam, Tb, FUtn,
                   mbload, benl, eq):
        """ """
        #
        # [load_name, member_load_title, load_type, load_system, 
        # beam_number, x, Fx, Fy, Fz]
        #
        hforce = [*self._plane.hforce]
        tforce = ['Psi', 'B', 'Tw'] # torsion part
        hload = ['axial', 'torsion', 'VM_inplane', 'VM_outplane']        
        #
        #mbload = lf_grp.get_group(key[:3])
        mbload = mbload.groupby(['element_name'])                    
        mnload = mbload.get_group((mname, ))
        mnload = mnload.groupby(['length'],
                                as_index=False)[hload].sum()
        #
        #benl = benlgrp.get_group(key)
        benl = benl.groupby(['element_name'])                    
        bnload = benl.get_group((mname, ))
        #
        # Torsion
        mtload = bnload.groupby(['node_name'],
                                as_index=False)[tforce].sum()
        mtload.set_index('node_name', inplace=True)                    
        #
        # Get beam end loads
        #
        nf_local = bnload.groupby(['node_name'],
                                   as_index=False)[hforce].sum()
        nf_local.set_index('node_name', inplace=True)                 
        #
        FUen = np.concatenate((nf_local.loc[nodes[0]],
                               nf_local.loc[nodes[1]]), axis=None)
        #
        Pu = Tb @ FUen - FUtn
        #
        # get reactions with correct sign for postprocessing
        #
        R0 = eq.R0(bload=Pu,
                   tload=mtload.loc[nodes[0]])
        #
        lbforce = [['local', mname, bstep.length,
                    *beam.response(x=bstep.length,
                                   R0=[R0.x, R0.t, R0.y, R0.z],
                                   Fx=[*bstep[2:]])]
                   for bstep in mnload.itertuples()]
        return lbforce
    #
    def _beam_no_load(self, mname: int|str,
                      beam, FUtn, eq, 
                      geometry, steps: int):
        """ """
        
        # Dummy Bending [V, M, theta, w]
        Fblank = [0, 0, 0, 0]
        # Dummy Torsion [T, Phi, Psi, B, Tw] - [T, theta, theta1, theta2, theta3]
        Fblank_t = [0, 0, 0, 0, 0]
        #
        # [Fx,Fy,Fz,Mx,My,Mz]
        Pu = -1 * FUtn
        R0 = eq.R0(bload=Pu,
                   tload=Tload(0, 0, 0))
        #
        # [beam_number, load_title, x, Fx, Fy, Fz, Mx, My, Mz]
        #
        Lsteps = linstep(d=geometry.d,
                         L=beam.L, steps=steps)
        #
        lbforce = [['local', mname,  xstep,
                    *beam.response(x=xstep, R0=[R0.x, R0.t, R0.y, R0.z],
                                   Fx=[Fblank, Fblank_t, Fblank, Fblank])]
                   for xstep in Lsteps]
        return lbforce
    #
    # -----------------------------------------------------------
    # DataFrame Operations
    #
    def _mif2df(self, member_load:list) -> DBframework.DataFrame:
        """
        Return:
            Beams' internal forces dataframe
        """
        header = ['load_name', 'mesh_name', 'load_level',
                  'system',
                  'element_name', 'length',
                  'Fx', 'blank1', 'blank2', 'x',  # axial
                  'Mx', 'rx', 'Psi',  'B',  'Tw', # torsion
                  'Fy', 'Mz', 'rz', 'y',          # bending in plane
                  'Fz', 'My', 'ry', 'z']          # bending out plane
        #
        db = DBframework()
        df_membf = db.DataFrame(data=member_load,
                                columns=header, index=None)
        # reorder columns
        df_membf['result_name'] = self._result_name
        df_membf.drop(columns=['blank1', 'blank2'])
        return df_membf[['mesh_name', 'result_name',
                         'load_name', 'load_level',
                         'element_name', 'length',
                         'system',
                         'Fx', 'Fy', 'Fz',
                         'Mx', 'My', 'Mz',
                         'x', 'y', 'z',
                         'rx', 'ry', 'rz',
                         'Psi', 'B', 'Tw']]
    #
    #
    def _mf2df(self, Fb:list[list],
               head_force:list[str]|None=None)-> DBframework.DataFrame:
        """
        Return:
            Beams' forces dataframe
        """
        if not head_force:
            head_force = ['Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']
        #
        db = DBframework()
        header = ['load_name', 'mesh_name',
                  'load_level', 'system',
                  'element_name', 'node_name', *head_force]
        mf_df = db.DataFrame(data=Fb, columns=header, index=None)
        #
        mf_df['result_name'] = self._result_name
        header = ['mesh_name', 'result_name',
                  'load_name', 'load_level',
                  'element_name', 'node_name',
                  'system', *head_force]
        return mf_df[header]
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
        #R0y[0] *= -1   # Vy
        R0y[1] *= -1  # Mz
        R0y[2] *= -1  # rz
        #R0y[3] *= -1   # dy
        #
        # Out plane
        R0z[0] *= -1   # Vz
        R0z[1] *= -1   # My
        R0z[2] *= -1   # ry
        R0z[3] *= -1   # dz
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
        T0 = -1 * bload[3]
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
Tload = namedtuple('TorsionLoad', ['Psi', 'B', 'Tw'])
#
# -----------------------------------------------------------
# Combination utils
# -----------------------------------------------------------
#
#
#

