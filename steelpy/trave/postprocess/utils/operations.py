# 
# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
from collections import namedtuple
from dataclasses import dataclass
import time
from typing import NamedTuple
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
    __slots__ = ['_plane', '_mesh', '_result_name']
    
    def __init__(self, mesh, result_name:int|str,) -> None:
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
        #
        if Pdelta:
            load = self._mesh._load._combination
            Un_set = Un_grp.get_group(('combination', ))
        else:
            load = self._mesh._load._basic
            Un_set = Un_grp.get_group(('basic', ))  
        #
        Fbeam = self.solve_beam_end_disp(elements, Un_set, Pdelta)        
        #
        df_beamf, df_nforce = self.solve_beam_forces(elements,
                                                     Fb=Fbeam,
                                                     load=load,
                                                     steps=steps,
                                                     Pdelta=Pdelta)
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
        #ndof:int = 6
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
                Tb = element.T3D()
                #
                # ---------------------------------------------
                # beam end-nodes global displacement {U}
                Un_global = np.concatenate((Un_set.loc[nodes[0]],
                                            Un_set.loc[nodes[1]]),
                                           axis=None)
                #
                # ---------------------------------------------
                # Convert global end-node disp in beam's local system
                # 
                Un_local = Tb @ Un_global
                #
                # --------------------------------------------
                # get beam [K] global
                Kb = getattr(element, Kglobal)(plane2D, Un_global)
                #
                # ---------------------------------------------
                # convert beam end-node disp to force ({F} = [K]{U}) in global system
                FU_global = Kb @ Un_global
                #
                # ---------------------------------------------
                # Calculate beam end-nodes force in local system
                #Fb_local = Tb @ FU_global
                FU_local = Kb @ Un_local
                #
                # ---------------------------------------------
                # Build list with results
                # [load_name, mesh_name, load_level, system, element_name,
                # node_name, Un_local, FU_local, FU_global]
                #Fb_temp.append([*key, e_name, nodes[0],
                #                Un_local[:ndof], FU_local[:ndof], FU_global[:ndof]])
                #Fb_temp.append([*key, e_name, nodes[1],
                #                Un_local[ndof:], FU_local[ndof:], FU_global[ndof:]])
                #
                Fb_temp.append([*key, e_name, nodes,
                                Un_local, FU_local, FU_global])                
        #
        # ---------------------------------------------
        #Fb_temp = self._mf2df(Fb=Fb_temp)
        #
        db = DBframework()
        header = ['load_name', 'mesh_name',
                  'load_level',
                  'element_name', 'nodes',
                  'Un_local', 'FU_local', 'FU_global']
        
        Fb_temp = db.DataFrame(data=Fb_temp, columns=header, index=None)        
        #
        return Fb_temp
    #
    def solve_beam_forces(self, elements,
                          Fb:DBframework.DataFrame,
                          load:DBframework.DataFrame,
                          #beam_fer:DBframework.DataFrame,
                          steps: int, Pdelta: bool)-> (DBframework.DataFrame, DBframework.DataFrame):
        """
        Beam's internal forces along its lenght
        """
        #1/0
        start_time = time.time()
        #
        beam_fer = load.FER_ENL()
        #
        ndof = 6 # self._plane.ndof
        #
        Fb_grp = Fb.groupby(['load_name', 'mesh_name',
                             'load_level'])
        #
        bfer_grp = beam_fer.groupby(['load_name', 'mesh_name',
                                     'load_level'])
        #
        #
        header = ['element_name', 'nodes', 'Un_local', 'FU_local', 'FU_global']
        #
        #plane2d = self._plane.plane2D
        dftemp: list = []
        mif_data: list = []
        for key, item in Fb_grp:
            # get node displacement variables
            Fb_set = item[header].set_index('element_name')
            #
            for beam_item in Fb_set.itertuples():
                bname = beam_item.Index
                nodes = beam_item.nodes
                #
                element = elements[bname]
                #Tb = element.T(plane2d)
                Tb = element.T3D()
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
                #
                #TODO: confirm change reactions sign
                Fbeq = FbEqMods(beam_item.Un_local,
                                self._plane, Pdelta)
                #
                # ---------------------------------------------
                # Beam load (udl/point)
                try:                   
                    # get FER beam end load
                    beam_fer_enl = bfer_grp.get_group(key)
                    beam_fer_enl = beam_fer_enl.groupby(['element_name'])                    
                    beam_fer_enl = beam_fer_enl.get_group((bname, ))
                    # calculate beam end 1 reactions
                    R0 = self._beam_R0(nodes=nodes,
                                       Tb=Tb,
                                       FUn=beam_item.FU_global, 
                                       fer_enl=beam_fer_enl,
                                       eq=Fbeq)
                    #
                    # get load function
                    load_item = load[key[0]]
                    #print(f'P0 = {Fbeq.P0()}')
                    #bload = load_item.beam()[bname]
                    bload_fuction = load_item.function(beam_name=bname,
                                                       load_name=key[0], 
                                                       steps=steps,
                                                       Pa=Fbeq.P0())
                    # ['load_name', 'mesh_name',
                    #  'load_comment', 'load_type',
                    #  'load_level', 'load_system',
                    #  'element_name', 'length',
                    #  'axial', 'torsion', 'VM_inplane', 'VM_outplane']
                    #
                    #print('---->')
                    lbforce = [['local', bname, bstep.length,
                                *beam.response(x=bstep.length,
                                               R0=[R0.x, R0.t, R0.y, R0.z],
                                               Fx=[bstep.axial, bstep.torsion,
                                                   bstep.VM_inplane, bstep.VM_outplane])]
                               for bstep in bload_fuction.itertuples()]
                    #
                except KeyError :
                    # --------------------------------------------
                    #1/0
                    #print('---->')
                    lbforce = self._beam_no_load(bname,
                                                 beam, beam_item.FU_global,
                                                 Fbeq, Tb, 
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
                gn_force = beam_item.FU_global
                #
                # ---------------------------------------------
                #
                dftemp.append([*key, 'global', bname, nodes[0], *gn_force[:ndof]])
                dftemp.append([*key, 'global', bname, nodes[1], *gn_force[ndof:]])                
                #
                # ---------------------------------------------
                # Member total force in local system
                #
                # Axial   [FP, blank, blank, Fu]
                # Torsion [T, Phi, Psi, B, Tw] - [T, theta, theta1, theta2, theta3]
                # Bending [V, M, theta, w]
                #
                mif_data.extend([[*key[:3], # load_name, component_name, load_level,
                                  *lbf[:3], # load_system, element_name, length
                                  *lbf[3],  # Axial
                                  *lbf[4],  # Torsion
                                  *lbf[5],  # Bending in plane
                                  *lbf[6]]  # Bending out plane
                                 for lbf in lbforce])
        #
        # ---------------------------------------------
        #
        df_nforce = self._mf2df(Fb=dftemp)
        # Member internal forces
        df_mif = self._mif2df(mif_data)
        #
        # ---------------------------------------------
        uptime = time.time() - start_time
        print(f"** Beam Force Calculation: {uptime:1.4e} sec")
        return df_mif, df_nforce
    #   
    #
    # -----------------------------------------------------------
    # Operations
    #
    def _beam_R0(self,
                 nodes: list,
                 Tb, FUn,
                 fer_enl, eq):
        """ """
        #
        # [load_name, member_load_title, load_type, load_system, 
        # beam_number, x, Fx, Fy, Fz]
        #
        #hforce = [*self._plane.hforce]
        hforce = ['Fx','Fy','Fz','Mx','My','Mz']
        tforce = ['Psi', 'B', 'Tw'] # torsion part       
        #
        # Torsion
        mtload = fer_enl.groupby(['node_name'],
                                 as_index=False)[tforce].sum()
        mtload.set_index('node_name', inplace=True)                    
        #
        # Get beam end loads
        #
        fer_enl = fer_enl.groupby(['node_name'],
                                  as_index=False)[hforce].sum()
        
        fer_enl.set_index('node_name', inplace=True)                 
        #
        FU_fer = np.concatenate((fer_enl.loc[nodes[0]],
                                 fer_enl.loc[nodes[1]]), axis=None)
        #
        # Calculate total end node force in local system
        Pu = Tb @ (FUn - FU_fer)
        # get reactions with correct sign for postprocessing
        #
        R0 = eq.R0(bm_load=Pu,
                   tload=mtload.loc[nodes[0]])
        return R0
    #    
    def _beam_no_load(self, beam_name: int|str,
                      beam, FUn, eq, Tb, 
                      geometry, steps: int):
        """ """
        
        # Dummy Bending [V, M, theta, w]
        Fblank = [0, 0, 0, 0]
        # Dummy Torsion [T, Phi, Psi, B, Tw] - [T, theta, theta1, theta2, theta3]
        Fblank_t = [0, 0, 0, 0, 0]
        #
        # [Fx,Fy,Fz,Mx,My,Mz]
        Pu = Tb @ FUn
        R0 = eq.R0(bm_load=Pu,
                   tload=Tload(0, 0, 0))
        #
        # [beam_number, load_title, x, Fx, Fy, Fz, Mx, My, Mz]
        #
        Lsteps = linstep(d=geometry.d,
                         L=beam.L, steps=steps)
        #
        lbforce = [['local', beam_name,  xstep,
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
#@dataclass
class FbEqMods(NamedTuple):
    """ Local system"""
    Fb_local: list
    plane: bool
    Pdelta: bool
    #
    def P0(self) -> float:
        """ """
        # [Fx, Fy, Fz, Mx, My, Mz]
        #Fx = self.Fb_local[0]
        if self.Pdelta:
            return self.Fb_local[0]
        return 0.0
    #
    def R0(self, bm_load: list, tload: tuple|list) -> tuple:
        """
        Axial   [FP, blank, blank, Fu]
        Bending [V, M, theta, w]
        Torsion [T, Phi, Psi, B, Tw]
        """
        # select node disp end 0
        nd0 = self.Fb_local[:self.plane.ndof]
        bload = bm_load[:self.plane.ndof]
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
        # force In plane
        #R0y[0] *= -1   # Vy
        R0y[1] *= -1  # Mz
        # displacement
        #R0y[2] *= -1  # rz
        R0y[3] *= -1   # dy
        #
        # force Out plane
        #R0z[0] *= -1   # Vz
        #R0z[1] *= -1   # My
        # displacement
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

