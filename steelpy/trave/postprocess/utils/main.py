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
    __slots__ = ['_plane', '_mesh', '_result_name', '_results']
    
    def __init__(self, mesh, result_name:int|str, 
                 results) -> None:
        """
        """
        self._mesh = mesh
        self._plane = mesh._plane
        self._result_name = result_name
        #self._db_file = db_file
        self._results = results
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
        Un_grp = Un.groupby(['load_level'])
        #
        #
        if Pdelta:
            load = self._mesh._load._combination
            Un_set = Un_grp.get_group(('combination', ))
        else:
            load = self._mesh._load._basic
            Un_set = Un_grp.get_group(('basic', ))  
        #
        Fb_end = self.solve_beam_end(elements, Un_set, Pdelta)
        #
        Fb_result, Qn = self.solve_beam_forces(elements,
                                               Fb=Fb_end,
                                               load=load,
                                               steps=steps,
                                               Pdelta=Pdelta)
        #
        # ---------------------------------
        # Calculate node reactions
        Qn = self.solve_node_reactions(Fb_end=Qn,
                                       load=load)
        #
        # ---------------------------------
        # load comb update
        if not Pdelta:
            # combination
            combination = self._mesh._load.combination()
            comb2basic = combination.to_basic()
            Fb_result, Qn = self._add_comb(Fb_result, Qn, comb2basic)
        #
        # ---------------------------------
        # Process
        # element force
        self._results._push_beam_result(Fb_result)
        # Node
        self._results._push_node_reaction(Qn)
        #
        return Fb_result #, Qn
    #
    # -----------------------------------------------------------
    # update load
    def _add_comb(self, beam_force, Qn, lcomb):
        """
        Update load with combinations

        Qn:
        beam_force:
        lcomb:
        """
        #db = DBframework()
        #
        values = ['Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz',
                  'x', 'y', 'z', 'rx', 'ry', 'rz',
                  'Psi', 'B', 'Tw']
        item = ['mesh_name', 'result_name',
                'load_name', 'load_level',
                'element_name', 'length', 'system']
        beam_force = update_combination(member=beam_force,
                                        combination=lcomb,
                                        item=item,
                                        values=values)
        #
        #try:
        #    bf_comb = bf_comb.groupby(item, as_index=False)[values].sum()
        #    beam_force = db.concat([beam_force, bf_comb], ignore_index=True)
        #except UnboundLocalError:
        #    pass
        #
        # End node force combination
        values = ['Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']
        item = ['mesh_name', 'result_name',
                'load_name', 'load_level',
                'node_name', 'system']
        Qn = update_combination(member=Qn,
                                combination=lcomb,
                                item=item,
                                values=values)
        #
        #
        #try:
        #    Qn_comb = Qn_comb.groupby(item, as_index=False)[values].sum()
        #    Qn = db.concat([Qn, Qn_comb], ignore_index=True)
        #except UnboundLocalError:
        #    pass
        #
        return beam_force, Qn
    #
    # -----------------------------------------------------------
    #
    def solve_beam_end(self, elements,
                       Un:DBframework.DataFrame,
                       Pdelta: bool)-> DBframework.DataFrame:
        """
        Convert beam end-node disp to force [F = Kd] in global system

        Return:
        Fb_local : ?
        Fb_global : ?
        """
        # Setup
        plane2D:bool = False
        Kglobal:str = 'Ke'
        if Pdelta:
            Kglobal:str = 'Kt'
        #
        Un_grp = Un.groupby(['load_name', 'mesh_name', 'load_level'])
        head_disp = ['node_name', 'x', 'y', 'z', 'rx', 'ry', 'rz']
        #head_disp = ['node_name', *self._mesh._plane.hdisp]
        #
        Fb_temp: list = []
        for key, Un_item in Un_grp:
            Un_set = Un_item[head_disp].set_index('node_name')
            #
            for e_name, element in elements.items():
                nodes = element.connectivity
                Tb = element.T3D()
                # ---------------------------------------------
                # beam end-nodes global displacement {U}
                Un_global = np.concatenate((Un_set.loc[nodes[0]],
                                            Un_set.loc[nodes[1]]),
                                           axis=None)
                # --------------------------------------------
                # get beam [K] global
                Kb = getattr(element, Kglobal)(plane2D, Un_global)
                # ---------------------------------------------
                # convert beam end-node disp to force in global system
                # {F} = [K]{U}
                Fb_global = Kb @ Un_global
                # ---------------------------------------------
                # Convert global end-node disp in beam's local system
                Un_local = Tb @ Un_global
                # ---------------------------------------------
                # Calculate beam end-nodes force in local system
                Fb_local = Kb @ Un_local
                #Fb_local2 = Tb @ Fb_global
                # ---------------------------------------------
                # Build list with results
                # [load_name, mesh_name, load_level, system, element_name,
                #  node_name, Un_local, FU_local, FU_global]
                #
                Fb_temp.append([*key, e_name, nodes,
                                Un_local, Fb_local, Fb_global])
        #
        # ---------------------------------------------
        #
        header = ['load_name', 'mesh_name',
                  'load_level',
                  'element_name', 'nodes',
                  'Un_local', 'Fb_local', 'Fb_global']
        Fb_temp = list2df(data=Fb_temp, header=header)
        return Fb_temp
    #
    def solve_beam_forces(self, elements,
                          Fb:DBframework.DataFrame,
                          load:DBframework.DataFrame,
                          #beam_fer:DBframework.DataFrame,
                          steps: int, Pdelta: bool)-> (DBframework.DataFrame, DBframework.DataFrame):
        """
        Beam's internal forces along its length
        """
        #1/0
        start_time = time.time()
        #
        beam_fer = load.FER_ENL()
        #
        ndof = 6 # self._plane.ndof
        #
        head_grp = ['load_name', 'mesh_name',
                    'load_level']
        Fb_grp = Fb.groupby(head_grp)
        bfer_grp = beam_fer.groupby(head_grp)
        #
        #
        header = ['element_name', 'nodes',
                  'Un_local', 'Fb_local', 'Fb_global']
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
                    Pa, R0, Fbi = self._beam_R0(nodes=nodes,
                                               Tb=Tb,
                                               Fb_local=beam_item.Fb_local,
                                               fer_global=beam_fer_enl,
                                               eq=Fbeq)
                    #
                    gn_force = -1 * Tb @ Fbi
                    # get load function
                    load_item = load[key[0]]
                    #print(f'P0 = {Pa}')
                    #bload = load_item.beam()[bname]
                    bload_fuction = load_item.function(beam_name=bname,
                                                       load_name=key[0], 
                                                       steps=steps,
                                                       Pa=Pa)
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
                    #print('---->')
                except KeyError :
                    # --------------------------------------------
                    #print('---->')
                    lbforce = self._beam_no_load(beam_name=bname,
                                                 beam=beam,
                                                 Fb_local=beam_item.Fb_local,
                                                 eq=Fbeq,
                                                 geometry=element.section.geometry,
                                                 steps=steps)
                    #
                    gn_force = -1 * beam_item.Fb_global
                    #
                #
                # ---------------------------------------------
                # convert beam end-node disp to force [F = Kd] in global system
                #
                #gn_force = beam_item.Fb_global
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
                mif_data.extend([[*key[:3], # load_name, mesh_name, load_level,
                                  self._result_name,
                                  *lbf[:3], # load_system, element_name, length
                                  *lbf[3],  # Axial
                                  *lbf[4],  # Torsion
                                  *lbf[5],  # Bending in plane
                                  *lbf[6]]  # Bending out plane
                                 for lbf in lbforce])
        #
        # ---------------------------------------------
        #
        #
        # Member internal forces
        #
        Fb_result = mif2df2(mif_data)
        df_nforce = self._mf2df(Fb=dftemp)
        #
        # ---------------------------------------------
        uptime = time.time() - start_time
        print(f"** Beam Force Calculation: {uptime:1.4e} sec")
        return Fb_result, df_nforce
    #
    #
    def solve_stress(self, beam_force:DBframework.DataFrame):
        """get elements stress"""
        #print("** Beam stress calculation")
        start_time = time.time()
        #1 / 0
        beamfdf = beam_force.copy()
        #mesh_name=self._mesh._name
        elements = self._mesh._elements
        memb_grp = beamfdf.groupby(['mesh_name', 'result_name',
                                    'load_name',  'load_level',
                                    'system', 'element_name'])
        #
        for key, memb_item in memb_grp:
            #print(key)
            member = elements[key[-1].tolist()]
            section = member.section
            material = member.material
            beam_stress = section.stress(beam_result=memb_item,
                                         E=material.E,
                                         G=material.G,
                                         poisson=material.poisson)
            #
            self._results._push_beam_stress(beam_stress)
        #
        uptime = time.time() - start_time
        print(f"** Beam Stress Calculation: {uptime:1.4e} sec")
    #
    #
    # -----------------------------------------------------------
    #
    def solve_node_reactions(self, Fb_end, load):
        """ """
        # TODO : maybe there is a better pandas way that the code below..
        header = ['load_name', 'mesh_name', 'load_level', 'system']
        supports = self._mesh._boundaries.support()
        supports = list(supports)
        #
        Fb_end = Fb_end[Fb_end['node_name'].isin(supports)]
        # Fb_end = Fb_end.query('node_name != @supports')
        Fb_end_grp = Fb_end.groupby(header)
        #
        #node_load = load.node().force
        Pn = load.Pn()
        if not Pn.empty:
            Pn = Pn[Pn['node_name'].isin(supports)]
            # Pn = Pn.query('node_name != @supports')
            Pn_grp = Pn.groupby(header)
        else:
            Pn_grp = Pn
        #
        header += ['node_name']
        values = header + ['Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']
        temp = []
        for key, item in Fb_end_grp:
            Fb_set = item[values].set_index(header)
            # get nodal load if available
            try:
                Pn_set = Pn_grp.get_group(key)
                Pn_set = Pn_set[values].set_index(header)
                temp.append(Fb_set.add(Pn_set).reset_index())
            except (KeyError, AttributeError) :
                temp.append(Fb_set.reset_index())
            #print(key)
        #
        df = DBframework()
        Qn = df.concat(temp, ignore_index=True)
        Qn['result_name'] = self._result_name
        #Qn['system'] = 'global'
        #print('--->')
        return Qn
    #
    # -----------------------------------------------------------
    # Operations
    #
    def _beam_R0(self,
                 nodes: list,
                 Tb, Fb_local,
                 fer_global, eq):
        """ """
        force_header = ['Fx','Fy','Fz','Mx','My','Mz']
        torsion_header = ['Psi', 'B', 'Tw'] # torsion part
        #
        # Torsion
        bt_load = fer_global.groupby(['node_name'],
                                     as_index=False)[torsion_header].sum()
        bt_load.set_index('node_name', inplace=True)
        #
        # Get beam end loads
        fer_global = fer_global.groupby(['node_name'],
                                        as_index=False)[force_header].sum()
        
        fer_global.set_index('node_name', inplace=True)
        # get node 1 and 2 FER forces
        Fb_fer = np.concatenate((fer_global.loc[nodes[0]],
                                 fer_global.loc[nodes[1]]), axis=None)
        #
        # Calculate total end node force in local system
        #Ft_local = Tb @ (Fb_global - FU_fer)
        Fb = Fb_local - Tb @  Fb_fer
        # get reactions with correct sign for postprocessing
        R0 = eq.R0(Fb_local=Fb,
                   torsion=bt_load.loc[nodes[0]])
        #
        # Get axial load
        P0 = 0.0
        if eq.Pdelta:
            P0 = Fb_local[0]
        #
        return P0, R0, Fb
    #    
    def _beam_no_load(self, beam_name: int|str,
                      beam, Fb_local, eq, #Tb,
                      geometry, steps: int)->list[list]:
        """
        Returns:
            [system, beam_name, beam_step, Fx, Fy, Fz, Mx, My, Mz]
        """
        # Dummy Bending [V, M, theta, w]
        Fb_blank = [0, 0, 0, 0]
        # Dummy Torsion [T, Phi, Psi, B, Tw] - [T, theta, theta1, theta2, theta3]
        Ft_blank = [0, 0, 0, 0, 0]
        #
        R0 = eq.R0(Fb_local=Fb_local,
                   torsion=Tload(0, 0, 0))
        #
        Lsteps = linstep(d=geometry.d,
                         L=beam.L, steps=steps)
        #
        lbforce = [['local', beam_name,  bstep,
                    *beam.response(x=bstep, R0=[R0.x, R0.t, R0.y, R0.z],
                                   Fx=[Fb_blank, Ft_blank, Fb_blank, Fb_blank])]
                   for bstep in Lsteps]
        return lbforce
    #
    #
    def get_reactions(self, node_force:DBframework.DataFrame)-> (DBframework.DataFrame, DBframework.DataFrame):
        """ get nodal reactions """
        #
        header = ['mesh_name', 'result_name',
                  'load_name', 'load_level',
                  'node_name', 'system']
        #
        nforce = (node_force.groupby(header,
                                     as_index=False)
                   [self._plane.hforce].sum())
        #
        # TODO : maybe there is a better pandas way that the code below..
        supports = self._mesh._boundaries.support()
        nname = list(supports)
        nrsupp = nforce.loc[nforce['node_name'].isin(nname)]
        db = DBframework()
        nrsupp = db.DataFrame(data=nrsupp)
        return nforce, nrsupp
    #
    # -----------------------------------------------------------
    # DataFrame Operations
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
# -----------------------------------------------------------
# DataFrame Operations
# -----------------------------------------------------------
#
#
def update_combination(member, combination,
                       item: list[str],
                       values:list[str]):
    """
    Update node displacements to include combination
    """
    db = DBframework()
    # group basic load by name
    grp = member.groupby('load_name')
    comb_grp = combination.groupby('load_name')
    temp = []
    for key, comb_factors in comb_grp:
        for row in comb_factors.itertuples():
            try:
                comb = grp.get_group(row.basic_load).copy()
            except KeyError:
                continue
            comb.loc[:, values] *= row.factor
            comb['load_level'] = 'combination'
            comb['load_name'] = row.load_name
            comb['load_id'] = row.load_id
            comb['load_title'] = row.load_title
            temp.append(comb)
    #
    try:
        temp = db.concat(temp, ignore_index=True)
        temp = temp.groupby(item, as_index=False)[values].sum()
        member = db.concat([member, temp], ignore_index=True)
    except (UnboundLocalError, ValueError):
        pass
    return member
#
#
def mif2df2(member_load:list[list]) -> DBframework.DataFrame:
    """
    Return:
        Beams' internal forces dataframe
    """
    header = ['load_name', 'mesh_name', 'load_level',
              'result_name',
              'system', 'element_name', 'length',
              'Fx', 'blank1', 'blank2', 'x',  # axial
              'Mx', 'rx', 'Psi', 'B', 'Tw',  # torsion
              'Fy', 'Mz', 'rz', 'y',  # bending in plane
              'Fz', 'My', 'ry', 'z']  # bending out plane
    Fb_result = list2df(data=member_load, header=header)
    # reorder columns
    # Fb_result['result_name'] = self._result_name
    Fb_result.drop(columns=['blank1', 'blank2'])
    # Reorganizing
    header = ['mesh_name', 'result_name',
              'load_name', 'load_level',
              'element_name', 'length',
              'system',
              'Fx', 'Fy', 'Fz',
              'Mx', 'My', 'Mz',
              'x', 'y', 'z',
              'rx', 'ry', 'rz',
              'Psi', 'B', 'Tw']
    return Fb_result[header]
#
#
def list2df(data:list[list], header:list) -> DBframework.DataFrame:
    """
    Returns:
        dataframe
    """
    db = DBframework()
    dftemp = db.DataFrame(data=data, columns=header, index=None)
    return dftemp
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
    Un_local: list
    plane: bool
    Pdelta: bool
    #
    def R0(self, Fb_local: list, torsion: tuple) -> tuple:
        """
        Axial   [FP, blank, blank, Fu]
        Bending [V, M, theta, w]
        Torsion [T, Phi, Psi, B, Tw]
        """
        # select node disp end 0
        Un0 = self.Un_local[:self.plane.ndof]
        bload = Fb_local[:self.plane.ndof]
        #
        # Axial
        axial = bload[0]
        if np.isclose(a=axial, b=0.0, atol=0.01):
            axial = 0
        #
        R0x = [1 * axial,            # Fx
               0, 0,                 # blank, blank
               -1 * Un0[0]]          # dx
        #
        if self.plane.plane2D:
            R0t, R0y, R0z = self.D2(Un0, bload)
        else:
            R0t, R0y, R0z = self.D3(Un0, bload, torsion)
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
    def D3(self, Un0:list, Fbeam:list,
           torsion:tuple):
        """
        3D plane
        Un0: [dx,dy,dz,rx,ry,rz]
        Fbeam: [Fx,Fy,Fz,mx,my,mz]
        torsion: [T, Phi, Psi, B, Tw]

        Returns:
            Reaction end 0 [torsion, Ry, Rz]

        """
        # In plane
        R0y = [Fbeam[1],   # Vy
               Fbeam[5],   # Mz
               Un0[5],     # rz
               Un0[1]]     # dy
        # Out plane
        R0z = [Fbeam[2],  # Vz
               Fbeam[4],  # My
               Un0[4],    # ry
               Un0[2]]    # dz
        #
        # Torsion
        T0 = -1 * Fbeam[3]
        rx = 1 * Un0[3]
        #
        # Thin walled sections (Ibeam, Channel & Z)
        thetas = [1 * torsion.Psi,
                  1 * torsion.B,
                  1 * torsion.Tw]
        #
        R0t = [T0, rx, *thetas]
        #
        return R0t, R0y, R0z
    #
    def D2(self, Un0:list, Fbeam:list):
        """
        2D plane
        Un0 : [dx, dy, rz]
        Fbeam : [Fx,Fy,mz]

        Returns:
            Reaction end 0 [torsion, Ry, Rz]
        """
        # In plane
        R0y = [Fbeam[1],   # Vy
               Fbeam[2],   # Mz
               Un0[2],     # rz
               Un0[1]]     # dy
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

