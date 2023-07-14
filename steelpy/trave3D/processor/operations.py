# 
# Copyright (c) 2009-2023 fem2ufo
#
from __future__ import annotations
# Python stdlib imports
#from collections.abc import Mapping
from dataclasses import dataclass
from typing import NamedTuple
import itertools as it
#from itertools import chain
#import pickle
#import time
#
# package imports
from steelpy.formulas.main import BeamBasic
from steelpy.process.math.operations import (zeros, to_matrix, mtxmul,
                                             trnsload, linspace)
from steelpy.process.dataframe.main import DBframework
import numpy as np
#
#
#
# -----------------------------------------------------------
#
# -----------------------------------------------------------
#
def max_bandwidth(elements,  jbc):
    """
    calculate max bandwidth
    ------------------------  
    npi : connectivity end 1
    npj : connectivity end 2
    jbc : nodes freedom
    nel: number of elements
    if we
    npi ,npj, jbc, nel
    """
    ibndm4 = [0]
    for key, element in elements.items():
        idof, jdof = element.DoF
        bc1 = jbc[idof]
        bc2 = jbc[jdof]
        ieqn = bc1 + bc2
        try:
            ibndm4.append(max([abs(ieqn1 - ieqn2)
                               for x, ieqn1 in enumerate(ieqn) if ieqn1 > 0
                               for ieqn2 in ieqn[x+1:] if ieqn2 > 0]))
        except ValueError:
            continue
    #
    return max(ibndm4) + 1
#
def bd_condition(nodes, boundaries):
    """
    set boundary conditions
    bcs: set default = 0 free (1 fix)
    """
    nnp = len(nodes)
    jbc = zeros(nnp, 6, code='I')
    supports = boundaries.supports()
    #for node_name, bd in boundaries.node.items():
    for node_name, bd in supports.items():
        ind = nodes[bd.node].index
        jbc[ind] = list(bd[:6])
    return jbc
#
def spclbc(elements, nodes, free_nodes, jbc):
    """
    Impose condition for  special GLOBAL shapes
    bcs: set default = 0 free (1 fix)
    """
    # free_nodes = elements.get_free_nodes()
    #
    #for _element in elements.values():
    #for element in elements.iter_elements:
    for key, element in elements.items():
        conn = element.connectivity
        pos_node = set(conn) - set(free_nodes)
        for node_name in pos_node:
            ind = nodes[node_name].index
            #ind
            # jbc[ind][3:] = [1, 1, 1]
            # jbc[ind][3] = 1
            # jbc[ind][4] = 1
        # if _element.type == "truss":
        #    
        #    # end 1
        #    ind = self._nodes._labels.index(conn[0])
        #    jbc[ind][3:] = [0, 0, 0]
        #    # end 2
        #    ind = self._nodes._labels.index(conn[1])
        #    jbc[ind][3:] = [0, 0, 0]
    return jbc
#
def shape_cond(elements, nodes, boundaries, free_nodes):
    """
    jcs: modify default = 0 free (1 fix)
    """
    jbc = bd_condition(nodes, boundaries)
    #jbcX = nodes.neq()
    # TODO : check this module
    #jbc = spclbc(elements, nodes, free_nodes, jbc)
    #
    # Number the equations  in jbc from 1 up to the order.
    # Start assigning equation numbers for zero dof's
    # from 1 up;  only zero given a number.        
    #
    jbc = list(it.chain.from_iterable(jbc))
    counter = it.count(start=1)
    jbc = [next(counter) if _item == 0 else 0
           for _item in jbc]
    neq = max(jbc)
    jbc = to_matrix(jbc, 6)
    return jbc, neq
#
def get_bandwidth(elements, nodes, boundaries):
    """ """
    #jbcc, neq = shape_cond(elements=elements,
    #                       nodes=nodes,
    #                       boundaries=boundaries,
    #                       free_nodes=free_nodes)
    #
    jbc, neq = nodes.neq(supports=boundaries._nodes)
    #
    iband = max_bandwidth(elements=elements, jbc=jbc)
    return jbc, neq, iband
#
#
#
# -----------------------------------------------------------
# 
# -----------------------------------------------------------
#
def beam_end_forceX(elements, basic_load, df_ndisp, df_nload):
    """Convert node displacement to beam end forces [F = Kd]"""
    # get jbc
    #file = open("stfmx.f2u", "rb")
    #df_jbc = pickle.load( file )
    #file.close()
    #
    # reset jbc to zeros
    #df_jbc.iloc[:] = 0
    #
    #
    dispgrp = df_ndisp.groupby(['load_name', 'load_number', 'load_title',
                                'load_type', 'load_system'])
    nlgrp = df_nload.groupby(['load_name', 'load_number', 'load_title',
                              'load_type', 'load_system'])
    #
    dummyf = np.array([0]*6)
    ntest = []
    for key, item in nlgrp:
        nload = item[['node_name', 'Fx', 'Fy', 'Fz',
                      'Mx', 'My', 'Mz']] #.values
        nload.set_index('node_name', inplace=True)
        #
        df1 = dispgrp.get_group(key)
        df1 = df1[['node_name','x', 'y', 'z', 'rx', 'ry', 'rz']]
        df1.set_index('node_name', inplace=True)
        #
        #df2 = df_jbc.copy()
        #df2.loc[df2.index.isin(df1.index)] = df1[['x', 'y', 'z', 'rx', 'ry', 'rz']]
        #
        # ----------------
        # get basic load case
        lcase = basic_load[key[0]]
        # TODO: can be a plate element
        beamload = lcase.beam()
        #
        for mname, element in elements.items():
            nodes = element.connectivity
            #
            bload = beamload[mname]
            # ---------------------------------------------
            # get beam end-node displacement in global system
            ndisp = []
            for node in nodes:
                try:
                    ndisp.append(df1.loc[node])
                except KeyError:
                    ndisp.append(dummyf)
            gndisp = np.concatenate(ndisp, axis=None)
            #gndisp = np.concatenate((df1.loc[nodes[0]],
            #                         df1.loc[nodes[1]]), axis=None)
            #
            # ---------------------------------------------
            # convert beam end-node disp to force [F = Kd] in global system
            gnforce = mtxmul(element.K, gndisp)
            # ---------------------------------------------
            # check if load on member
            try:
                bltotal = 1 / (len(bload.line) + len(bload.point))
                # ---------------------------------------------
                # node force global system
                node1 = nload.loc[nodes[0]]
                node2 = nload.loc[nodes[1]]
                try:
                    node1.shape[1]
                    node1 = node1.sum()
                except IndexError:
                    pass
                
                try:
                    node2.shape[1]
                    node2 = node2.sum()
                except IndexError:
                    pass                
                
                #if len(node1.index) > 1:
                #gnload = np.concatenate((node1.sum(), node2.sum()), axis=None)
                #except IndexError:
                gnload = np.concatenate((node1, node2), axis=None)
                #
                nreac0 = gnforce[:6] + gnload[:6]
                nreac1 = gnforce[6:] + gnload[6:]
                #
                ntest.append([*key, mname, nodes[0], *nreac0])
                ntest.append([*key, mname, nodes[1], *nreac1])
            except ZeroDivisionError:
                ntest.append([*key, mname, nodes[0], *gnforce[:6]])
                ntest.append([*key, mname, nodes[1], *gnforce[6:]])
    #
    # get df 
    db = DBframework()
    header: list[str] = ['load_name', 'load_number', 'load_title',
                         'load_type', 'load_system', 
                         'element_name', 'node_name',
                         'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']
    df_nforce = db.DataFrame(data=ntest, columns=header, index=None)
    return df_nforce
#
#
def beam_end_force(elements, df_ndisp,
                   m2D:bool = False):
    """Convert node displacement to beam end forces [F = Kd]"""
    # setup 
    ndof:int = 6
    headdisp = ['node_name','x', 'y', 'z', 'rx', 'ry', 'rz']
    headforce = ['node_name','Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']
    if m2D:
        ndof:int = 3
        headdisp = ['node_name','x', 'y', 'rz']
        headforce = ['node_name','Fx', 'Fy', 'Mz']
    #
    dummyf = np.array([0]*ndof)
    dispgrp = df_ndisp.groupby(['load_name', 'load_number', 'load_title',
                                'load_type', 'load_system'])
    #
    ntest = []
    for key, item in dispgrp:
        df1 = item[headdisp]
        df1.set_index('node_name', inplace=True)
        # start element big loop
        for mname, element in elements.items():
            nodes = element.connectivity
            # ---------------------------------------------
            # get beam end-node displacement in global system
            #ndisp = []
            #for node in nodes:
            #    try:
            #        ndisp.append(df1.loc[node])
            #    except KeyError:
            #        ndisp.append(dummyf)
            #gndisp = np.concatenate(ndisp, axis=None)
            gndisp = np.concatenate((df1.loc[nodes[0]],
                                     df1.loc[nodes[1]]), axis=None)
            #
            # ---------------------------------------------
            # convert beam end-node disp to force [F = Kd] in global system
            gnforce = mtxmul(element.K(m2D=m2D), gndisp)
            # ---------------------------------------------
            #
            ntest.append([*key, mname, nodes[0], *gnforce[:ndof]])
            ntest.append([*key, mname, nodes[1], *gnforce[ndof:]])
    #
    # get df 
    db = DBframework()
    header: list[str] = ['load_name', 'load_number', 'load_title',
                         'load_type', 'load_system', 
                         'element_name', *headforce]
    df_nforce = db.DataFrame(data=ntest, columns=header, index=None)
    return df_nforce
#
#
# -----------------------------------------------------------
# 
# -----------------------------------------------------------
#
#
def beam_int_force(elements, basic_load,
                   df_ndisp, df_nforce,
                   df_nload, 
                   steps:int = 10,
                   m2D:bool = False):
    """get beam forces for basic loads"""
    print("** Calculating Member Forces")
    #
    ndof:int = 6
    headforce = ['Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']
    headdisp = ['node_name', 'x', 'y', 'z', 'rx', 'ry', 'rz']
    headload = ['axial', 'torsion', 'VM_inplane', 'VM_outplane']
    #trns = trnsload
    if m2D:
        ndof:int = 3
        #df_jbc = df_jbc[['x', 'y', 'rz']]
        headforce = ['Fx', 'Fy', 'Mz']
        headdisp = ['node_name','x', 'y', 'rz']
        #trns = trnsload_2Dv
    #
    bload_func = basic_load.process(elements=elements, steps=steps)
    #
    # -------------
    #comb0 = df_nforce.groupby(['load_name', 'load_number','load_type',
    #                           'load_title', 'load_system','node_name'],
    #                            as_index=False)[headforce].sum()
    #comb1 = comb0.set_index('node_name')
    # -------------
    #
    ndgrp = df_ndisp.groupby(['load_name', 'load_type', 'load_number', 'load_title'])
    nfgrp = df_nforce.groupby(['load_name',  'load_type', 'load_number', 'load_title'])
    nlgrp = df_nload.groupby(['load_name',  'load_type', 'load_number', 'load_title'])
    blgrp = bload_func.groupby(['load_name',  'load_type'])
    #
    # Dummy Bending [V, M, theta, w]
    Fblank = [0, 0, 0, 0]
    # Dummy Torsion [T, B, Psi, Phi, Tw]
    Fblank_t = [0, 0, 0, 0, 0]
    dummyf = np.array([0]*ndof)
    #
    member_load: list = []
    #ndbasic = ndgrp.get_group('basic')
    for key, noded in ndgrp:
        #
        ndisp = noded[headdisp]
        ndisp.set_index('node_name', inplace=True)        
        #
        # get elements
        nfval = nfgrp.get_group(key)
        enfgrp = nfval.groupby('element_name')
        #
        # check if basic load
        try:
            mbload = blgrp.get_group(key[:2])
            mbload = mbload.groupby(['element_name'])
            #mbload = bload_func[key[0]]
            #
            nlval = nlgrp.get_group(key)
            nlval = nlval.groupby(['element_name'])
        except KeyError:
            mbload = {}
        #
        for mname, nodef in enfgrp:
            member = elements[mname]
            Tlg = member.T(m2D=m2D)           
            #
            material = member.material
            section = member.section.properties()
            beam = BeamBasic(L=member.L, area=section.area, 
                             Iy=section.Iy, Iz=section.Iz,
                             J=section.J, Cw=section.Cw, 
                             E=material.E, G=material.G)         
            #
            #in1, in2 = node_id.index(nodes[0]),  node_id.index(nodes[1]) 
            #n1, n2 = member.connectivity
            #nlitem = []
            #for node in member.connectivity:
            #    try:
            #        nlitem.append(nlval.loc[node].values)
            #    except KeyError:
            #        nlitem.append(dummyf)
            #
            #nlitem = np.concatenate(nlitem, axis=None)
            #
            nodes = member.connectivity
            nditem = np.concatenate((ndisp.loc[nodes[0]],
                                     ndisp.loc[nodes[1]]), axis=None)
            #
            #nfitem = np.concatenate((comb1.loc[nodes[0]][headforce].values,
            #                         comb1.loc[nodes[1]][headforce].values), axis=None)             
            # ---------------------------------------------
            # get beam end-node displacement in global system
            #
            gndisp = np.concatenate(nditem, axis=None)
            #
            # convert global end-node disp in beam's local system
            #lndisp = trns(gndisp, Tlg)
            #lndisp2 = gndisp @ Tlg
            lndisp = Tlg @ gndisp
            # displacement end 0
            #lndisp0 = lndisp[:ndof]
            #lndisp1 = lndisp[6:]
            #
            # --------------------------------------------
            # convert beam end-node disp to force [F = Kd] in global system
            #
            #gnforce = mtxmul(member.K(m2D=m2D), gndisp)
            lnforce = mtxmul(member.k(m2D=m2D), lndisp)
            #gnforce = np.concatenate(nfitem, axis=None)
            #gnforce += nlitem            
            #
            # FIXME : select node by number
            #gnforce = np.concatenate(nodef[headforce].values)
            # covert global nodal force in beam's local system
            #lnforce = trns(gnforce, Tlg)
            #lnforce2 = Tlg @ gnforce.T
            #lnforce3 = Tlg @ gnforce
            # force end 0
            #lnforce0 = lnforce[:ndof]
            #
            # --------------------------------------------
            # set beam to general response expresions --> R0
            # [V, M, theta, w]
            #
            #TODO: confirm change reactions sign
            eq = NodeGenRespEq(lnforce, lndisp, ndof, m2D)
            #
            #
            # ---------------------------------------------
            #
            try:
                # Beam load (udl/point)
                # [load_name, member_load_title, load_type, load_system, 
                # beam_number, x, Fx, Fy, Fz]                
                #mnload = mbload[mname]
                mnload = mbload.get_group(mname)
                mnload = mnload.groupby(['node_end'], as_index=False)[headload].sum()
                #
                #
                #nlval = nlgrp.get_group(key)
                #nlval = nlval.groupby(['element_name'])
                nodeloads = nlval.get_group(mname)
                nodeloads = nodeloads.groupby(['node_name'], as_index=False)[headforce].sum()
                nodeloads = nodeloads.set_index('node_name')                 
                #
                blitem = np.concatenate((nodeloads.loc[nodes[0]],
                                         nodeloads.loc[nodes[1]]), axis=None)                
                #
                #blitem2 = trns(blitem, member.T(m2D=m2D))
                #Tlg = member.T(m2D=m2D)
                #blitem = Tlg @ blitem.T
                blitem = Tlg @ blitem
                #
                R0 = eq.R0(bload=blitem)
                #R0 = eq.R0()
                #
                #
                #lbforce = []
                #for bstep in  mnload.itertuples():
                #    #bstep
                #    lbforce.append(['local', mname, bstep.node_end, 
                #                    *beam.response(x=bstep.node_end, R0=[R0.x, R0.t, R0.y, R0.z],
                #                                   Fx=[*bstep[2:]])])
                #
                # [load_title, load_system, beam_number,  x, Fx, Fy, Fz, Mx, My, Mz]
                #lbforce = [[bstep[1], *bstep[3:6],
                #            *beam.response(x=bstep[5], R0=[R0.x, R0.t, R0.y, R0.z],
                #                           Fx=[*bstep[6:]])]
                #           for bstep in mnload]
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
                Lsteps = linspace(start=0, stop=member.L, num=steps+1, endpoint=True)
                lbforce = [['local', mname,  xstep,
                            *beam.response(x=xstep, R0=[R0.x, R0.t, R0.y, R0.z],
                                           Fx=[Fblank, Fblank_t, Fblank, Fblank])]
                           for xstep in Lsteps]
            #
            # ---------------------------------------------
            # Member total force in local system
            #
            # Axial   [FP, blank, blank, Fu]
            # torsion [T, B, Psi, Phi, Tw]
            # Bending [V, M, theta, w]
            #
            # ---------------------------------------------
            #
            member_load.extend([[*key,    # load_name
                                 *lbf[:3],# 'load_number', 'load_type'
                                 *lbf[3], # axial
                                 *lbf[4], # torsion
                                 *lbf[5], # bending in plane
                                 *lbf[6]] # bending out plane
                                for lbf in lbforce])
    #
    # --------------------------------------------
    # Seting member df
    # --------------------------------------------
    #
    header = ['load_name', 'load_type', 'load_number','load_title',
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
    df_membf = df_membf[['load_name', 'load_number',
                         'load_type', 'load_title', 
                         'load_system',
                         'element_name', 'node_end',
                         'F_Vx', 'F_Vy', 'F_Vz',
                         'F_Mx', 'F_My', 'F_Mz',
                         'F_wx', 'F_wy', 'F_wz',
                         'F_phix', 'F_thetay', 'F_thetaz']]
    #
    #
    return df_membf
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
#
# -----------------------------------------------------------
# Combination process
# -----------------------------------------------------------
#
#
def comb_update(df_comb, df_nload, df_ndisp, df_nforce, df_membf,
                m2D: bool):
    """Update df to include load combinations"""
    # ---------------------------------------
    # Update df to include load combinations
    # ---------------------------------------
    #
    disphead = ['x', 'y', 'z', 'rx', 'ry', 'rz']
    forcehead = ['Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']
    if m2D:
        forcehead = ['Fx', 'Fy', 'Mz']
        disphead = ['x', 'y', 'rz']
    #
    df_ndisp = update_ndf(dfnode=df_ndisp, dfcomb=df_comb,
                          values=disphead)
    #
    df_nload = update_ndf(dfnode=df_nload, dfcomb=df_comb,
                          values=forcehead)
    #
    df_nforce = update_memberdf(dfmemb=df_nforce, dfcomb=df_comb,
                                values=forcehead)
    #
    df_membf = update_memberdf2(dfmemb=df_membf, dfcomb=df_comb,
                                values=['F_Vx', 'F_Vy', 'F_Vz',
                                        'F_Mx', 'F_My', 'F_Mz',
                                        'F_wx', 'F_wy', 'F_wz',
                                        'F_phix', 'F_thetay', 'F_thetaz'])
    #
    return df_nload, df_ndisp, df_nforce, df_membf
#    
#
# -----------------------------------------------------------
#
def update_ndf(dfnode, dfcomb, 
               values:list[str]): #['x', 'y', 'z', 'rx', 'ry', 'rz']
    """
    Update node displacements to include lcomb
    """
    db = DBframework()
    # group basic load by name
    ndgrp = dfnode.groupby('load_name')
    # get combinations with basic loads 
    #
    combgrp = dfcomb.groupby('load_name')
    for key, combfactors in combgrp:
        for row in combfactors.itertuples():
            comb = ndgrp.get_group(row.basic_load).copy()
            comb.loc[:, values] *= row.factor
            comb['load_type'] = 'combination'
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
        #check = dftemp.groupby(['node_name', 'c']).sum().reset_index()
        dftemp = dftemp.groupby(['load_name', 'load_number','load_type',
                                 'load_title', 'load_system','node_name'],
                                  as_index=False)[values].sum()
        #test
        dfnode = db.concat([dfnode, dftemp], ignore_index=True)
    except UnboundLocalError:
        pass
    #
    return dfnode #, memb_comb
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
            comb['load_type'] = 'combination'
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
        dftemp = dftemp.groupby(['load_name', 'load_number','load_type',
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
def update_memberdf2(dfmemb, dfcomb, 
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
            comb['load_type'] = 'combination'
            comb['load_name'] = row.load_name
            comb['load_number'] = row.load_number
            comb['load_title'] = row.load_title
            #
            try:
                dftemp = db.concat([dftemp, comb], ignore_index=True)
            except UnboundLocalError:
                dftemp = comb
        #
        #comb = dftemp.groupby(['load_name', 'load_number','load_type',
        #                       'load_title', 'load_system',
        #                       'element_name' ,'node_end'],
        #                        as_index=False)[values].sum()
        #test
        #dfmemb2 = db.concat([dfmemb2, comb], ignore_index=True)
    #
    try:
        dftemp = dftemp.groupby(['load_name', 'load_number','load_type',
                                 'load_title', 'load_system',
                                 'element_name' ,'node_end'],
                                  as_index=False)[values].sum()    
        #
        dfmemb = db.concat([dfmemb, dftemp], ignore_index=True)
    except UnboundLocalError:
        pass
    #
    return dfmemb
#
#
# -----------------------------------------------------------
#
# -----------------------------------------------------------
#
#
def get_reactions(boundaries, nforce, m2D: bool):
    """ get nodal reactions """
    forcehead = ['Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']
    if m2D:
        forcehead = ['Fx', 'Fy', 'Mz']
    #
    supports = boundaries.supports()
    nname = list(supports)
    #
    nrsupp = nforce.loc[nforce['node_name'].isin(nname)]
    #nrsupp2 =  nreacs.loc[nreacs['node_name'].isin(nname),
    #                      ['Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']].sum()
    # 
    nreaction = (nrsupp.groupby(['load_name', 'load_number', 'load_type',
                                  'load_system', 'load_title', 'node_name'])
                 [forcehead].sum())
    # reset groupby to create new df
    nreaction = nreaction.reset_index(names=['load_name', 'load_number', 'load_type',
                                             'load_system', 'load_title', 'node_name'])
    #
    return nreaction
#
#
# -----------------------------------------------------------
#
# -----------------------------------------------------------
#
#
def get_deflection(nodes, basic_load, load_combination,
                   basic_res, comb_res):
    """Expand displacement to full vector"""
    print("** assigning node displacement")
    # basic load
    dftemp = []
    #bload = LoadResult()
    # TODO: node name instead number
    node_name = list(nodes.keys())
    for name, basic in basic_res.items():
        blitem = basic_load[name]
        # convert list to node
        ndisp = to_matrix(basic, 6)
        # expand list of node to include number, title & system
        filldata = [[name, blitem.number, blitem.title, 'basic', 'global',
                     nname, *ndisp[x]]
                    for x, nname in enumerate(node_name)]
        dftemp.extend(filldata)
    #dfload = pd.DataFrame(data=dftemp, columns=header, index=None)
    #
    # load combination
    for name, comb in comb_res.items():
        lcitem = load_combination[name]
        ndisp = to_matrix(comb, 6)
        #
        filldata = [[name, lcitem.number, lcitem.title, 'combination', 'global',
                     nname, *ndisp[x]]
                    for x, nname in enumerate(node_name)]
        dftemp.extend(filldata)
    #
    header = ['load_name','load_number', 'load_title',
              'load_type', 'load_system',
              'node_name', 'x', 'y', 'z', 'rx', 'ry', 'rz']    
    dfload = pd.DataFrame(data=dftemp, columns=header, index=None)
    return dfload
#