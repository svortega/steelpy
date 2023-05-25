# 
# Copyright (c) 2009-2023 fem2ufo
#
from __future__ import annotations
# Python stdlib imports
#import time
#
#
# package imports
from steelpy.process.math.operations import zeros, to_matrix, linspace, mtxmul
from steelpy.process.dataframe.main import DBframework
#
from steelpy.process.math.operations import trns_3Dv, mtxmul, linspace
from steelpy.formulas.main import BeamBasic
#from steelpy.process.dataframe.main import DBframework
#
import numpy as np
#
# -----------------------------------------------------------
#
# -----------------------------------------------------------
#
#
def get_reactions(boundaries, nforce):
    """ get nodal reactions """
    #
    #
    supports = boundaries.supports()
    nname = list(supports)
    #
    nrsupp = nforce.loc[nforce['node_name'].isin(nname)]
    #nrsupp2 =  nreacs.loc[nreacs['node_name'].isin(nname),
    #                      ['Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']].sum()
    # 
    nreaction = (nrsupp.groupby(['load_name', 'load_number', 'load_type',
                                  'system', 'load_title', 'node_name'])
                 [['Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']].sum())
    # reset groupby to create new df
    nreaction = nreaction.reset_index(names=['load_name', 'load_number', 'load_type',
                                             'system', 'load_title', 'node_name'])
    #
    return nreaction
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
    header = ['load_name','load_number', 'load_title', 'load_type', 'system',
              'node_name', 'x', 'y', 'z', 'rx', 'ry', 'rz']    
    dfload = pd.DataFrame(data=dftemp, columns=header, index=None)
    return dfload
#
# -----------------------------------------------------------
# Combination process
# -----------------------------------------------------------
#
#
def updape_ndf(dfnode, dfcomb, 
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
            #check = dftemp.groupby(['node_name', 'c']).sum().reset_index()
            comb = dftemp.groupby(['load_name', 'load_number','load_type',
                                   'load_title', 'system','node_name'],
                                    as_index=False)[values].sum()
            #test
            dfnode = db.concat([dfnode, comb], ignore_index=True)
        #
        return dfnode #, memb_comb
#
#
def updape_memberdf(dfmemb, dfcomb, 
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
            comb = dftemp.groupby(['load_name', 'load_number','load_type',
                                   'load_title', 'system','element_name' ,'node_name'],
                                    as_index=False)[values].sum()
            #test
            dfmemb = db.concat([dfmemb, comb], ignore_index=True)
        #
        return dfmemb
#
#
def updape_memberdf2(dfmemb, dfcomb, 
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
            comb = dftemp.groupby(['load_name', 'load_number','load_type',
                                   'load_title', 'system','element_name' ,'node_end'],
                                    as_index=False)[values].sum()
            #test
            dfmemb = db.concat([dfmemb, comb], ignore_index=True)
        #
        return dfmemb
#
#
#
# -----------------------------------------------------------
# 
# -----------------------------------------------------------
#
def beam_end_force(elements, basic_load, df_ndisp, df_nload):
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
                                'load_type', 'system'])
    nlgrp = df_nload.groupby(['load_name', 'load_number', 'load_title',
                              'load_type', 'system'])
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
                #if len(node1.index) > 1:
                    gnload = np.concatenate((node1.sum(), node2.sum()), axis=None)
                except IndexError:
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
    header: list[str] = ['load_name', 'load_number', 'load_title', 'load_type', 'system', 
                         'element_name', 'node_name', 'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']
    df_nforce = db.DataFrame(data=ntest, columns=header, index=None)
    return df_nforce
#
#
def beam_force(elements, basic_load,
               df_ndisp, df_nforce,
               steps:int = 10):
    """get beam forces for basic loads"""
    #
    print("** Calculating Member Forces")
    bload_func = basic_load.process(elements=elements, steps=steps)
    #    
    ndgrp = df_ndisp.groupby(['load_name', 'load_number', 'load_type', 'load_title'])
    nfgrp = df_nforce.groupby(['load_name', 'load_number', 'load_type', 'load_title'])
    #
    Fblank = [0, 0, 0, 0]
    #Fblank2 = [0, 0, 0, 0]
    dummyf = np.array([0]*6)
    #
    member_load: list = []
    #ndbasic = ndgrp.get_group('basic')
    for key, noded in ndgrp:
        #
        ndisp = noded[['node_name', 'x', 'y', 'z',
                       'rx', 'ry', 'rz']]
        ndisp.set_index('node_name', inplace=True)        
        #
        # get elements
        nfval = nfgrp.get_group(key)
        enfgrp = nfval.groupby('element_name')
        #
        # check if basic load
        try:
            mbload = bload_func[key[0]]
        except KeyError:
            mbload = {}
        #
        for mname, nodef in enfgrp:
            member = elements[mname]
            #beam = member.beam()           
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
            nitem = []
            for node in member.connectivity:
                try:
                    nitem.append(ndisp.loc[node])
                except KeyError:
                    nitem.append(dummyf)
            #
            # ---------------------------------------------
            # get beam end-node displacement in global system
            #
            gndisp = np.concatenate(nitem, axis=None)
            #
            # convert global end-node disp in beam's local system
            #lndisp = beam.transformation(vector=gndisp)
            lndisp = trns_3Dv(gndisp, member.T)
            # displacement end 0
            lndisp0 = lndisp[:6]
            #lndisp1 = lndisp[6:]
            #
            # --------------------------------------------
            # convert beam end-node disp to force [F = Kd] in global system
            #gnforce2 = mtxmul(a=member.K, b=gndisp)
            # FIXME : select node by number
            gnforce = np.concatenate(nodef[['Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']].values)
            # covert global nodal force in beam's local system
            lnforce = trns_3Dv(gnforce, member.T)
            # force end 0
            lnforce0 = lnforce[:6]
            #
            # set beam to general response expresions --> R0
            # [V, M, theta, w]
            #TODO: confirm change reactions sign
            R0x = [1 * lnforce0[0], 0, 0, 1 * lndisp0[0]]
            R0y = [1 * lnforce0[1], 1 * lnforce0[5], -1 * lndisp0[5], 1 * lndisp0[1]]
            R0z = [1 * lnforce0[2], 1 * lnforce0[4], -1 * lndisp0[4], 1 * lndisp0[2]]
            R0t = [1 * lnforce0[3], 0, 0, 1 * lndisp0[3]]
            #
            #Fblank2[0] = -1 * lnforce0[0]  # axial load
            #Fblank2[1] = -1 * lnforce0[3]  # torsion
            #Fblank2[2] = 1 * lndisp0[3]   # theta
            #Fblank2[3] = 1 * lndisp0[0]   # displacement
            #
            # ---------------------------------------------
            #
            try:
                # [load_name, member_load_title, load_type, load_system, 
                # beam_number, x, Fx, Fy, Fz]                
                mnload = mbload[mname] 
                #
                R0y = [-1 * lnforce0[1], 1 * lnforce0[5], -1 * lndisp0[5], 1 * lndisp0[1]]
                R0z = [-1 * lnforce0[2], 1 * lnforce0[4], -1 * lndisp0[4], 1 * lndisp0[2]]
                # TODO : axial function needed here
                #Fblank2[0] = -1 * lnforce0[0]
                #Fblank2[1] =  1 * lnforce0[3]
                #
                # [load_title, load_system, beam_number,  x, Fx, Fy, Fz, Mx, My, Mz]
                lbforce = [[bstep[1], *bstep[3:6],
                            *beam.response(x=bstep[5], R0=[R0x, R0t, R0y, R0z],
                                           Fx=[*bstep[6:]])]
                           for bstep in mnload]
            
            except (KeyError, AttributeError):
                # [beam_number, load_title, x, Fx, Fy, Fz, Mx, My, Mz]
                Lsteps = linspace(start=0, stop=member.L, num=steps+1, endpoint=True)
                lbforce = [[None,'local', mname,  xstep,
                            *beam.response(x=xstep, R0=[R0x, R0t, R0y, R0z],
                                           Fx=[Fblank, Fblank, Fblank, Fblank])]
                           for xstep in Lsteps]
            #
            # Axial   [FP, blank, blank, Fu]
            # torsion [T, B, Psi, Phi, Tw]
            # Bending [V, M, theta, w]
            #
            # ---------------------------------------------
            #
            # member local system
            #
            member_load.extend([[*key,  *lbf[:4],
                                 *lbf[4], # axial
                                 *lbf[5], # torsion
                                 *lbf[6], # bending in plane
                                 *lbf[7]] # bending out plane
                                for lbf in lbforce])
    #
    # --------------------------------------------
    # Seting member df
    # --------------------------------------------
    #
    header = ['load_name', 'load_number', 'load_type', 'load_title',
              'element_load_name', 'system',
              'element_name', 'node_end',
              'F_Vx', 'blank1', 'blank2', 'F_wx',  # axial
              'F_Mx', 'F_B', 'F_psi', 'F_phix',    # torsion
              'F_Vy', 'F_Mz', 'F_thetaz', 'F_wy',  # bending in plane
              'F_Vz', 'F_My', 'F_thetay', 'F_wz']  # bending out plane
    #
    db = DBframework()
    df_membf = db.DataFrame(data=member_load, columns=header, index=None)
    # reorder columns
    df_membf = df_membf[['load_name', 'load_number', 'load_type', 'load_title', 'system',
                         'element_name', 'element_load_name', 'node_end',
                         'F_Vx', 'F_Vy', 'F_Vz', 'F_Mx', 'F_My', 'F_Mz',
                         'F_wx', 'F_wy', 'F_wz', 'F_phix', 'F_thetay', 'F_thetaz']]
    return df_membf
#
#
# -----------------------------------------------------------
# 
# -----------------------------------------------------------
#
