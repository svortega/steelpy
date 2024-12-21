#
# Copyright (c) 2009 steelpy
#
from __future__ import annotations
# Python stdlib imports
#from array import array
#from copy import copy
#from math import fsum
#import pickle
from dataclasses import dataclass
#from typing import NamedTuple
#from itertools import chain
import time
#
# package imports
from steelpy.utils.math.operations import to_matrix
from steelpy.utils.dataframe.main import DBframework
#from steelpy.utils.math.operations import remove_column_row
#
from steelpy.trave.postprocess.main import PostProcess
#
import numpy as np
#from scipy.linalg import cholesky_banded, cho_solve_banded
from scipy.sparse.linalg import spsolve


#
#
# ---------------------------------------------
# solver using numpy
# ---------------------------------------------
#
def solver_np(a, b):
    """
    a : Global stiffness matrix
    b : Global force vector
    """
    #nloads = df2.stack() #.values
    ndisp = np.linalg.solve(a, b)
    #if not np.allclose(np.dot(stf, ndisp), nload): #, rtol=1e-05, atol=1e-08
    #    raise RuntimeError('Solution fail')
    #
    return ndisp


#
def solver_sparse(a, b):
    """
    """
    ndisp = spsolve(a.tocsr(), b)
    return ndisp


#
#
# ---------------------------------------------
#
@dataclass
class StaticSolver:
    """ Linear static solver class"""
    __slots__ = ['_mesh', '_method', '_postprocess',
                 '_log', '_result_name', 'db_file',
                 'second_order', 'ineleastic']

    def __init__(self, mesh,
                 result_name: int | str,
                 db_file: str, log: bool,
                 second_order: bool = False,
                 ineleastic: bool = False,
                 ) -> None:
        """
        plane : Plane system (3D/2D)
        """
        self._mesh = mesh
        self._result_name = result_name
        self.db_file = db_file
        self._postprocess = PostProcess(mesh=self._mesh,
                                        result_name=result_name,
                                        db_file=db_file)
        self._log = log
        self.second_order = second_order
        self.ineleastic = ineleastic

    #
    def Ke(self, sparse: bool):
        """Stiffness matrix"""
        return self._mesh.Ke(sparse=sparse)

    #
    # ------------------------------------------
    #
    def Dn(self, load_level: str) -> DBframework.DataFrame|None:
        """
        load_level: basic/combination

        Return:
            Nodal global displacement dataframe
        """
        if load_level in ['combination']:
            Dn = self._mesh._load._combination.ND_global(plane=self._mesh._plane)
        else:  # basic
            Dn = self._mesh._load._basic.ND_global(plane=self._mesh._plane)
        #
        if Dn.empty:
            return Dn
        #
        col_disp = self._mesh._plane.hdisp
        jbc = self._mesh.jbc()
        dfjbc = jbc[jbc.any(axis=1)]
        #
        Di = Dn.loc[Dn['node_name'].isin(dfjbc.index)]
        Di = Di.loc[Di[col_disp].any(axis=1)]
        #
        col_grp = ['load_name', 'load_id',
                   'load_level', 'load_title',
                   'load_system', 'mesh_name',
                   'node_name', 'node_index',
                   *col_disp]
        Di = Di[col_grp]
        return Di

    #
    def Fn(self, load_level: str) -> DBframework.DataFrame|None:
        """
        load_level: basic/combination

        Return:
            Nodal global force dataframe
        """
        if load_level in ['combination']:
            Fn = self._mesh._load._combination.NF_global(plane=self._mesh._plane)

        else:  # basic
            Fn = self._mesh._load._basic.NF_global(plane=self._mesh._plane)
        #
        if Fn.empty:
            return Fn
        #
        col_force = self._mesh._plane.hforce
        jbc = self._mesh.jbc()
        dfjbc = jbc[jbc.any(axis=1)]
        #
        Fi = Fn.loc[Fn['node_name'].isin(dfjbc.index)]
        #Fi = Fi.loc[Fi[col_force].any(axis=1)]
        #
        if Fi.empty:
            Fi = Fn
        #
        columns = ['load_name', 'load_id',
                   'load_level', 'load_title',
                   'load_system', 'mesh_name',
                   'node_name', 'node_index',
                   *col_force]
        Fi = Fi[columns]
        return Fi

    #
    # ------------------------------------------
    #
    def _Msplit(self, Km, jbcflat: list):
        """
        Km: The un-partitioned matrix (or vector) to be partitioned.
        
        Return:
            m1 : 
            m2 : 
        """
        jbcc = jbcflat.values
        #
        # List of the indices for degrees of freedom with unknown displacements.
        D1 = [i for i, item in enumerate(jbcc)
              if item != 0]
        # List of the indices for degrees of freedom with known displacements.
        D2 = [i for i, item in enumerate(jbcc)
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
    #
    def _update_Un_comb(self, Un: DBframework.DataFrame) -> DBframework.DataFrame:
        """
        Update node displacements to include lcomb

        Un : node displacement dataframe
        Return:
            Un + load combinations
        """
        disp_head = self._mesh._plane.hdisp
        # group basic load by name
        Un_grp = Un.groupby('load_name')
        #
        # get combinations with basic loads
        load_combination = self._mesh._load.combination()
        lcomb2basic = load_combination.to_basic()
        comb_grp = lcomb2basic.groupby('load_name')
        #
        df_temp = []
        for key, comb_factors in comb_grp:
            for row in comb_factors.itertuples():
                comb_item = Un_grp.get_group(row.basic_load).copy()
                comb_item.loc[:, disp_head] *= row.factor
                comb_item['load_level'] = 'combination'
                comb_item['load_name'] = row.load_name
                comb_item['load_id'] = row.load_id
                comb_item['load_title'] = row.load_title
                comb_item['mesh_name'] = row.mesh_name
                df_temp.append(comb_item)
        #
        db = DBframework()
        columns = ['load_name', 'load_id', 'load_level',
                   'load_title', 'load_system',
                   'mesh_name', 'result_name',
                   'node_name']
        try:
            df_temp = db.concat(df_temp, ignore_index=True)
            df_temp = df_temp.groupby(columns, as_index=False)[disp_head].sum()
            Un = db.concat([Un, df_temp], ignore_index=True)
        except ValueError:
            pass
        #
        return Un

    #
    #
    def _linear(self, load_level: str,
                #col_grp: list[str],
                sparse: bool = False):
        """
        load_level: basic/combination
        columns
        sparse: matrix operation (True/False)
        """
        # ------------------------------
        # Get global nodal load
        #col_grp = ['load_name', 'load_id',
        #           'load_level', 'load_system', 'load_title',
        #           'mesh_name']
        # Get basic data
        # ------------------------------          
        Ke = self._mesh.Ke(sparse=sparse)
        #
        # ------------------------------       
        # and displacement (combination)
        Dn = self.Dn(load_level)
        Fn = self.Fn(load_level)
        # ------------------------------
        # linear solution
        Un = self.PMT(Ke=Ke,
                      Fn=Fn, Dn=Dn,
                      sparse=sparse)
        return Un , Ke, Fn, Dn

    #
    # ------------------------------------------
    #
    def solve_Linear(self,
                     sparse: bool = False,
                     max_iter: int | None = None):
        #load_level: str = 'basic'):
        """
        Linear Static Analysis

        Input:
        Ke  : Global stiffness matrix
        Fn  : FER Node load dataframe
        jbc : Node boundary condition dataframe

        Return:
        Udf : Node displacement global system dataframe
        """
        start_time = time.time()
        #
        # ------------------------------
        # Get global nodal load
        #col_grp = ['load_name', 'load_id',
        #           'load_level', 'load_system',
        #           'load_title', 'mesh_name']
        #     
        # ------------------------------
        # Step 1 : linear solution
        # Un, Fn_grp, Dn_grp =
        Un, Ke, Fn, Dn = self._linear(load_level='basic',
                                      #col_grp=col_grp,
                                      sparse=sparse)
        #
        # ------------------------------
        # Post-processing load combinations
        Un = self._update_Un_comb(Un=Un)
        #
        uptime = time.time() - start_time
        print(f"** {{F}} = [Ke] {{U}} Solution: {uptime:1.4e} sec")
        return Un

    #
    def solve_Pdelta(self,
                     sparse: bool = False,
                     max_iter: int = 30):
        """
        Linear Static Analysis (2nd Order)
        Approximate method: Two cycles iterative method.
        
        Input: 
        Ks  : Global stiffness matrix
        Fn  : FER Node load dataframe
        jbc : Node boundary condition dataframe
        
        Return: 
        Udf : Node displacement global system dataframe
        """
        # ------------------------------
        # start process
        # ------------------------------
        start_time = time.time()
        #
        # Get global nodal load
        col_grp = ['load_name', 'load_id',
                   'load_level', 'load_system', 'load_title',
                   'mesh_name']
        #       
        # ------------------------------
        # Step 1 : linear solution
        # ------------------------------
        # 
        Un, Ke, Fn, Dn = self._linear(load_level='combination',
                                      #col_grp=col_grp,
                                      sparse=sparse)
        #
        # ------------------------------
        # Step 2 : Pdelta solution
        # ------------------------------
        # 
        Un_grp = Un.groupby(col_grp)
        #
        Un_temp = []
        for key, Un_step in Un_grp:
            # ------------------------------
            if Fn.empty:
                Fn_set = Fn
            else:
                Fn_set = Fn.groupby(col_grp).get_group(key)
            #
            if Dn.empty:
                Dn_set = Dn
            else:
                Dn_set = Dn.groupby(col_grp).get_group(key)
            # ------------------------------
            # Assembly matrix including node force
            ke = self._mesh.Kt(Dn=Un_step)
            # ------------------------------
            Un_temp.append(self.PMT(Ke=ke,
                                    Fn=Fn_set,
                                    Dn=Dn_set,
                                    sparse=sparse))
            # ------------------------------
        #
        # Results to dataframe
        db = DBframework()
        Us = db.concat(Un_temp, ignore_index=True)
        # ------------------------------
        # end process
        # ------------------------------
        uptime = time.time() - start_time
        print(f"** {{F}} = [Kt] {{U}} Solution: {uptime:1.4e} sec")
        return Us

    #
    #    
    # ------------------------------------------
    #
    def _get_solver(self, Ke, sparse: bool):
        """
        """
        if sparse:
            solver = solver_sparse
            # Ke = Ke.tolil()
            Kfactor = np.max(Ke.data.max())
        else:
            solver = solver_np
            Kfactor = Ke.max()

        return Kfactor, solver

    #
    def _Us(self, Ke, Fs, solver,
            jbc: list, ndof: float):
        """ """
        # jbc Matrix to vector
        jbcflat = jbc.stack()
        # FIXME: Matrix condensed
        K11, K12 = self._Msplit(Ke, jbcflat)
        # Solve displacements U
        #try:
        Us = iter(solver(K11, Fs))
        #except Warning:
        #    print('-->')
        # reshape vector in matrix form [row, col]
        Us = [next(Us) if ieqnum != 0 else ieqnum
              for ieqnum in jbcflat]
        # bak to matrix form
        Us = to_matrix(Us, ndof)
        return Us

    #
    def _mask_DOF(self, jbc, rename: bool = True):
        """ """
        # remove rows with zeros
        dfjbc = jbc.copy()
        if rename:
            dfjbc.rename(columns=self._mesh._plane.colrename,
                         inplace=True)
        dfjbc = dfjbc[jbc.any(axis=1)]
        dfjbc = dfjbc.replace(float(0.0), np.nan)
        #
        dfbool = dfjbc.copy()
        dfbool = dfbool.notnull()
        dfbool = dfbool.stack(future_stack=True)
        # Copy dataframe
        # dfzeros = dfjbc.copy()
        # dfzeros.iloc[:] = float(0.0)
        #
        dfjbc.iloc[:] = float(0.0)
        # for col in dfjbc.columns:
        #    dfjbc[col].values[:] = float(0.0)
        #
        return dfbool.astype('bool'), dfjbc.astype('float64')

    #
    def _to_df(self, Us):
        """displacement dataframe"""
        db = DBframework()
        header = ['load_name', 'load_id', 'load_level',
                  'load_system', 'load_title', 'mesh_name',
                  'result_name', 'node_name',
                  *self._mesh._plane.hdisp]
        dftemp = db.DataFrame(data=Us, columns=header, index=None)
        return dftemp  #[header]

    #
    # ------------------------------------------
    #
    def _get_FnDn(self, Fn, Dn,
                  colgrp: list[str]|None=None):
        """
        Return:
            group: 
            Fn_group
            Dn_group
        """
        if not colgrp:
            colgrp = ['load_name', 'load_id',
                      'load_level', 'load_system',
                      'load_title', 'mesh_name']
        #
        # Force
        Fn_grp = None
        if not Fn.empty:
            Fn_grp = Fn.groupby(colgrp)
            Fn_keys = list(Fn_grp.groups.keys())
        #
        # Displacement
        if Dn.empty:
            if Fn.empty:
                #return Fn, None, None
                raise IOError('{Fn} and {Dn} missing')
            Dn_grp = None
        else:
            Dn_grp = Dn.groupby(colgrp)
            Dn_keys = list(Dn_grp.groups.keys())
            if Fn.empty:
                Fn_keys = Dn_keys
            else:
                Fn_idx = [item[0] for item in Fn_keys]
                Fn_keys.extend([item for item in Dn_keys
                                if item[0] not in Fn_idx])
        #
        return Fn_keys, Fn_grp, Dn_grp

    #
    def _K_penalty(self, K_stiffness: float,
                   factor: float = 10_000):
        """
        Return:
            Kp : dof[stiffness]
            head_stiff: dict
        """
        K_trans = K_stiffness * factor
        rot_factor = factor * 10
        K_rot = K_stiffness * rot_factor
        col_disp = self._mesh._plane.hdisp
        #
        K_item = {'x': K_trans, 'y': K_trans, 'z': K_trans,
                  'rx': K_rot, 'ry': K_rot, 'rz': K_rot, }

        head_stiff = {key: K_item[key] for key in col_disp}
        Kp = np.array([K_item[key] for key in col_disp])
        #
        return Kp, head_stiff

    #
    def PMT(self, Ke, Fn, Dn,
            sparse: bool = False):
        """
        Penalty Method
        
        Ke : Global stiffness matrix
        Fn : Global force  df
        Dn : Global displacement df
        jbc = Nodes with boundary
        
        Return:
        U : Global node displacement
        """
        jbc = self._mesh.jbc()
        K_factor, solver = self._get_solver(Ke, sparse)
        #
        # TODO : must be better ways to manipulate dfs
        #
        ndof = self._mesh._plane.ndof
        col_force = self._mesh._plane.hforce
        head_force = ['node_name', *col_force]
        col_disp = self._mesh._plane.hdisp
        head_disp = ['node_name', 'node_index', *col_disp]
        col_rename = self._mesh._plane.colrename
        #
        # Select nodes' free DOF
        Fn_bool, Fn_zeros = self._mask_DOF(jbc)
        #
        # group global node load and displacement
        load_keys, Fn_grp, Dn_grp = self._get_FnDn(Fn, Dn)
        #
        # --------------------------------------
        #
        #Ke_sup = np.zeros((ndof, ndof), dtype=np.float64)
        #
        Kp, head_stiff = self._K_penalty(K_factor)
        #
        # --------------------------------------
        #
        Utemp = []
        for key in load_keys:
            Fs_item = Fn_zeros.copy()
            Ke_item = Ke.copy()
            #
            # Load Section
            try:
                Fn_item = Fn_grp.get_group(key)
                Fn_set = Fn_item[head_force].set_index(['node_name'])
                Fs_item.loc[Fn_set.index] = Fn_set.astype('float64')
            except (KeyError, AttributeError):
                pass
            #
            # Displacement Section
            try:
                Dn_item = Dn_grp.get_group(key)
                FDs = Dn_item[head_disp].set_index(['node_name'])
                #
                # calculate force {F} = {D} x Kmod
                FDs[col_disp] = FDs[col_disp].mul(head_stiff)
                FDs.rename(columns=col_rename, inplace=True)
                #
                Fs_item.loc[FDs.index] = Fs_item.loc[FDs.index].add(FDs[col_force])
                #
                Ke_item = self._update_K_FDn(Ke_item=Ke_item,
                                             FDs=FDs,
                                             Kp=Kp,
                                             col_force=col_force)
                #
                #for nodeid, idof in zip(FDs.index.to_list(),
                #                        FDs['node_index'].to_list()):
                #    # DOF
                #    niqi = idof * ndof
                #    niqj = niqi + ndof
                #    #
                #    nload = FDs[col_force].loc[nodeid]
                #    kc = Kp.copy()
                #    kc[nload == 0] = 0
                #    #
                #    Ke_sup.flat[0::ndof + 1] = kc
                #    # update global K matrix with dummy stiffness
                #    Ke_item[niqi:niqj, niqi:niqj] += Ke_sup
            except (KeyError, AttributeError):
                pass
            #
            # get load matrix as vector
            Fs_item = Fs_item.stack()
            Fs_item = Fs_item.loc[Fn_bool]
            # Solve displacements U
            Us = self._Us(Ke_item, Fs_item, solver, jbc, ndof)
            # pack basic load data for df's dumping 
            Utemp.extend([[*key, self._result_name, nname, *Us[x]]
                          for x, nname in enumerate(jbc.index)])
        #
        Us = self._to_df(Utemp)
        return Us

    #
    def _update_K_FDn(self, Ke_item, FDs,
                      Kp: list, ndof: int,
                      col_force: list[str]):
        """
        """
        Ke_sup = np.zeros((ndof, ndof), dtype=np.float64)
        for nodeid, idof in zip(FDs.index.to_list(),
                                FDs['node_index'].to_list()):
            # DOF
            niqi = idof * ndof
            niqj = niqi + ndof
            #
            nload = FDs[col_force].loc[nodeid]
            kc = Kp.copy()
            kc[nload == 0] = 0
            #
            Ke_sup.flat[0::ndof + 1] = kc
            # update global K matrix with dummy stiffness
            Ke_item[niqi:niqj, niqi:niqj] += Ke_sup
        #
        return Ke_item

    #
    # ------------------------------------------
    # tb removed
    #
    def PMTX(self, load_keys, Ke, Fn_grp, Dn_grp,
            sparse: bool = False):
        """
        Penalty Method

        Ke : Global stiffness matrix
        Fn : Global force  df
        Dn : Global displacement df
        jbc = Nodes with boundary

        Return:
        U : Global node displacement
        """
        jbc = self._mesh.jbc()
        K_factor, solver = self._get_solver(Ke, sparse)
        #
        # TODO : must be better ways to manipulate dfs
        #
        ndof = self._mesh._plane.ndof
        col_force = self._mesh._plane.hforce
        head_force = ['node_name', *col_force]
        col_disp = self._mesh._plane.hdisp
        head_disp = ['node_name', 'node_index', *col_disp]
        col_rename = self._mesh._plane.colrename
        #
        # Select nodes' free DOF
        Fn_bool, Fn_zeros = self._mask_DOF(jbc)
        #
        # group global node load and displacement
        # load_keys, Fn_grp, Dn_grp = self._Fn_items(Fn, Dn)
        #
        # --------------------------------------
        #
        Ke_sup = np.zeros((ndof, ndof), dtype=np.float64)
        #
        Kp, head_stiff = self._K_penalty(K_factor)
        #
        # --------------------------------------
        #
        Utemp = []
        for key in load_keys:
            Fs = Fn_zeros.copy()
            Kitem = Ke.copy()
            #
            # Load Section
            try:
                Fn_item = Fn_grp.get_group(key)
                Fn_set = Fn_item[head_force].set_index(['node_name'])
                Fs.loc[Fn_set.index] = Fn_set.astype('float64')
            except (KeyError, AttributeError):
                # print(f'--> Force {key}')
                pass
            #
            # Displacement Section
            try:
                Dn_item = Dn_grp.get_group(key)
                FDs = Dn_item[head_disp].set_index(['node_name'])
                #
                # calculate force {F} = {D} x Kmod
                FDs[col_disp] = FDs[col_disp].mul(head_stiff)
                FDs.rename(columns=col_rename, inplace=True)
                #
                Fs.loc[FDs.index] = Fs.loc[FDs.index].add(FDs[col_force])
                #
                # idx = {name: i for i, name in enumerate(list(FDs), start=1)}
                # node_name = Ds.index.to_list()
                for nodeid, idof in zip(FDs.index.to_list(),
                                        FDs['node_index'].to_list()):
                    # for item in FDs.itertuples():
                    # idof = item.node_index
                    # nodeid = item.Index
                    #
                    # DOF
                    niqi = idof * ndof
                    niqj = niqi + ndof
                    #
                    nload = FDs[col_force].loc[nodeid]
                    kc = Kp.copy()
                    kc[nload == 0] = 0
                    #
                    Ke_sup.flat[0::ndof + 1] = kc
                    # update global K matrix with dummy stiffness
                    Kitem[niqi:niqj, niqi:niqj] += Ke_sup
                    # print('---> ', niqi, niqj)
            except (KeyError, AttributeError):
                # print(f'--> node disp {key}')
                pass
            #
            # get load matrix as vector
            Fs = Fs.stack()
            Fs = Fs.loc[Fn_bool]
            # Solve displacements U
            Us = self._Us(Kitem, Fs, solver, jbc, ndof)
            # pack basic load data for df's dumping
            Utemp.extend([[*key, self._result_name,
                           # Fn_item['load_title'].iloc[0],
                           nname, *Us[x]]
                          for x, nname in enumerate(jbc.index)])
        #
        Us = self._to_df(Utemp)
        return Us

    #
    def _PTM_DnX(self, Kitem, Fs,
                Fn_item, Dn_set, Dn_zeros,
                ndof: int, K_factor: float):
        """
        """
        # set up
        # Dn_index = Dn_zeros.index.isin(Fn_item['node_name'])
        # Ds = Dn_zeros.copy()
        # Ds.loc[Dn_index] = Dn_set.astype('float64')
        #
        Ke_sup = np.zeros((ndof, ndof), dtype=np.float64)
        #
        K_trans = K_factor * 10_000
        K_rot = K_factor * 100_000
        ditem = {'x': K_trans, 'y': K_trans, 'z': K_trans,
                 'rx': K_rot, 'ry': K_rot, 'rz': K_rot, }
        #
        col_rename = self._mesh._plane.colrename
        head_disp = self._mesh._plane.hdisp
        head_stiff = {key: ditem[key] for key in head_disp}
        kp = np.array([ditem[key] for key in head_disp])
        #
        # Displacement check and postprocess
        #
        # calculate force {F} = {D} x Kmod
        FDs = Ds.mul(head_stiff)
        FDs.rename(columns=col_rename, inplace=True)
        # Update global load vector with displacement load
        # print(f'---> {Ds} {Fs}')
        for nodeid, idof in zip(Fn_item['node_name'], Fn_item['node_index']):
            # check if node with displacement
            try:
                nload = FDs.loc[nodeid]  # .to_numpy()
                if not nload.any():
                    continue
                Fs.loc[nodeid] = Fs.loc[nodeid].add(nload, axis='rows')
            except KeyError:
                continue
            #
            # DOF
            niqi = idof * ndof
            niqj = niqi + ndof
            #
            kc = kp.copy()
            kc[nload == 0] = 0
            #
            Ke_sup.flat[0::ndof + 1] = kc
            # update global K matrix with dummy stiffness
            Kitem[niqi:niqj, niqi:niqj] += Ke_sup
            # print('--->', niqi, niqj)
        #
        return Kitem, Fs

    #
    #
    # ------------------------------------------
    #
    def solve(self,
              sparse: bool = True,
              max_iter: int = 30,
              beam_steps: int = 10):
        """
        Solves the static system by the Direct Stiffness Method (DSM)

        method : banded, frontal
        """
        if self._mesh:
            if self.second_order:
                order = "2nd"
                self._postprocess._Pdelta = True
                run = self.solve_Pdelta
            else:
                order = "1st"
                self._postprocess._Pdelta = False
                run = self.solve_Linear
            #
            print(f"** Solving Linear Static [{order} order] ")
            #
            Un = run(sparse=sparse,
                     max_iter=max_iter)
            #
            self._postprocess.Un.df = Un
        else:
            raise IOError('** error: mesh missing')
        # run postprocessing
        self._postprocess.run(beam_steps=beam_steps)
        #
        return self._postprocess._results
#
