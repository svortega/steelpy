#
# Copyright (c) 2009 steelpy
#
from __future__ import annotations
from dataclasses import dataclass
import time

# package imports
#
from steelpy.utils.math.operations import to_matrix
from steelpy.utils.dataframe.main import DBframework
#
import numpy as np
#from scipy.linalg import cholesky_banded, cho_solve_banded
from scipy.sparse.linalg import spsolve

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

@dataclass
class LinearStatic:
    """ """
    __slots__ = ['_mesh', '_result_name', '_nonlinear']

    def __init__(self, mesh, result_name:int|str,
                 nonlinear: bool = False):
        self._mesh = mesh
        self._result_name = result_name
        self._nonlinear = nonlinear
    #
    # ------------------------------------------
    #
    def solve(self,
              sparse: bool = False,
              max_iter: int | None = None):
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
        Un, Ke, Fn, Dn = self._linear(load_level='basic',
                                      sparse=sparse)
        #
        # ------------------------------
        # Post-processing load combinations
        Un = self._update_Un_comb(Un=Un)
        #
        #self._postprocess.Un = Un
        #
        uptime = time.time() - start_time
        print(f"** {{F}} = [Ke] {{U}} Solution: {uptime:1.4e} sec")
        return Un
    #
    # ------------------------------------------
    #
    def _linear(self, load_level: str,
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
        # ------------------------------
        #
        Fn = self._mesh.Fn()
        Fn = Fn[Fn['load_level'] == load_level]
        Dn = self._mesh.Dn()
        Dn = Dn[Dn['load_level'] == load_level]
        #if load_level in ['comb', 'combination']:
        #   load = self._mesh._load._combination
        #   Fn = Fn[Fn['load_level']==load_level]
        #else:
        #    load = self._mesh._load._basic
        # get load data
        #Dn = load.Dnt()
        #Fn = load.Fnt()
        # ------------------------------
        # linear solution
        Un = self.PMT(Ke=Ke,
                      Fn=Fn, Dn=Dn,
                      sparse=sparse)
        return Un , Ke, Fn, Dn
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
        comb2basic = load_combination.to_basic()
        comb_grp = comb2basic.groupby('load_name')
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
                   'load_title', 'system',
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
    # ------------------------------------------
    #
    def PMT(self, Ke, Fn, Dn,
            sparse: bool = False):
        """
        Penalty Method

        Ke : Global stiffness matrix
        Fn : Global force  df
        Dn : Global displacement df
        sparse = True/False

        Return:
        U : Global node displacement
        """
        jbc = self._mesh.jbc()
        K_factor, solver = self._solver(Ke, sparse)
        #
        # TODO : must be better ways to manipulate dfs
        # --------------------------------------
        #
        ndof = self._mesh._plane.ndof
        col_force = self._mesh._plane.hforce
        head_force = ['node_name', *col_force]
        col_disp = self._mesh._plane.hdisp
        head_disp = ['node_name', 'node_index', *col_disp]
        col_rename = self._mesh._plane.colrename
        #
        # --------------------------------------
        # Select nodes' free DOF
        Fn_bool, Fn_zeros = self._mask_DOF(jbc)
        # group global node load and displacement
        load_keys, Fn_grp, Dn_grp = self._get_FnDn(Fn, Dn)
        #
        # --------------------------------------
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
                                             Kp=Kp, ndof=ndof,
                                             col_force=col_force)
            except (KeyError, AttributeError):
                pass
            #
            # get load matrix as vector
            Fs_item = Fs_item.stack()
            Fs_item = Fs_item.loc[Fn_bool]
            # Solve displacements {U} = [k]{Fu}
            Us = self._KeFu(Ke_item, Fs_item, solver, jbc)
            # pack basic load data for df's dumping
            Utemp.extend([[*key, self._result_name, item, *Us[x]]
                          for x, item in enumerate(jbc.index)])
        #
        # --------------------------------------
        #
        header = ['load_name', 'load_id', 'load_level',
                  'system', 'load_title', 'mesh_name',
                  'result_name', 'node_name',
                  *col_disp]
        Us = list2df(data=Utemp, header=header)
        return Us
    # ------------------------------------------
    #
    def _solver(self, Ke, sparse: bool):
        """
        Ke : stiffness matrix
        sparse: True/False

        Returns:
            solver : sparse/
            Kfactor : max value in matrix
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
    def _KeFu(self, Ke, Fs, solver,
              jbc: list):
        """ """
        dof = self._mesh._plane.ndof
        # jbc Matrix to vector
        jbcflat = jbc.stack()
        # FIXME: Matrix condensed
        K11, K12 = self._mesh._matrix_partition(Ke)
        # Solve displacements U
        #try:
        Us = iter(solver(K11, Fs))
        #except Warning:
        #    print('-->')
        # reshape vector in matrix form [row, col]
        Us = [next(Us) if item != 0 else item
              for item in jbcflat]
        # bak to matrix form
        Us = to_matrix(Us, dof)
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
                      'load_level', 'system',
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

@dataclass
class PDelta:
    """ """
    __slots__ = ['_mesh', '_static']

    def __init__(self, mesh, static:LinearStatic):
        self._mesh = mesh
        self._static = static
    #
    def TCM(self, sparse: bool = False,
            max_iter: int | None = None):
        """
        Two Cycles Iterative Method
        [Stability Design of Steel Frames - Chen & Lui]

        """
        # ------------------------------
        # start process
        # ------------------------------
        start_time = time.time()
        #
        # Get global nodal load
        col_grp = ['load_name', 'load_id',
                   'load_level', 'system', 'load_title',
                   'mesh_name']
        #
        # ------------------------------
        # Step 1 : linear Static solution
        # ------------------------------
        #
        Un, Ke, Fn, Dn = self._static._linear(load_level='combination',
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
            Un_temp.append(self._static.PMT(Ke=ke,
                                            Fn=Fn_set,
                                            Dn=Dn_set,
                                            sparse=sparse))
            # ------------------------------
        #
        # Results to dataframe
        db = DBframework()
        Un = db.concat(Un_temp, ignore_index=True)
        # ------------------------------
        # end process
        # ------------------------------
        #
        uptime = time.time() - start_time
        print(f"** {{F}} = [Kt] {{U}} Solution: {uptime:1.4e} sec")
        return Un
#
#
#
# ---------------------------------------------
#
def list2df(data:list[list], header:list):
    """ dataframe"""
    db = DBframework()
    dftemp = db.DataFrame(data=data, columns=header, index=None)
    return dftemp
#
#