#
# Copyright (c) 2009 steelpy
#
from __future__ import annotations
from dataclasses import dataclass
import time
from operator import sub, add
from math import dist

# package imports
#
#from steelpy.utils.math.operations import to_matrix
from steelpy.utils.dataframe.main import DBframework
from steelpy.ufo.utils.beam import BeamItemDeformed
from steelpy.utils.math.pymatsolver import Solver
from steelpy.utils.math.rotation import R_from_drot, transform_unit
#


#
import numpy as np
#from scipy.linalg import cholesky_banded, cho_solve_banded
from scipy.sparse.linalg import spsolve
#from scipy.sparse import coo_matrix
from scipy.sparse import dok_matrix

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
class StaticBasic:
    __slots__ = ['_mesh', '_result_name', '_nonlinear']
    
    def __init__(self, mesh, result_name:int|str,
                 nonlinear: bool):
        self._mesh = mesh
        self._result_name = result_name
        self._nonlinear = nonlinear
    #
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
            #Kfactor = np.max(Ke.data.max())
            Kfactor = np.max(Ke.tolil().data.max())
        else:
            solver = solver_np
            Kfactor = Ke.max()

        return Kfactor, solver
    #
    # ------------------------------------------
    #
    def linear(self, load_level: str):
        """
        Linear Static Analysis

        load_level: basic/combination
        columns
        sparse: matrix operation (True/False)
        """
        # ---------------------------------
        # 
        jbc = self._mesh.jbc()
        # ---------------------------------
        # Assembly stiffness matrix [Ke]
        Ke = self._mesh.Ke()
        #
        # ---------------------------------
        # get nodal force and prescribed displacements
        Fn = self._mesh.Fn()
        Fn = Fn[Fn['load_level'] == load_level]
        Dn = self._mesh.Dn()
        #Dn = Dn[Dn['load_level'] == load_level]
        #
        # ---------------------------------
        # Solution by DSM
        Un, Rn = self.DSM(Ke=Ke, Fn=Fn,
                          Dn=Dn, jbc=jbc)
        #
        return Un, Rn, Ke, Fn, Dn
    #
    # ------------------------------------------
    #
    def DSM(self, Ke,
            Fn: DBframework, Dn: DBframework,
            jbc: DBframework, factor: float = 1.0):
        """
        Direct Stiffness Method (DSM)

        Ke : Global stiffness matrix
        Fn : Global force  df
        Dn : Global displacement df

        Return:
        {Un} : Global node displacements
        {Rn} : Global node reactions
        """
        #jbc = self._mesh.jbc()
        #K_factor, solver = self._solver(Ke, sparse)
        #
        # TODO : must be better ways to manipulate dfs
        # --------------------------------------
        #
        #ndof = self._mesh._plane.ndof
        col_force = self._mesh._plane.hforce
        head_force = ['node_name', *col_force]
        col_disp = self._mesh._plane.hdisp
        #head_disp = ['node_name', 'node_index', *col_disp]
        head_disp = ['node_name', *col_disp]
        #col_rename = self._mesh._plane.colrename
        #
        header = ['load_name', 'load_id', 'load_level',
                  'system', 'load_title',
                  'mesh_name','result_name']        
        #
        # --------------------------------------
        # Select nodes' free DOF
        #Fn_bool, Fn_zeros = self._mask_F(jbc)
        #Dn_bool, Dn_zeros = self._mask_D(jbc)
        # group global node load and displacement
        load_keys, Fn_grp, Dn_grp = self._get_FnDn(Fn, Dn)
        #load_keys, Fn_grp = self._get_FnDn(Fn)
        #
        # --------------------------------------
        #
        #Kp, head_stiff = self._K_penalty(K_factor, factor=1)
        #
        Dn_zeros = self._zeros_df(df=jbc)
        #
        renane=self._mesh._plane.colrename
        Fn_zeros = self._zeros_df(df=jbc, rename=renane)
        #
        # --------------------------------------
        #
        Rtemp = []
        Utemp = []
        for key in load_keys:
            Ke_item = Ke.copy()
            Fs_item = Fn_zeros.copy()
            #
            # --------------------------------------
            # Node Displacement Section
            # --------------------------------------
            Dn_item = Dn_zeros.copy()
            Dn_set = Dn_grp.get_group((key[-1], ))
            Dn_set = Dn_set[head_disp].set_index(['node_name'])
            Dn_item.loc[Dn_set.index] = Dn_set.astype('float32')
            #Dn_item *= factor
            #
            # --------------------------------------
            # Nodal Load Section
            # --------------------------------------
            try:
                Fn_set = Fn_grp.get_group(key)
                Fn_set = Fn_set[head_force].set_index(['node_name'])
                Fs_item.loc[Fn_set.index] = Fn_set.astype('float32')
                #Fs_item *= factor
            except (KeyError, AttributeError):
                pass            
            #
            #try:
            #    Dn_item = Dn_grp.get_group(key)
            #    Dn_item = Dn_item[head_disp].set_index(['node_name'])
            #    #Dn_item.loc[Dn_set.index] = Dn_set.astype('float64')
            #except (KeyError, AttributeError):
            #    pass
            #
            # --------------------------------------
            #
            #Ke_item, Fs_item = self._get_KeFs(key, Fn_grp, Dn_grp,
            #                                  Fs_item, Ke_item,
            #                                  Kp, head_stiff)
            # --------------------------------------
            # get load matrix as vector
            #Fs_item = Fs_item.stack()
            #Fs_item = Fs_item.loc[Fn_bool]
            #
            #Dn_item = Dn_item.stack()
            #Dn_item = Dn_item.loc[Dn_bool]
            #
            #1 / 0
            #
            # --------------------------------------
            # Solve displacements {U} = [k]{Fu}
            # --------------------------------------
            Us, Rf = self._DMCiRC(Ke_item, Fs_item,
                                  Dn_item, jbc=jbc,
                                  factor=factor)
            #
            #Dn_item[col_disp].add(Us)
            #
            #Dnn = Dn_item.stack()
            #Us[Dn_bool] = np.add(Us[Dn_bool], Dnn.values)
            #Us = to_matrix(Us, ndof)
            # --------------------------------------
            # Add load reference to displacement
            # --------------------------------------
            #Utemp.extend([[*key, self._result_name, item, *Us[x]]
            #              for x, item in enumerate(jbc.index)])
            Us.reset_index(inplace=True)
            Us[header] = [*key, self._result_name]
            Utemp.append(Us)
            #
            # --------------------------------------
            # Add load reference to reactions
            # --------------------------------------
            Rf.reset_index(inplace=True)
            Rf[header] = [*key, self._result_name]
            Rtemp.append(Rf)
        #
        db = DBframework()
        # --------------------------------------
        #
        #header = ['load_name', 'load_id', 'load_level',
        #          'system', 'load_title', 'mesh_name',
        #          'result_name', 'node_name',
        #          *col_disp]
        header_df = [*header, 'node_name', *col_disp]
        #Us = list2df(data=Utemp, header=header_df)
        #
        Us = db.concat(Utemp)
        Us = Us[header_df]
        #
        # --------------------------------------
        #
        Rf = db.concat(Rtemp)
        header_df = [*header, 'node_name', *col_force]
        Rf = Rf[header_df]
        #
        return Us, Rf
    #
    def _DMCiRC(self, Ke, Fs:DBframework, Ds:DBframework,
                jbc:DBframework, factor: float = 1.0):
        """
        Direct Stiffness Method including Restrained Coordinates (DSMiRC).
        Matrix Analysis of Structure - Aslam Kassimali - Chapter 9

        [Ke] : Global stiffness matrix
        Fn : Global force  df
        Dn : Global displacement df

        Return:
        Un : Global node displacements df
        Rn : Global node reactions df

        """
        #dof = self._mesh._plane.ndof
        # jbc Matrix to vector
        #jbcflat = jbc.stack()
        # FIXME: Matrix condensed
        #
        D1, D2 = self._mesh._D_partition(jbc)
        K11, K12, K21, K22 = self._K_partition(Ke, D1, D2)
        #1 / 0
        # Solve displacements U
        #try:
        #for idx in Ds.index:
        #    Fs.loc[idx] = np.subtract(Fs.loc[idx], K12 @ Ds.loc[idx])
        #
        #Ds3 = Ds.stack().values
        P = Fs.stack(future_stack=True) * factor
        #P.iloc[D1] = float(0)
        #
        Dr = Ds.stack(future_stack=True) * factor
        #Dr = Dr[D2]
        #K12Dr = K12 @ Dr[D2]
        #F = np.subtract(Fs, K12 @ Ds)
        Pf = np.subtract(P.iloc[D1], K12 @ Dr.iloc[D2]) 
        #F = self._vector_partition(F, jbc)
        #
        #Us = iter(solver(K11, Pf))
        #except Warning:
        #    print('-->')
        # reshape vector in matrix form [row, col]
        #Us = np.array([next(Us) if item != 0 else item
        #               for item in jbcflat])
        #
        Dr.replace(np.nan, 0, inplace=True)
        #Dr.iloc[D1] += solver_sparse(K11, Pf)
        Kinv = Solver(K11.tocsc())
        Dr.iloc[D1] += Kinv * Pf
        #
        # add prescribed displacement
        #Us.iloc[D2] = Us.iloc[D2].add(Dr.iloc[D2])
        #
        # Reactions
        #R = np.add(K21 @ Us[D1] + K22 @ np.nan_to_num(Dr[D2]), P[D2])
        #R = Fs.stack(future_stack=True).mul(-float(1))
        #R.iloc[D2] = R.iloc[D2].add(K21 @ Us[D1] + K22 @ np.nan_to_num(Dr[D2])).astype('float64')
        #R.iloc[D2] += (K21 @ Us[D1] + K22 @ np.nan_to_num(Dr[D2])) #.astype('float64')
        #
        #P = P.copy().mul(float(-1))
        P *= float(-1)
        P.iloc[D2] += (K21 @ Dr.iloc[D1]
                       + K22 @ np.nan_to_num(Dr.iloc[D2])).astype('float32')
        P.iloc[D1] = float(0)
        #
        P = P.unstack() #.reset_index()
        Dr = Dr.unstack() #.reset_index()
        # bak to matrix form
        #Us = to_matrix(Us, dof)
        return Dr, P
    #
    # ------------------------------------------
    #
    def _K_partition(self, Km, D1, D2):
        """
        Partitions a matrix into sub-matrices based on degree of freedom boundary conditions
        
        Km: The un-partitioned matrix (or vector) to be partitioned.

        Return:
            m1 :
            m2 :
        """
        #jbc = self._mesh.jbc()
        #jbcc = jbcflat.values
        #jbc = jbc.stack()
        #
        # List of the indices for degrees of freedom with unknown displacements.
        #D1 = [i for i, item in enumerate(jbc)
        #      if item != 0]
        # List of the indices for degrees of freedom with known displacements.
        #D2 = [i for i, item in enumerate(jbc)
        #      if item == 0]
        #
        #1 / 0
        # 1D vectors
        #if Km.shape[1] == 1:
        #    # Partition the vector into 2 subvectors
        #    m1 = Km[D1, :]
        #    m2 = Km[D2, :]
        #    return m1, m2
        #
        # 2D matrices
        #else:
        # Partition the matrix into 4 submatrices
        m11 = Km[D1, :][:, D1] # dof x dof
        m12 = Km[D1, :][:, D2] # dof x R
        m21 = Km[D2, :][:, D1] # R x dof
        m22 = Km[D2, :][:, D2] # R x R
        return m11, m12, m21, m22
    #
    def _D_partitionX(self, jbc:DBframework):
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
    def _zeros_df(self, df:DBframework, rename: list | None=None)->DBframework:
        """ """
        df_zeros = df.copy().astype('float32')
        if rename:
            df_zeros.rename(columns=rename, inplace=True)
        #
        df_zeros.iloc[:] = float(0.0)
        return df_zeros
    #
    # ------------------------------------------
    # Penalty Method
    # ------------------------------------------
    #
    def PMTXX(self, Ke, Fn, Dn,
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
            Us = self._DMCiRC(Ke_item, Fs_item, solver, jbc)
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
    #
    #
    def _penalty_function(self, key, Dn_grp,
                          Ke_item, Fs_item,
                          Kp, head_stiff):
        """ """
        ndof = self._mesh._plane.ndof
        col_force = self._mesh._plane.hforce
        col_rename = self._mesh._plane.colrename
        #
        col_disp = self._mesh._plane.hdisp
        head_disp = ['node_name', 'node_index', *col_disp]
        #
        try:
            Dn_item = Dn_grp.get_group(key)
            FDs = Dn_item[head_disp].set_index(['node_name'])
            #
            # calculate force {F} = {D} x Kmod
            FDs[col_disp] = FDs[col_disp].mul(head_stiff)
            FDs.rename(columns=col_rename, inplace=True)
            #
            Fs_item.loc[FDs.index] = Fs_item.loc[FDs.index].add(FDs[col_force])
            #Fs_item.loc[FDs.index] = Fs_item.loc[FDs.index].subtract(FDs[col_force])
            #
            Ke_item = self._update_K_FDn(Ke_item=Ke_item,
                                         FDs=FDs,
                                         Kp=Kp, ndof=ndof,
                                         col_force=col_force)
        except (KeyError, AttributeError):
            pass
        #
        return Ke_item, Fs_item
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
    def _get_KeFs(self, key, Fn_grp, Dn_grp,
                  Fs_item, Ke_item, Kp, head_stiff):
        """ """
        # elements = self._mesh._elements
        # set headers
        col_force = self._mesh._plane.hforce
        head_force = ['node_name', *col_force]
        ndof = self._mesh._plane.ndof
        col_force = self._mesh._plane.hforce
        col_rename = self._mesh._plane.colrename
        #
        col_disp = self._mesh._plane.hdisp
        head_disp = ['node_name', 'node_index', *col_disp]
        #
        # Load Section
        try:
            Fn_item = Fn_grp.get_group(key)
            Fn_set = Fn_item[head_force].set_index(['node_name'])
            Fs_item.loc[Fn_set.index] = Fn_set.astype('float32')
        except (KeyError, AttributeError):
            pass
        #
        # Displacement Section
        Ke_item, Fs_item = self._penalty_function(key, Dn_grp,
                                                  Ke_item, Fs_item,
                                                  Kp, head_stiff)
        # try:
        #    Dn_item = Dn_grp.get_group(key)
        #    FDs = Dn_item[head_disp].set_index(['node_name'])
        #    #
        #    # calculate force {F} = [K] {u}
        #    for nodeid, idof in zip(FDs.index.to_list(),
        #                            FDs['node_index'].to_list()):
        #        ndisp = FDs[col_disp].loc[nodeid]
        #        # DOF
        #        niqi = idof * ndof
        #        niqj = niqi + ndof
        #        K_item = Ke_item[niqi:niqj, niqi:niqj]
        #        #bname = 5
        #        #beam = elements[bname]
        #        #K_item = beam.Ke(plane2D=self._mesh._plane.plane2D)
        #        #
        #        #ndisp = np.concatenate(([0]*6, ndisp))
        #        #
        #        #FDs[col_disp] = FDs[col_disp].mul(K_item)
        #        #FDs.rename(columns=col_rename, inplace=True)
        #        Fb = K_item @ ndisp
        #        #Fb[ndisp == 0] = 0
        #        #1/0
        #        Fs_item.loc[FDs.index] = Fs_item.loc[FDs.index].add(Fb)
        # except (KeyError, AttributeError):
        #    pass

        return Ke_item, Fs_item

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
    # ------------------------------------------
    #
    #
    def _mask_F(self, jbc, rename: bool):
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
        return dfbool.astype('bool'), dfjbc.astype('float32')
    #
    def _mask_D(self, jbc):
        """ """
        # remove rows with zeros
        dfjbc = jbc.copy()
        #dfjbc = dfjbc[jbc.any(axis=1)]
        dfjbc = dfjbc[jbc.eq(0).any(axis=1)]
        #dfjbc = dfjbc.replace(float(0.0), np.nan)
        #
        dfbool = dfjbc.copy()
        dfbool = (dfbool == 0)
        dfbool = dfbool.stack(future_stack=True)
        # Copy dataframe
        # dfzeros = dfjbc.copy()
        # dfzeros.iloc[:] = float(0.0)
        #
        dfjbc.iloc[:] = float(0.0)
        # for col in dfjbc.columns:
        #    dfjbc[col].values[:] = float(0.0)
        #
        return dfbool.astype('bool'), dfjbc.astype('float32')
    #
    #
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
            raise IOError('boundary conditions missing')
        #    if Fn.empty:
        #        #return Fn, None, None
        #        raise IOError('{Fn} and {Dn} missing')
        #    Dn_grp = None
        else:
            Dn_grp = Dn.groupby(['mesh_name'])
        #    Dn_grp = Dn.groupby(colgrp)
        #    Dn_keys = list(Dn_grp.groups.keys())
        #    if Fn.empty:
        #        Fn_keys = Dn_keys
        #    else:
        #        Fn_idx = [item[0] for item in Fn_keys]
        #        Fn_keys.extend([item for item in Dn_keys
        #                        if item[0] not in Fn_idx])
        #
        return Fn_keys, Fn_grp, Dn_grp
    #
    #
    # ------------------------------------------
    #
#
#
@dataclass
class LinearStatic(StaticBasic):
    """ """
    __slots__ = ['_mesh', '_result_name', '_nonlinear']

    def __init__(self, mesh, result_name:int|str,
                 nonlinear: bool = False):
        super().__init__(mesh=mesh,
                         result_name=result_name,
                         nonlinear=nonlinear)
    #
    # ------------------------------------------
    #
    def solve(self, max_iter: int | None = None):
        """
        Linear Static Analysis

        Input:
        Ke  : Global stiffness matrix
        Fn  : FER Node load dataframe
        jbc : Node boundary condition dataframe

        Return:
        Udf : Node displacement global system dataframe
        """
        print(f"** Solving Linear Static [1st order] ")
        start_time = time.time()
        # ------------------------------
        if self._nonlinear:
            raise NotImplementedError("Non-linear solution not implemented")
        else:
            # Linear solution
            # ------------------------------
            # Step 1 : linear solution basic load
            Un, Rn, Ke, Fn, Dn = self.linear(load_level='basic')
            #
            # ------------------------------
            # Step 2 : Update Un load combinations
            #          with superposition method
            header = self._mesh._plane.hdisp
            Un = self.update_combination(Un, header)
            header = self._mesh._plane.hforce
            Rn = self.update_combination(Rn, header)
        #
        # ------------------------------
        #self._postprocess.Un = Un
        #
        uptime = time.time() - start_time
        print(f"** {{F}} = [Ke] {{U}} Solution: {uptime:1.4e} sec")
        return Un, Rn
    #    
    #
    def update_combination(self, Un: DBframework,
                           header: list[str]) -> DBframework:
        """
        Update node displacements to include lcomb

        Un : node displacement dataframe
        Return:
            Un + load combinations
        """
        #disp_head = self._mesh._plane.hdisp
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
                comb_item.loc[:, header] *= row.factor
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
            df_temp = df_temp.groupby(columns, as_index=False)[header].sum()
            Un = db.concat([Un, df_temp], ignore_index=True)
        except ValueError:
            pass
        #
        return Un

    #
    #


@dataclass
class PDelta(StaticBasic):
    """ """
    __slots__ = ['_mesh', '_result_name', '_nonlinear', 'tol']

    def __init__(self, mesh,
                 result_name:int|str,
                 nonlinear: bool = False, 
                 method: str = 'simple_step'):
        """
        nonlinear: Nonlinear Analysis True/False
        method: solution method: simple_step/iterative
        """
        super().__init__(mesh=mesh,
                         result_name=result_name,
                         nonlinear=nonlinear)        
        self.tol = {'u': 1e-2, 'r': 1e-2}
    #
    # ------------------------------------------
    #
    def solve(self, 
              max_iter: int = 30,
              method: str = 'simple_step'):
        """ """
        print(f"** Solving Linear Static [2nd order] ")
        start_time = time.time()
        #
        # ------------------------------
        #
        if method.lower() in ['simple_step']:
            Un, Rn = self.TCM()
        else:
            #Un, Rn = self.LCM_Euler(max_iter=max_iter)
            Un, Rn = self.LCM_RungeKutta(max_iter=max_iter)
            
        #
        # ------------------------------
        #
        uptime = time.time() - start_time
        print(f"** {{F}} = [Ke] {{U}} Solution: {uptime:1.4e} sec")
        return Un, Rn        
    #
    # ------------------------------------------
    #
    def TCMxx(self, sparse: bool = False,
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
        Un, Ke, Fn, Dn = self.linear(load_level='combination',
                                     sparse=sparse)
        #
        # ------------------------------
        # Step 2 : Pdelta solution
        # ------------------------------
        #
        #Un_grp = Un.groupby(col_grp)
        ##
        #Un_temp = []
        #for key, Un_step in Un_grp:
        #    # ------------------------------
        #    if Fn.empty:
        #        Fn_set = Fn
        #    else:
        #        Fn_set = Fn.groupby(col_grp).get_group(key)
        #    #
        #    if Dn.empty:
        #        Dn_set = Dn
        #    else:
        #        Dn_set = Dn.groupby(col_grp).get_group(key)
        #    # ------------------------------
        #    # Assembly matrix including node force
        #    ke = self._mesh.Kt(Dn=Un_step)
        #    # ------------------------------
        #    Un_temp.append(self._static.PMT(Ke=ke,
        #                                    Fn=Fn_set,
        #                                    Dn=Dn_set,
        #                                    sparse=sparse))
        #    # ------------------------------
        #
        #
        # Results to dataframe
        #db = DBframework()
        #Un = db.concat(Un_temp, ignore_index=True)
        #
        Un = self._loop_PMT(Un, Fn, Dn,
                            col_grp=col_grp,
                            sparse=sparse)
        #
        # ------------------------------
        # end process
        # ------------------------------
        #
        uptime = time.time() - start_time
        print(f"** {{F}} = [Kt] {{U}} Solution: {uptime:1.4e} sec")
        return Un
    #
    #
    def TCM(self):
        """
        Two Cycles Iterative Method
        [Stability Design of Steel Frames - Chen & Lui]

        """
        # ------------------------------
        # start process
        # ------------------------------
        start_time = time.time()
        #
        #
        # ------------------------------
        # Step 1 : linear Static solution
        # ------------------------------
        #
        Un, Rn, Ke, Fn, Dn = self.linear(load_level='combination')
        #
        # ------------------------------
        # Step 2 : Pdelta solution
        # ------------------------------
        #
        #elements = self._mesh._elements
        #Pdelta = False
        #Fb = self.solve_beam_end(elements, Un, Rn, Pdelta)
        #
        # ---------------------------------
        # 
        jbc = self._mesh.jbc()        
        #
        Un, Rn = self._loop_PMT(Un, Rn, Fn, Dn, jbc)
        #
        # ------------------------------
        # end process
        # ------------------------------
        #
        uptime = time.time() - start_time
        print(f"** {{F}} = [Kt] {{U}} Solution: {uptime:1.4e} sec")
        return Un, Rn
    #
    #
    def _loop_PMT(self, Un, Rn, Fn, Dn, jbc):
        """
        """
        elements = self._mesh._elements
        #jbc = self._mesh.jbc()
        col_disp = self._mesh._plane.hdisp
        head_disp = ['node_name', *col_disp]        
        #
        #Fn.drop(columns=['system', 'load_title'], inplace=True)
        # Get global nodal load
        col_grp = ['load_name', 'load_id',
                   'load_level', 'system', 'load_title',
                   'mesh_name']        
        #
        Un_grp = Un.groupby(col_grp)
        Fn_grp = Fn.groupby(col_grp)
        Rn_grp = Rn.groupby(col_grp)
        Dn_grp = Dn.groupby(['mesh_name'])
        #
        head_Fb = ['element_name',
                   #'L', 'uv', 'R',
                   'Fb_local']
        #
        Un_temp = []
        Rn_temp = []
        for key, Un_step in Un_grp:
            #Fb_set = Fb_step[head_Fb].set_index('element_name')
            Un_set = Un_step[head_disp].set_index(['node_name'])
            # ------------------------------
            # get node displacements
            Dn_set = Dn_grp.get_group((key[-1], ))
            # ------------------------------
            # get nodal force
            #try:
            Fn_set = Fn_grp.get_group(key)
            #except AttributeError:
            #    1 / 0
                #Fn_set = Fn_grp.get_group(key)
            #
            #try:
            #Rn_set = Rn_grp.get_group(key)
            #    #Dn_set = Dn_grp.get_group(key)
            #except AttributeError:
            #    1 / 0
            #    Dn_set = Dn_grp
            #
            # ------------------------------
            #
            Kt, Fb = self.solve_beam_end(elements, Un_set, #Rn_set,
                                         jbc, Pdelta=True)
            #Fb_set = Fb[head_Fb].set_index('element_name')
            #
            # ------------------------------
            # Assembly tangent matrix [Kt]
            #Kt, Rb = self._mesh.Kt_R(Fb=Fb_set)
            #Kg, Rb = self._mesh.Kg(Fb=Fb_set)
            #Kt = Ke + Kg
            #
            # ------------------------------
            # Solve
            Us, Rs = self.DSM(Ke=Kt, Fn=Fn_set,
                              Dn=Dn_set, jbc=jbc)
            # ------------------------------
            Un_temp.append(Us)
            Rn_temp.append(Rs)
        #
        # Results to dataframe
        db = DBframework()
        Ut = db.concat(Un_temp, ignore_index=True)
        Rt = db.concat(Rn_temp, ignore_index=True)
        return Ut, Rt
    #
    # ------------------------------------------
    #
    def IterMethodXX(self, sparse: bool = False,
                   max_iter:int = 30):
        """ """
        # ------------------------------
        # start process
        # ------------------------------
        #start_time = time.time()
        #
        ndof = self._mesh._plane.ndof
        col_rename = self._mesh._plane.colrename
        result_name = self._static._result_name
        #
        # --------------------------------------        
        #
        #col_force = ['Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']
        col_force = self._mesh._plane.hforce
        # Get global nodal load
        col_grp = ['load_name', 'load_id',
                   'load_level', 'system', 'load_title',
                   'mesh_name']
        #
        head_force = ['node_name', *col_force]
        #
        col_disp = self._mesh._plane.hdisp
        head_disp = ['node_name', 'node_index', *col_disp]
        #
        # ------------------------------
        # Step 1 : linear Static solution
        # ------------------------------
        #
        load_level='combination'
        Un, Ke, Fn, Dn = self._static.linear(load_level=load_level,
                                             sparse=sparse)
        #
        #
        jbc = self._mesh.jbc()
        K_factor, solver = self._static._solver(Ke, sparse)
        #
        # --------------------------------------
        #
        Kp, head_stiff = self._static._K_penalty(K_factor)
        Fn_bool, Fn_zeros = self._static._mask_DOF(jbc)
        Dn_zeros = Fn_zeros.copy()
        Dn_zeros.rename(columns=self._mesh._plane.colrename_f2u,
                        inplace=True)
        #        
        #
        # ------------------------------
        # Step 2 : Iter solution
        # ------------------------------
        #
        dfheader = ['load_name', 'load_id', 'load_level',
                    'system', 'load_title', 'mesh_name',
                    'result_name', 'node_name',
                    *col_disp]        
        #
        head_disp2 = ['node_name',  *col_disp]
        #1 / 0
        Un_grp = Un.groupby(col_grp)
        #
        for key, Un_step in Un_grp:
            Fs_item = Fn_zeros.copy()
            Dn_item = Dn_zeros.copy()
            Ke_item = Ke.copy()
            # ------------------------------
            try:
                Fn_item = Fn.groupby(col_grp).get_group(key)
                Fn_set = Fn_item[head_force].set_index(['node_name'])
                Fs_item.loc[Fn_set.index] = Fn_set.astype('float64')
            except (KeyError, AttributeError):
                pass            
            #
            try:
                Dn_item = Dn.groupby(col_grp).get_group(key)
                FDs = Dn_item[head_disp].set_index(['node_name'])
                #
                # calculate force {F} = {D} x Kmod
                FDs[col_disp] = FDs[col_disp].mul(head_stiff)
                FDs.rename(columns=col_rename, inplace=True)
                #
                Fs_item.loc[FDs.index] = Fs_item.loc[FDs.index].add(FDs[col_force])
                #
                Ke_item = self.linear._update_K_FDn(Ke_item=Ke_item,
                                                    FDs=FDs,
                                                    Kp=Kp, ndof=ndof,
                                                    col_force=col_force)
            except (KeyError, AttributeError):
                pass            
            # ------------------------------
            #
            # get load matrix as vector
            Fi = Fs_item.stack()
            Fi = Fi.loc[Fn_bool]
            #
            #
            Ke = Ke_item
            #Fi = Fs_item
            #Di = Dn_item
            #
            u = Un_step.copy()
            #u[col_disp] = 0.0
            #ui =  ui[head_disp2].set_index(['node_name'])
            #
            f = Fi.copy()
            #
            #
            # ------------------------------
            #
            # copy previous force level (used to scale residual for convergence check)
            f_prev = f
            # force at load step or time increment
            f = Fi.copy()         
            # force increment
            df = f #- f_prev
            #
            #
            # deform nodes in part given by u => new f_int and K from elements
            # ??
            # total displacement during increment
            du_inc = Un_step.copy()
            du_inc[col_disp] = 0.0
            du_inc = du_inc[head_disp2].set_index(['node_name'])
            #
            # Calculate internal forces and residual force
            # f_int =
            #
            # residual force
            #r = f - f_int
            #
            #1 / 0
            #iter_count = 1
            #while True:
            for iter_count in range(max_iter):
                #
                # ------------------------------
                # Assembly matrix including node force
                Ki = self._mesh.Kt(Dn=u)
                #
                # Iteration, new displacement (NR iteration)
                #
                # Solve displacements {U} = [k]{Fu}
                du = self._static._DMCiRC(Ki, r, solver, jbc)
                #
                du = [[*key, result_name, item, *du[x]]
                      for x, item in enumerate(jbc.index)]
                du = list2df(data=du, header=dfheader)
                du = du[head_disp2].set_index(['node_name'])
                #
                u =  u[head_disp2].set_index(['node_name'])
                u = u + du  # add to u, NR
                du_inc = du_inc + du
                #
                # Update residual
                #
                # deform nodes in part given by u => new f_int and K from elements
                # new internal (stiffness) force
                # calc residual force
                #
                # Check convergence
                converged = is_converged([np.linalg.norm(du), np.linalg.norm(r)],
                                         [self.tol['u'], self.tol['r']], 
                                         scaling=[np.linalg.norm(du_inc), np.linalg.norm(Fi)])                
                #
                if converged:
                    #print(f'- Tension/compression-only analysis converged after ' + str(iter_count) + ' iteration(s)')
                    print(f'analysis converged after {iter_count} iteration(s)')
                    break
                else:
                    #print('- Tension/compression-only analysis did not converge on this iteration')
                    print('- analysis did not converge on this iteration')
                    print('- Stiffness matrix will be adjusted')
                    print('- P-Delta analysis will be restarted')
                    1 / 0
                #
                #
                # Check for tension/compression-only divergence
                if iter_count > max_iter:
                    #divergence = True
                    raise Exception('Model diverged during tension/compression-only analysis')
                #factor = step / max_iter
                #Fni = Fn.copy()
                #Fni[col_force] = Fni[col_force] * factor
                #print(factor, Fni[self._mesh._plane.hforce])
                # Keep track of the number of tension/compression only iterations
                #iter_count += 1
                #
                # Assemble tangent stiffness if a new iteration is needed
                #u.reset_index()
                #r = Fi - F_int
        #
        1 / 0
    #
    #
    #
    def LCM_Euler(self, step: int = 10,
                  max_iter:int = 30):
        """
        Load Control Method
        """
        # ------------------------------
        # start process
        # ------------------------------
        #
        elements = self._mesh._elements
        jbc = self._mesh.jbc()
        #D1, D2 = self._mesh._D_partition(jbc)
        #
        # ------------------------------
        # start process
        # ------------------------------
        #
        # ------------------------------
        # Step 1 : linear Static solution
        # ------------------------------
        #
        #Un, Rn, Ke, Fn, Dn = self.linear(load_level='combination')
        # ---------------------------------
        # 
        jbc = self._mesh.jbc()
        # ---------------------------------
        # Assembly stiffness matrix [Ke]
        #Ke = self._mesh.Ke()
        #
        # ---------------------------------
        # get nodal force and prescribed displacements
        load_level='combination'
        Fn = self._mesh.Fn()
        Fn = Fn[Fn['load_level'] == load_level]
        Dn = self._mesh.Dn()
        #        
        #        
        # ------------------------------
        #        
        elements = self._mesh._elements
        #
        col_disp = self._mesh._plane.hdisp
        head_disp = ['node_name', *col_disp]        
        #
        # Get global nodal load
        col_grp = ['load_name', 'load_id',
                   'load_level', 'system', 'load_title',
                   'mesh_name']        
        #
        #Un_grp = Un.groupby(col_grp)
        Fn_grp = Fn.groupby(col_grp)
        Dn_grp = Dn.groupby(['mesh_name'])
        #
        # ------------------------------
        #
        #col_force = self._mesh._plane.hforce
        #head_force = ['node_name', *col_force]        
        #
        Dn_zeros = self._zeros_df(df=jbc)
        #
        #renane=self._mesh._plane.colrename
        #Fn_zeros = self._zeros_df(df=jbc, rename=renane)        
        #
        steps =  np.linspace(0.0, 1.0, num=step+1, retstep=True)#.astype('float32')
        #time = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
        #time = [1.0]
        #
        # ------------------------------
        #
        Un_temp = []
        Rn_temp = []
        for key, Fn_set in Fn_grp:
            #Un_set = Un_step[head_disp].set_index(['node_name'])
            #Fn_item = Fn_zeros.copy()
            #Dn_item = Dn_zeros.copy()
            # ------------------------------
            # get node displacements
            Dn_set = Dn_grp.get_group((key[-1], ))
            #Dn_set = Dn_set[head_disp].set_index(['node_name'])
            #Dn_item.loc[Dn_set.index] = Dn_set.astype('float32')             
            # ------------------------------
            # get nodal force
            #Fn_set = Fn_grp.get_group(key)
            #Fn_set = Fn_set[head_force].set_index(['node_name'])
            #Fn_item.loc[Fn_set.index] = Fn_set.astype('float32')
            #
            # ------------------------------
            # Solve linear
            #Us, Rs = self.DSM(Ke=Ke, Fn=Fn_set,
            #                  Dn=Dn_set, jbc=jbc)
            #Un_set = Us[head_disp].set_index(['node_name'])
            #
            # ------------------------------
            # Start first loop
            #fn_prev = Fn_zeros.copy()
            un_prev = Dn_zeros.copy()
            #
            #
            for dlambda in steps[0]:
                #
                factor = dlambda
                #
                # ------------------------------
                # get [Kt] matrix
                Kt, Fb = self.solve_beam_end(elements,
                                             Un=un_prev,
                                             jbc=jbc,
                                             Pdelta=True)                
                #
                # ------------------------------
                # Solve
                Us, Rs = self.DSM(Ke=Kt, Fn=Fn_set,
                                  Dn=Dn_set, jbc=jbc,
                                  factor=factor)
                #
                # ------------------------------
                #
                un_prev = Us[head_disp].set_index(['node_name'])
                #
            #
            # ------------------------------
            Un_temp.append(Us)
            Rn_temp.append(Rs)
            #
        #
        header = ['load_name', 'load_id', 'load_level',
                  'system', 'load_title', 'mesh_name',
                  'result_name']          
        # Results to dataframe
        db = DBframework()
        Ut = db.concat(Un_temp, ignore_index=True)
        #Ut = Ut.groupby(header)[col_disp].sum()
        #Ut.reset_index(inplace=True)
        #
        Rt = db.concat(Rn_temp, ignore_index=True)
        #Rt = Rt.groupby(header)[col_force].sum()
        #Rt.reset_index(inplace=True)
        #
        return Ut, Rt
    #
    #
    def LCM_RungeKutta(self, step: int = 10,
                       max_iter:int = 30):
        """
        Load Control Method
        """
        # ------------------------------
        # start process
        # ------------------------------
        #
        elements = self._mesh._elements
        jbc = self._mesh.jbc()
        #D1, D2 = self._mesh._D_partition(jbc)
        #
        # ------------------------------
        # start process
        # ------------------------------
        #
        # ------------------------------
        # Step 1 : linear Static solution
        # ------------------------------
        #
        #Un, Rn, Ke, Fn, Dn = self.linear(load_level='combination')
        # ---------------------------------
        # 
        jbc = self._mesh.jbc()
        # ---------------------------------
        # Assembly stiffness matrix [Ke]
        Kt = self._mesh.Ke()
        #
        # ---------------------------------
        # get nodal force and prescribed displacements
        load_level='combination'
        Fn = self._mesh.Fn()
        Fn = Fn[Fn['load_level'] == load_level]
        Dn = self._mesh.Dn()    
        #        
        # ------------------------------
        #        
        elements = self._mesh._elements
        #
        col_disp = self._mesh._plane.hdisp
        head_disp = ['node_name', *col_disp]        
        #
        # Get global nodal load
        col_grp = ['load_name', 'load_id',
                   'load_level', 'system', 'load_title',
                   'mesh_name']        
        #
        #Un_grp = Un.groupby(col_grp)
        Fn_grp = Fn.groupby(col_grp)
        Dn_grp = Dn.groupby(['mesh_name'])
        #
        mu = 1 / 2.0
        #
        # ------------------------------
        #
        #col_force = self._mesh._plane.hforce
        #head_force = ['node_name', *col_force]        
        #
        Dn_zeros = self._zeros_df(df=jbc)
        #
        #renane=self._mesh._plane.colrename
        #Fn_zeros = self._zeros_df(df=jbc, rename=renane)        
        #
        steps =  np.linspace(0.0, 1.0, num=step+1, retstep=True)#.astype('float32')
        #time = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
        #time = [1.0]
        #
        # ------------------------------
        #
        Un_temp = []
        Rn_temp = []
        for key, Fn_set in Fn_grp:
            # ------------------------------
            # get node displacements
            Dn_set = Dn_grp.get_group((key[-1], ))            
            # ------------------------------
            #
            # ------------------------------
            # Start first loop
            #fn_prev = Fn_zeros.copy()
            un_prev = Dn_zeros.copy()
            um_prev = Dn_zeros.copy()
            #
            factor = 0
            #
            for dlambda in steps[0]:
                #
                factor = dlambda *  mu
                #
                # ------------------------------
                # Solve
                Um, Rm = self.DSM(Ke=Kt,
                                  Fn=Fn_set, Dn=Dn_set,
                                  jbc=jbc,
                                  factor=factor)                
                #
                #                
                um_prev = Um[head_disp].set_index(['node_name'])
                #
                # ------------------------------
                # get [Km] matrix
                Km, Fb = self.solve_beam_end(elements,
                                             Un=um_prev,
                                             jbc=jbc,
                                             Pdelta=True)
                #
                # ------------------------------
                # Solve
                Us, Rs = self.DSM(Ke=Km,
                                  Fn=Fn_set, Dn=Dn_set,
                                  jbc=jbc,
                                  factor=dlambda)
                #
                # ------------------------------
                #
                un_prev = Us[head_disp].set_index(['node_name'])
                #
                # get [Kt] matrix
                Kt, Fb = self.solve_beam_end(elements,
                                             Un=un_prev,
                                             jbc=jbc,
                                             Pdelta=True)
            #
            # ------------------------------
            #
            #
            # ------------------------------
            Un_temp.append(Us)
            Rn_temp.append(Rs)
            #
        #
        #header = ['load_name', 'load_id', 'load_level',
        #          'system', 'load_title', 'mesh_name',
        #          'result_name']          
        # Results to dataframe
        db = DBframework()
        Ut = db.concat(Un_temp, ignore_index=True)
        #Ut = Ut.groupby(header)[col_disp].sum()
        #Ut.reset_index(inplace=True)
        #
        Rt = db.concat(Rn_temp, ignore_index=True)
        #Rt = Rt.groupby(header)[col_force].sum()
        #Rt.reset_index(inplace=True)
        #
        return Ut, Rt
    #     
    #
    def LCM_Euler2(self, max_iter:int = 30):
        """
        Load Control Method
        """
        # ------------------------------
        # start process
        # ------------------------------
        #
        elements = self._mesh._elements
        jbc = self._mesh.jbc()
        #D1, D2 = self._mesh._D_partition(jbc)
        #
        # ------------------------------
        # start process
        # ------------------------------
        #start_time = time.time()
        #
        # ------------------------------
        # Step 1 : linear Static solution
        # ------------------------------
        #
        Un, Rn, Ke, Fn, Dn = self.linear(load_level='combination')
        #        
        # ------------------------------
        #        
        #        
        #
        #result_name = self._result_name
        #load_level='combination'
        #jbc = self._mesh.jbc()
        #ndof = self._mesh._plane.ndof
        col_force = self._mesh._plane.hforce
        head_force = ['node_name', *col_force]
        col_disp = self._mesh._plane.hdisp
        #head_disp = ['node_name', 'node_index', *col_disp]
        head_disp = ['node_name', *col_disp]
        #col_rename = self._mesh._plane.colrename        
        #
        # ---------------------------------
        # Assembly stiffness matrix [Ke]
        # 
        #Ke = self._mesh.Ke()
        #K_factor, solver = self._solver(Ke)
        #
        # --------------------------------------
        #
        #Kp, head_stiff = self._K_penalty(K_factor)        
        #
        # ---------------------------------
        # get nodal force and prescribed displacements
        #Fn = self._mesh.Fn()
        #Fn = Fn[Fn['load_level'] == load_level]
        #Dn = self._mesh.Dn()
        #Dn = Dn[Dn['load_level'] == load_level]        
        #
        # --------------------------------------
        # Select nodes' free DOF
        #Fn_bool, Fn_zeros = self._mask_DOF(jbc)
        #
        Dn_zeros = self._zeros_df(df=jbc)
        #
        renane=self._mesh._plane.colrename
        Fn_zeros = self._zeros_df(df=jbc, rename=renane)       
        #
        # group global node load and displacement
        #load_keys, Fn_grp, Dn_grp = self._get_FnDn(Fn, Dn)
        #
        # --------------------------------------
        #
        time = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
        #
        col_grp = ['load_name', 'load_id',
                   'load_level', 'system', 'load_title',
                   'mesh_name']        
        #
        Un_grp = Un.groupby(col_grp)
        Fn_grp = Fn.groupby(col_grp)
        Rn_grp = Rn.groupby(col_grp)
        Dn_grp = Dn.groupby(['mesh_name'])
        #
        #head_Fb = ['element_name', 'L', 'uv', 'R', 'Fb_local']
        header = ['load_name', 'load_id', 'load_level',
                  'system', 'load_title', 'mesh_name',
                  'result_name']         
        #
        #mu = 1.0 / 2.0
        #dlambda = 1.0 / 3.0
        #time = [1.0]
        #
        #Un_temp = []
        #Rn_temp = []
        for key, Un_step in Un_grp:
            Kt_prev = Ke.copy()
            Fn_item = Fn_zeros.copy()
            Dn_item = Dn_zeros.copy()
            #
            # ------------------------------
            # get node displacements
            Dn_set = Dn_grp.get_group((key[-1], ))
            Dn_set = Dn_set[head_disp].set_index(['node_name'])
            Dn_item.loc[Dn_set.index] = Dn_set.astype('float32')            
            # ------------------------------
            # get nodal force
            Fn_set = Fn_grp.get_group(key)
            Fn_set = Fn_set[head_force].set_index(['node_name'])
            Fn_item.loc[Fn_set.index] = Fn_set.astype('float32')            
            #
            # --------------------------------------
            # get reactions
            #Rn_set = Rn_grp.get_group(key)
            #
            # --------------------------------------
            #
            Un_set = Un_step[head_disp].set_index(['node_name'])
            #Kt, Fb = self.solve_beam_end(elements, Un=Un_set, 
            #                             jbc=jbc, Pdelta=True)
            #
            # --------------------------------------
            # Start first loop
            f_prev = Fn_zeros.copy()
            #un_prev = Dn_item.copy()
            #du = Dn_zeros.copy()
            un_prev = Dn_zeros.copy()
            #rn = Fn_zeros.copy()
            #
            dn_temp = []
            ru_temp = []
            #R = None
            for dlambda in time:
                #
                # --------------------------------------
                # Euler Method
                # --------------------------------------
                #
                # Force increment
                ##f_prev = f.copy()
                f = Fn_item * dlambda
                f_delta = f - f_prev
                #
                #un_prev = un.copy()
                #un = Dn_item * dlambda
                #un_delta = un - un_prev                
                #
                # Solve displacements {U} = [k]{Fu}
                du, ru = self._DMCiRC(Kt_prev, f_delta, un_prev, jbc)
                #
                # get [Kt] 
                Kt, Fb = self.solve_beam_end(elements, Un=du,
                                             jbc=jbc, Pdelta=True)
                #
                f_prev = f
                Kt_prev = Kt
                #un_prev = -du
                #
            #
            # --------------------------------------
            # Add load reference to displacement
            # --------------------------------------
            #
            #duu = du.reset_index()
            du.reset_index(inplace=True)
            du[header] = [*key, self._result_name]
            dn_temp.append(du)
            #
            #un = un[head_disp].set_index(['node_name'])
            #
            # --------------------------------------
            # Add load reference to reactions
            # --------------------------------------
            #
            ru.reset_index(inplace=True)
            ru[header] = [*key, self._result_name]
            ru_temp.append(ru)                
            #
            # --------------------------------------
        #
        # --------------------------------------
        #
        header = ['load_name', 'load_id', 'load_level',
                  'system', 'load_title', 'mesh_name',
                  'result_name','node_name']
        #
        db = DBframework()     
        #1 / 0
        #header_df = [*header, 'node_name']
        Un = db.concat(dn_temp, ignore_index=True)
        Un = Un.groupby(header)[col_disp].sum()
        Un.reset_index(inplace=True)
        #
        #header_df = [*header, 'node_name', *col_force]
        Rf = db.concat(ru_temp, ignore_index=True)
        Rf = Rf.groupby(header)[col_force].sum()
        Rf.reset_index(inplace=True)
        #
        # ------------------------------
        # end process
        # ------------------------------
        #
        return Un, Rf        
    #    
    #
    def LCM(self, max_iter:int = 30):
        """
        Load Control Method
        """
        # ------------------------------
        # start process
        # ------------------------------
        #
        elements = self._mesh._elements
        jbc = self._mesh.jbc() 
        #
        # ------------------------------
        # start process
        # ------------------------------
        #start_time = time.time()
        #
        # ------------------------------
        # Step 1 : linear Static solution
        # ------------------------------
        #
        Un, Rn, Ke, Fn, Dn = self.linear(load_level='combination')
        #        
        # ------------------------------
        #        
        #        
        #
        #result_name = self._result_name
        #load_level='combination'
        #jbc = self._mesh.jbc()
        #ndof = self._mesh._plane.ndof
        col_force = self._mesh._plane.hforce
        head_force = ['node_name', *col_force]
        col_disp = self._mesh._plane.hdisp
        #head_disp = ['node_name', 'node_index', *col_disp]
        head_disp = ['node_name', *col_disp]
        #col_rename = self._mesh._plane.colrename        
        #
        # ---------------------------------
        # Assembly stiffness matrix [Ke]
        # 
        #Ke = self._mesh.Ke()
        #K_factor, solver = self._solver(Ke)
        #
        # --------------------------------------
        #
        #Kp, head_stiff = self._K_penalty(K_factor)        
        #
        # ---------------------------------
        # get nodal force and prescribed displacements
        #Fn = self._mesh.Fn()
        #Fn = Fn[Fn['load_level'] == load_level]
        #Dn = self._mesh.Dn()
        #Dn = Dn[Dn['load_level'] == load_level]        
        #
        # --------------------------------------
        # Select nodes' free DOF
        #Fn_bool, Fn_zeros = self._mask_DOF(jbc)
        #
        Dn_zeros = self._zeros_df(df=jbc)
        #
        renane=self._mesh._plane.colrename
        Fn_zeros = self._zeros_df(df=jbc, rename=renane)       
        #
        # group global node load and displacement
        #load_keys, Fn_grp, Dn_grp = self._get_FnDn(Fn, Dn)
        #
        # --------------------------------------
        #
        time = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
        #
        col_grp = ['load_name', 'load_id',
                   'load_level', 'system', 'load_title',
                   'mesh_name']        
        #
        Un_grp = Un.groupby(col_grp)
        Fn_grp = Fn.groupby(col_grp)
        Rn_grp = Rn.groupby(col_grp)
        Dn_grp = Dn.groupby(['mesh_name'])
        #
        #head_Fb = ['element_name', 'L', 'uv', 'R', 'Fb_local']
        header = ['load_name', 'load_id', 'load_level',
                  'system', 'load_title', 'mesh_name',
                  'result_name']         
        #
        mu = 1.0 / 2.0
        #dlambda = 1.0 / 3.0
        #
        #Un_temp = []
        #Rn_temp = []
        for key, Un_step in Un_grp:
            Kt_prev = Ke.copy()
            Fn_item = Fn_zeros.copy()
            Dn_item = Dn_zeros.copy()
            #
            # ------------------------------
            # get node displacements
            Dn_set = Dn_grp.get_group((key[-1], ))
            Dn_set = Dn_set[head_disp].set_index(['node_name'])
            Dn_item.loc[Dn_set.index] = Dn_set.astype('float32')            
            # ------------------------------
            # get nodal force
            Fn_set = Fn_grp.get_group(key)
            Fn_set = Fn_set[head_force].set_index(['node_name'])
            Fn_item.loc[Fn_set.index] = Fn_set.astype('float32')            
            #
            # --------------------------------------
            # get reactions
            #Rn_set = Rn_grp.get_group(key)
            #
            # --------------------------------------
            #
            #Un_set = Un_step[head_disp].set_index(['node_name'])
            #Kt, Fb = self.solve_beam_end(elements, Un=Un_set, 
            #                             jbc=jbc, Pdelta=True)
            #
            # --------------------------------------
            # Start first loop
            f_prev = Fn_zeros.copy()
            un_prev = Dn_item.copy()
            #du = Dn_zeros.copy()
            #un = Dn_zeros.copy()
            #rn = Fn_zeros.copy()
            #
            dn_temp = []
            ru_temp = []
            #R = None
            for dlambda in time:
                #
                # --------------------------------------
                # Euler Method
                # --------------------------------------
                #
                # Force increment
                ##f_prev = f.copy()
                #f = Fn_item * dlambda
                #f_delta = f - f_prev
                #
                # Solve displacements {U} = [k]{Fu}
                #du, ru = self._DMCiRC(Kt_prev, f_delta, un_prev, jbc)
                #
                # get [Kt] 
                #Kt, Fb = self.solve_beam_end(elements, Un=du,
                #                             jbc=jbc, Pdelta=True)
                #
                #f_prev = f
                #Kt_prev = Kt
                #
                # --------------------------------------
                # Runge-Kutta Mid-point method
                # --------------------------------------
                #
                # Force increment
                #f_prev = f.copy()
                f = Fn_item * dlambda * mu
                f_delta = f - f_prev                
                #
                # Solve displacements {U} = [k]{Fu}
                dum, rum = self._DMCiRC(Kt_prev, f_delta, un_prev, jbc)                
                #
                # get [Kt] 
                Km, Fb = self.solve_beam_end(elements, Un=dum,
                                             jbc=jbc, Pdelta=True)                
                #
                #
                # Predictor step
                f = Fn_item * dlambda
                Pmi = f - f_prev 
                #
                du, ru = self._DMCiRC(Km, Pmi, un_prev, jbc)
                #
                # Corrector step
                Kmi, Fbmi = self.solve_beam_end(elements, Un=du,
                                                jbc=jbc, Pdelta=True)
                #
                Kt_prev = Kmi # + Km
                f_prev = f
                #
                #f = f + f_prev
                #du = du + dum
                #ru = ru + rum
                #
                #
                # --------------------------------------
                # Add load reference to displacement
                # --------------------------------------
                #
                duu = du.reset_index()
                #du.reset_index(inplace=True)
                duu[header] = [*key, self._result_name]
                dn_temp.append(duu)
                #
                #un = un[head_disp].set_index(['node_name'])
                #
                # --------------------------------------
                # Add load reference to reactions
                # --------------------------------------
                #
                ru.reset_index(inplace=True)
                ru[header] = [*key, self._result_name]
                ru_temp.append(ru)                
                #
                # --------------------------------------
        #
        # --------------------------------------
        #
        header = ['load_name', 'load_id', 'load_level',
                  'system', 'load_title', 'mesh_name',
                  'result_name','node_name']
        #
        db = DBframework()     
        #1 / 0
        #header_df = [*header, 'node_name']
        Un = db.concat(dn_temp, ignore_index=True)
        Un = Un.groupby(header)[col_disp].sum()
        Un.reset_index(inplace=True)
        #
        #header_df = [*header, 'node_name', *col_force]
        Rf = db.concat(ru_temp, ignore_index=True)
        Rf = Rf.groupby(header)[col_force].sum()
        Rf.reset_index(inplace=True)
        #
        # ------------------------------
        # end process
        # ------------------------------
        #
        #uptime = time.time() - start_time
        #print(f"** {{F}} = [Kt] {{U}} Solution: {uptime:1.4e} sec")        
        #1 / 0
        return Un, Rf        
    #
    #
    def tb_removed(self):
        pass
        #df =  DBframework()
        #
        #header = ['node_name', *col_disp]
        #header = ['load_name', 'load_id', 'load_level',
        #          'system', 'load_title', 'mesh_name',
        #          'result_name']          
        #Utemp = []
        #for key in load_keys:
            #Fs_item = Fn_zeros.copy()
            #Ke_item = Ke.copy()
            #
            #Ke_item = Ke.copy()
            #Fs_item = Fn_zeros.copy()
            #
            # --------------------------------------
            # Node Displacement Section
            # --------------------------------------
            #Dn_item = Dn_zeros.copy()
            #Dn_set = Dn_grp.get_group((key[-1], ))
            #Dn_set = Dn_set[head_disp].set_index(['node_name'])
            #Dn_item.loc[Dn_set.index] = Dn_set.astype('float32')
            #
            # --------------------------------------
            # Nodal Load Section
            # --------------------------------------
            #try:
            #    Fn_set = Fn_grp.get_group(key)
            #    Fn_set = Fn_set[head_force].set_index(['node_name'])
            #    Fs_item.loc[Fn_set.index] = Fn_set.astype('float32')
            #except (KeyError, AttributeError):
            #    pass             
            #
            # Load Section
            #try:
            #    Fn_item = Fn_grp.get_group(key)
            #    Fn_set = Fn_item[head_force].set_index(['node_name'])
            #    Fs_item.loc[Fn_set.index] = Fn_set.astype('float64')
            #except (KeyError, AttributeError):
            #    pass
            #
            # Displacement Section
            #try:
            #    Dn_item = Dn_grp.get_group(key)
            #    FDs = Dn_item[head_disp].set_index(['node_name'])
            #    #
            #    # calculate force {F} = {D} x Kmod
            #    FDs[col_disp] = FDs[col_disp].mul(head_stiff)
            #    FDs.rename(columns=col_rename, inplace=True)
            #    #
            #    Fs_item.loc[FDs.index] = Fs_item.loc[FDs.index].add(FDs[col_force])
            #    #
            #    Ke_item = self._update_K_FDn(Ke_item=Ke_item,
            #                                 FDs=FDs,
            #                                 Kp=Kp, ndof=ndof,
            #                                 col_force=col_force)
            #except (KeyError, AttributeError):
            #    pass
            #
            # get load matrix as vector
            #Fs_item = Fs_item.stack()
            #Fs_item = Fs_item.loc[Fn_bool]
            #
            #f = Fn_zeros.copy()
            #un = Dn_zeros.copy()
            #f = Fs_item * 0.0
            # Solve displacements {U} = [k]{Fu}
            #du = self._static._KeFu(Ke_item, f, solver, jbc)
            #du = [[item, *du[x]]
            #       for x, item in enumerate(jbc.index)]
            #du = list2df(data=du, header=header)
            #Kt = Ke_item
            #
            #du_temp = []
            #ru_temp = []
            #R = None
            #for tk in time:
            ##for tk in range(10):
            #    f_prev = f.copy()
            #    f = Fs_item * tk  # 0.10
            #    f_delta = f - f_prev   # force increment
            #    #
            #    #
            #    un_prev = un.copy()
            #    un = Dn_item * tk
            #    un_delta = un - un_prev
            #    #
            #    # ------------------------------
            #    # Assembly tangent matrix [Kt]
            #    #
            #    #Kt, R = self._mesh.Kt_R(Dn=du, Rb=R)
            #    #
            #    #Kg = self._mesh.Kg(Dn=du)
            #    #Kt = Ke_item + Kg
            #    #
            #    # Solve displacements {U} = [k]{Fu}
            #    #du, ru = self._DMCiRC(Kt, df, jbc)
            #    #
            #    # --------------------------------------
            #    # Solve displacements {U} = [k]{Fu}
            #    # --------------------------------------
            #    du, ru = self._DMCiRC(Kt, f_delta, un_delta, jbc)                
            #    #
            #    # pack results for df's dumping
            #    #u.extend([[*key, result_name, item, *du[x]]
            #    #          for x, item in enumerate(jbc.index)])
            #    #
            #    #du, Ke_item = self._step_du(f, Ke_item, jbc, solver, header)
            #    #
            #    # ------------------------------
            #    # Solve displacements {U} = [k]{Fu}
            #    #du = self._KeFu(Ke_item, df, solver, jbc)
            #    #du = [[item, *du[x]]
            #    #       for x, item in enumerate(jbc.index)]
            #    #du = list2df(data=du, header=header)
            #    #
            #    # ------------------------------
            #    Kt, R = self._mesh.Kt_R(Dn=du, Rb=R)                
            #    #
            #    # --------------------------------------
            #    # Add load reference to displacement
            #    # --------------------------------------
            #    #Utemp.extend([[*key, self._result_name, item, *Us[x]]
            #    #              for x, item in enumerate(jbc.index)])
            #    du[header] = [*key, self._result_name]
            #    du_temp.append(du)
            #    #
            #    # --------------------------------------
            #    # Add load reference to reactions
            #    # --------------------------------------
            #    ru[header] = [*key, self._result_name]
            #    ru_temp.append(ru)                
            #    #
        #
        # --------------------------------------
        #
        #header = ['load_name', 'load_id', 'load_level',
        #          'system', 'load_title', 'mesh_name',
        #          'result_name','node_name', *col_disp]
        #Un = list2df(data=u, header=header)
        #header = ['load_name', 'load_id', 'load_level',
        #          'system', 'load_title', 'mesh_name',
        #          'result_name','node_name']
        header = ['load_name', 'load_id', 'load_level',
                  'system', 'load_title', 'mesh_name',
                  'result_name','node_name']
        #
        db = DBframework()     
        #1 / 0
        #header_df = [*header, 'node_name']
        Un = db.concat(dn_temp, ignore_index=True)
        Un = Un.groupby(header)[col_disp].sum()
        Un.reset_index(inplace=True)
        #
        #header_df = [*header, 'node_name', *col_force]
        Rf = db.concat(ru_temp, ignore_index=True)
        Rf = Rf.groupby(header)[col_force].sum()
        Rf.reset_index(inplace=True)
        #
        # ------------------------------
        # end process
        # ------------------------------
        #
        #uptime = time.time() - start_time
        #print(f"** {{F}} = [Kt] {{U}} Solution: {uptime:1.4e} sec")        
        #1 / 0
        return Un, Rf
    #
    def _step_du(self, f, Ke_item, jbc, solver, header):
        """ """
        # Solve displacements {U} = [k]{Fu}
        du0 = self._DMCiRC(Ke_item, f, solver, jbc)
        du0 = [[item, *du0[x]]
               for x, item in enumerate(jbc.index)]
        du0 = list2df(data=du0, header=header)
        # ------------------------------
        # Assembly tangent matrix [Kt]
        Kt = self._mesh.Kt(Dn=du0)
        #
        # Solve displacements {U} = [k]{Fu}
        du = self._DMCiRC(Kt, f, solver, jbc)
        du = [[item, *du[x]]
               for x, item in enumerate(jbc.index)]
        du = list2df(data=du, header=header)
        #
        return du, Kt
    #
    # -----------------------------------------------------------
    #
    def _get_Fb(self, beams,
                #Un0: DBframework.DataFrame,
                Un: DBframework.DataFrame,
                Pdelta: bool) -> DBframework.DataFrame:
        """
        Convert beam end-node disp to force [F = Kd] in global system

        Return:
        Fb_local : ?
        Fb_global : ?
        """
        # Setup
        #dim = 3
        plane2D: bool = False
        Kglobal: str = 'Ke'
        if Pdelta:
            Kglobal: str = 'Kt'
        #
        Un_grp = Un.groupby(['load_name', 'mesh_name', 'load_level'])
        head_disp = ['node_name', 'x', 'y', 'z', 'rx', 'ry', 'rz']
        # head_disp = ['node_name', *self._mesh._plane.hdisp]
        #
        Fb_temp: list = []
        for key, Un_item in Un_grp:
            Un_set = Un_item[head_disp].set_index('node_name')
            #
            #rot = {}
            #for item in Un_set.itertuples():
            #    rot[item.Index] = R_from_drot(item[4:]) @ np.eye(3)
            #
            for name, beam in beams.items():
                nodes = beam.connectivity
                #Tb = beam.T3D()
                #
                # ---------------------------------------------
                # beam end-nodes global displacement {U}
                Un_global = np.concatenate((Un_set.loc[nodes[0]],
                                            Un_set.loc[nodes[1]]),
                                           axis=None)
                #
                #                
                #
                # ---------------------------------------------
                #
                #
                #T0 = Tb[:dim, :dim]
                #t0 = T0[2, :]
                #
                #e1,e2,e3 = beam.unit_vector
                #
                #
                # assigns .tmat (new e2 --> change), .L, .e, .psi from methods
                #
                #
                # ---------------------------------------------
                # update geometry
                #
                beam_mod = beam.deformed(du=Un_global)
                Tb = beam_mod.T3D()
                #
                # update base vectors of element from rotations of nodes
                #
                #R1 = rot[nodes[0]]
                #R2 = rot[nodes[1]]
                #e3_temp = R1 @ t0 + R2 @ t0
                #
                #e2 = np.cross(e3_temp, e1)/np.linalg.norm(np.cross(e3_temp, e1))
                #
                #T0 = transform_unit(e1, e2=e2)
                #
                #Tb = blkdiag(T0, n=4)
                #
                # update K
                #1 / 0
                #
                # --------------------------------------------
                # get beam [K] global
                Kb = getattr(beam_mod, Kglobal)(plane2D, Un_global)
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
                # Fb_local2 = Tb @ Fb_global
                # ---------------------------------------------
                # Build list with results
                # [load_name, mesh_name, load_level, system, element_name,
                #  node_name, Un_local, FU_local, FU_global]
                #
                Fb_temp.append([*key, name, nodes,
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
    #
    # -----------------------------------------------------------
    #
    def solve_beam_end(self, elements,
                       Un:DBframework, #Rn:DBframework,
                       jbc, Pdelta: bool)-> DBframework: #jbc:DBframework, 
        """
        Convert beam end-node disp to force [F = Kd] in global system

        Return:
        Fb_local : ?
        Fb_global : ?
        """
        #
        #jbc = self._mesh.jbc(plane2D=False)
        #D1, D2 = self._mesh._D_partition(jbc)
        #1/0
        # Setup
        plane2D:bool = False
        #plane2D = self._mesh._plane.plane2D
        #Kglobal:str = 'Ke'
        #if Pdelta:
        #    Kglobal:str = 'Kt'
        #
        #col_grp = ['load_name', 'load_id', 'mesh_name', 'load_level']
        #Un_grp = Un.groupby(col_grp)
        #Rn_grp = Rn.groupby(col_grp)
        #head_disp = ['node_name', 'x', 'y', 'z', 'rx', 'ry', 'rz']
        ##head_disp = ['node_name', *self._mesh._plane.hdisp]
        #head_force = ['node_name', 'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']
        #head_force = ['node_name', *self._mesh._plane.hforce]
        #
        # ---------------------------------------------
        #
        nn = len(jbc)
        ndof = 6
        Ka = dok_matrix((nn * ndof, nn * ndof), dtype=np.float32)
        #Ka = np.zeros((nn * ndof, nn * ndof), dtype=np.float32)
        #beam_temp = []
        Fb_temp: list = []
        #for key, Un_item in Un_grp:
        #    Un_set = Un_item[head_disp].set_index('node_name')
        #    Rn_item = Rn_grp.get_group(key)
        #    Rn_set = Rn_item[head_force].set_index('node_name')
        #    #
        #
        #Un_set = Un[head_disp].set_index(['node_name'])
        #Rn_set = Rn[head_force].set_index('node_name')
        #
        for e_name, element in elements.items():
            nodes = element.connectivity
            #Tbx = element.T(plane2D=plane2D)
            # ---------------------------------------------
            # beam end-nodes global displacement {U}
            #
            #jbc_set = jbc.loc[nodes].stack(future_stack=True)
            #D1 = (jbc_set == 0).tolist()
            #D2 = (jbc_set != 0).tolist()
            #
            Un_global = Un.loc[nodes].stack(future_stack=True)
            #
            #
            #Rn_global = Rn_set.loc[nodes].stack(future_stack=True)
            #
            #L, e1, rot, Tb = self.get_rot(du=Un_global.to_list(),
            #                              beam=element)            
            #
            #
            #Kb = getattr(element, 'Ke')(plane2D)
            if Pdelta:
                # --------------------------------------------
                beam_def = element.deformed2(Un=Un_global.to_list())
                Tb = beam_def.T(plane2D)                
                #Fb = Kb @ Un_local
                #Kb = getattr(element, 'Kt')(plane2D, Fb)
                Kb = beam_def.Kt(plane2D)
            else:
                Tb = element.T(plane2D=plane2D)
                # get beam [K] global
                Kb = element.Ke(plane2D)
            #
            # ---------------------------------------------
            # Convert global end-node disp in beam's local system
            Un_local = Tb @ Un_global
            #
            #Un_localx = Tbx @ Un_global            
            #    #Kg = getattr(element, 'Kg')(plane2D, Fb)
            #    #Kb += Kg
            #Kb = getattr(element, Kglobal)(plane2D, Un_global)
            # ---------------------------------------------
            # convert beam end-node disp to force in global system
            # {F} = [K]{U}
            Fb_global = Kb @ Un_global
            #
            #Fb_res = Kb @ Un_global
            #Fb_global = Rn_global.copy()
            #Fb_global[:] = float(0.0)
            #Fb_global.iloc[D1] += Rn_global.iloc[D1]
            #Fb_global.iloc[D2] += (Fb_res[D2]).astype('float32')
            #
            # ---------------------------------------------
            # Calculate beam end-nodes force in local system
            Fb_local = Tb @ Fb_global
            # ---------------------------------------------
            # Build list with results
            # [load_name, mesh_name, load_level, system, element_name,
            #  node_name, Un_local, FU_local, FU_global]
            #
            #L, e1, rot = self.get_rot(du=Un_global.to_list(),
            #                          beam=element)
            #
            #
            Fb_temp.append([e_name, #L, e1, rot,
                            Un_local, Fb_local])
            #
            #
            assembly(Kb, Ka, ndof, element.DoF)            
        #
        # ---------------------------------------------
        #
        header = [#'load_name', 'load_id',
                  #'mesh_name', 'load_level',
                  'element_name',
                  ##'L', 'uv', 'R', 
                  'Un_local', 'Fb_local']
        Fb_temp = list2df(data=Fb_temp, header=header)
        #
        print('--> beam ends')
        #Ka =  coo_matrix(Ka)
        return Ka, Fb_temp
    #
    #
    def solve_beam_endXX(self, elements,
                       Un:DBframework, Rn:DBframework,
                       Pdelta: bool)-> DBframework: #jbc:DBframework, 
        """
        Convert beam end-node disp to force [F = Kd] in global system

        Return:
        Fb_local : ?
        Fb_global : ?
        """
        #
        jbc = self._mesh.jbc(plane2D=False)
        D1, D2 = self._mesh._D_partition(jbc)
        #1/0
        # Setup
        plane2D:bool = False
        #plane2D = self._mesh._plane.plane2D
        Kglobal:str = 'Ke'
        if Pdelta:
            Kglobal:str = 'Kt'
        #
        col_grp = ['load_name', 'load_id', 'mesh_name', 'load_level']
        Un_grp = Un.groupby(col_grp)
        Rn_grp = Rn.groupby(col_grp)
        head_disp = ['node_name', 'x', 'y', 'z', 'rx', 'ry', 'rz']
        #head_disp = ['node_name', *self._mesh._plane.hdisp]
        head_force = ['node_name', 'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']
        #head_force = ['node_name', *self._mesh._plane.hforce]
        #
        # ---------------------------------------------
        #
        Fb_temp: list = []
        for key, Un_item in Un_grp:
            Un_set = Un_item[head_disp].set_index('node_name')
            Rn_item = Rn_grp.get_group(key)
            Rn_set = Rn_item[head_force].set_index('node_name')
            #
            for e_name, element in elements.items():
                nodes = element.connectivity
                Tb = element.T(plane2D=plane2D)
                # ---------------------------------------------
                # beam end-nodes global displacement {U}
                #
                jbc_set = jbc.loc[nodes].stack(future_stack=True)
                D1 = (jbc_set == 0).tolist()
                D2 = (jbc_set != 0).tolist()
                #
                Un_global = Un_set.loc[nodes].stack(future_stack=True)
                #
                Rn_global = Rn_set.loc[nodes].stack(future_stack=True)
                #
                # --------------------------------------------
                # get beam [K] global
                Kb = getattr(element, Kglobal)(plane2D, Un_global)
                # ---------------------------------------------
                # convert beam end-node disp to force in global system
                # {F} = [K]{U}
                Fb_res = Kb @ Un_global
                Fb_global = Rn_global.copy()
                Fb_global[:] = float(0.0)
                Fb_global.iloc[D1] += Rn_global.iloc[D1]
                Fb_global.iloc[D2] += (Fb_res[D2]).astype('float32')
                #
                # ---------------------------------------------
                # Convert global end-node disp in beam's local system
                Un_local = Tb @ Un_global
                # ---------------------------------------------
                # Calculate beam end-nodes force in local system
                Fb_local = Tb @ Fb_global
                # ---------------------------------------------
                # Build list with results
                # [load_name, mesh_name, load_level, system, element_name,
                #  node_name, Un_local, FU_local, FU_global]
                #
                L, e1, rot = self.get_rot(du=Un_global.to_list(),
                                          nodes=element.nodes)
                #
                Fb_temp.append([*key, e_name, L, e1, rot,
                                Un_local, Fb_local, Fb_global])
        #
        # ---------------------------------------------
        #
        header = ['load_name', 'load_id',
                  'mesh_name', 'load_level',
                  'element_name',
                  'L', 'e1', 'rotation', 
                  'Un_local', 'Fb_local', 'Fb_global']
        Fb_temp = list2df(data=Fb_temp, header=header)
        #Fb_temp['system'] = 'global'
        return Fb_temp
    #    
    #
    def get_rot(self, du, beam):
        """ """
        nodes = beam.nodes
        #dim = 3
        coord1_du = du[:3] # x,y,z
        coord2_du = du[6:9]
        #
        rot1 = du[3:6] # rx,ry,rz
        rot2 = du[9:]
        #
        node1, node2 = nodes
        #
        coord1 = list(map(add, node1[:3], coord1_du))
        coord2 = list(map(add, node2[:3], coord2_du))
        L = dist(coord1, coord2)
        e1 = list(map(add, coord1, coord2))
        e1 = [item / L for item in e1]
        #print(L)
        #1 / 0
        # -------------------------------------------------
        #
        dim = 3
        #
        Tb = beam.T3D()
        T0 = Tb[:dim, :dim]
        t0 = T0[2, :]
        #
        # -------------------------------------------------
        #
        R = [np.eye(3), np.eye(3)]
        #
        # update base vectors of element from rotations of nodes
        R1 = R_from_drot(rot1) @ R[0]
        #R1 = R1 @ R1
        R2 = R_from_drot(rot2) @ R[1]
        #R2 = R2 @ R2
        #
        #
        e3_temp = R1 @ t0 + R2 @ t0
        e2 = np.cross(e3_temp, e1) / np.linalg.norm(np.cross(e3_temp, e1))
        uv = transform_unit(e1, e2=e2)        
        #
        # -------------------------------------------------
        # Return deformed beam item class
        beam_def = BeamItemDeformed(beam_name=beam.name,
                                    material=beam.material,
                                    section=beam.section,
                                    #du=Fb_local, 
                                    L=L, DoF=beam.DoF, 
                                    unit_vector=uv,
                                    R = [R1, R2])       
        #
        Tb = beam_def.T(plane2D=False)
        #
        return L, uv, [R1, R2], Tb
        
#
# ---------------------------------------------
#
def assembly(keg, Ka: np.array,
             ndof:int, DoF):
    """
    element : Class
    Ka : The stiffness matrix's container
    ndof : node's degree of freedom
    matrix_type: Ke, Km, Kg, Kt
    plane2D: True/False
    Un: node global displacement

    Returns:
        The update stiffness matrix's container [Ka]
    """
    # TODO : check applicable to all element type
    #keg = getattr(element, matrix_type)(plane2D, Un)
    idof, jdof = DoF
    # node and corresponding dof (start, end)
    niqi, niqj = idof*ndof, idof*ndof + ndof
    njqi, njqj = jdof*ndof, jdof*ndof + ndof
    # assemble global stiffness matrix, quadrant 1 to 4
    Ka[niqi:niqj, niqi:niqj] += keg[:ndof, :ndof]             # 2nd
    Ka[niqi:niqj, njqi:njqj] += keg[:ndof, ndof:2*ndof]       # 1st
    Ka[njqi:njqj, niqi:niqj] += keg[ndof:2*ndof, :ndof]       # 3rd
    Ka[njqi:njqj, njqi:njqj] += keg[ndof:2*ndof, ndof:2*ndof] # 4th
    #return Ka
#
#
#
def transform_unitXX(e1:list, e2:list|None=None, e3:list|None=None,
                   warnings:bool=False):
    '''
    Establish transformation matrix from e1 and temporary e2 or e3 vectors.

    Arguments
    -----------
    e1 : unit vector describing element longitudinal direction
    e2 : temporary unit vector describing a chosen vector that's perpendicular to the longitudinal direction (approximate y-direction)
    e3 : temporary unit vector describing a chosen vector that's perpendicular to the longitudinal direction (approximate z-direction)
        if both e2 and e3 are different from None, e2 is used (e3 disregarded)

    Returns:
        T : transformation matrix
    '''

    e1 = np.array(e1).flatten()

    if (e2 is not None) and (e3 is not None):
        e3 = None

    if e2 is not None:
        e2 = np.array(e2).flatten()
        e3_tmp = np.cross(e1, e2)  # Direction of the third unit vector
        if np.all(e3 == 0):
            if warnings:
                print('Warning: e1 and e2 identical. Check orientations carefully!')

            e2 = e2 + 0.1
            e3 = np.cross(e1, e2)

        e2 = np.cross(e3_tmp, e1)  # Direction of the second unit vector

        e1 = e1 / np.linalg.norm(e1)  # Normalize the direction vectors to become unit vectors
        e2 = e2 / np.linalg.norm(e2)
        e3 = np.cross(e1, e2)

    elif e3 is not None:
        e3 = np.array(e3).flatten()
        e2_tmp = np.cross(e3, e1)  # Direction of the third unit vector
        e3 = np.cross(e1, e2_tmp)
        e1 = e1 / np.linalg.norm(e1)  # Normalize the direction vectors to become unit vectors
        e3 = e3 / np.linalg.norm(e3)
        e2 = np.cross(e3, e1)

    if e2 is None and e3 is None:
        raise ValueError('Specify either e2 or e3')

    T = np.vstack([e1, e2, e3])

    return T
#
# ---------------------------------------------
#
def list2df(data:list[list], header:list):
    """ dataframe"""
    db = DBframework()
    dftemp = db.DataFrame(data=data, columns=header, index=None)
    return dftemp
#
# ---------------------------------------------
#
def is_converged(values, tols, scaling=None):
    """
    Check whether multiple values are below specified tolerances.  (value/scaling)
    
    Arguments
    -----------------
    values : double
        list of values to check
    tols : double
        corresponding list of tolerances to compare values to
    scaling : double, optional
        corresponding list of scaling of values


    Returns
    -----------------
    ok : boolean
        converged or not?

    Notes
    --------------------
    If entry in tols is None, the corresponding value is assumed to pass tolerance criterion. 
    If entry in tols is different than None, the value pass if value <= tol * scaling.
            

    """

    if scaling is None:
        scaling = np.ones([len(tols)])
    
    for ix, value in enumerate(values):
        if (tols[ix] is not None) and (value>(tols[ix]*scaling[ix])):
            return False
        
    return True
#
#