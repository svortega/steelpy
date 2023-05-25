#
# Copyright (c) 2019-2023 steelpy
#

# Python stdlib imports
from __future__ import annotations
from dataclasses import dataclass
from collections.abc import Mapping
#

# package imports
#
from steelpy.process.dataframe.main import DBframework
import numpy as np
#
#
#
# ----------------------------------------
#      Standard Sections Profiles
# ----------------------------------------
#
#
#
#
#
@dataclass(kw_only=True)
class ShapeBasic:
    name: str|int
    build: str = 'welded'
    #
    # -------------------------------------
    # Operations
    # -------------------------------------
    #
    def stress(self, actions=None, stress=None, df=None):
        """return cross section stress"""
        #print('-->')
        try:
            df.columns
            dfres = DBframework()
            # -------------------------------------
            try:
                stress_df = df[['node_end', 'tau_x', 'tau_y', 'tau_z',
                                'sigma_x', 'sigma_y', 'sigma_z']]
                stress = self._stress(stress=stress_df)
            except KeyError:
                actions_df = df[['load_title', 'node_end', 'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']]
                stress = self._stress(actions=actions_df)
                # -------------------------------------
                header = ['load_name', 'load_title', 
                          'element_name','node_end','stress_points', 'y', 'z']
                # -------------------------------------
                coord = stress.stress_point
                items = [[row.load_name, row.load_title,
                          row.element_name, row.node_end,
                          x+1, coord.y[x], coord.z[x]]
                         for x in range(len(coord.y))
                         for row in df.itertuples()]
                df_stress = dfres.DataFrame(data=items, columns=header, index=None)
                # -------------------------------------
                # axial stress
                df_stress['tau_x'] = np.array(stress.tau_x).flatten()
                # In Plane shear stress
                df_stress['tau_y'] = np.array(stress.tau_y).flatten()
                # Out Plane shear stress
                df_stress['tau_z'] = np.array(stress.tau_z).flatten()
                # torsion stress
                df_stress['sigma_x'] = np.array(stress.sigma_x).flatten()
                # In plane bending stress
                df_stress['sigma_y'] = np.array(stress.sigma_y).flatten()
                # Out plane bending stress
                df_stress['sigma_z'] = np.array(stress.sigma_z).flatten()
                # return dataframe
                return df_stress

        except AttributeError:
            if stress:
                print("stress")
            elif actions:
                print("actions")
            else:
                print('--> ??')
                1 / 0
    #
    #@property
    def properties(self):
        """
        --------------------------
        General Beam Element Data
        --------------------------
        
        Parameters  
        ----------
        area: Section area
        Zc  : Elastic neutral centre
        Yc  : Elastic neutral centre
        
        Iy  : Second moment of area about mayor axis
        Zy : Elastic modulus about mayor axis
        Sy : Plastic modulus about mayor axis
        Avy : Shear area mayor axis
        ry  : Radius of gyration about mayor Axis
        
        Iz  : Second moment of area about minor axis
        Zz : Elastic modulus about minor axis
        Sz : Plastic modulus about minor axis
        Avz : Shear area minor axis
        rz  : Radius of gyration about minor Axis
        
        SCz  : Shear centre about z axis
        SCy  : Shear centre about y axis
        
        Cw  : Warping constant
        """
        return self._properties()
    #
    #
    def push_property(self):
        """ """
        self.properties
    #    
    
#
#
#
class SectionBasic(Mapping):
    #__slots__ = ['_labels', '_number', '_title', '_type',
    #             '_tubular', '_solid', '_ibeam', '_box',
    #             '_channel', '_tee', '_angle', '_default']

    def __init__(self):
        """
        """
        self._default: str | None = None
        self._labels: list[str|int] = []
        self._number: list[int] = []
        self._title: list[str] = []
        self._type: list = []
    #
    #
    def __len__(self):
        return len(self._labels)

    def __iter__(self):
        return iter(self._labels)

    def __contains__(self, value):
        return value in self._labels
    #
    def __delitem__(self, shape_name: str | int) -> None:
        try:
            _nodes_empty = []
            _index = self._labels.index(shape_name)
            1/0
        except ValueError:
            raise KeyError(f'    *** warning section {shape_name} does not exist')
    #
    def properties(self):
        """
        """
        return self._properties()    
    #
    def get_number(self, start: int = 1):
        """
        """
        try:
            n = max(self._number) + 1
        except ValueError:
            n = start
        #
        while True:
            yield n
            n += 1
    #
    #
    def df(self):
        """ """
        db = DBframework()
        header = ['number', 'name', 'title', 'diameter', 'wall_thickness',
                  'height', 'web_thickness',
                  'top_flange_width', 'top_flange_thickness',
                  'bottom_flange_width', 'bottom_flange_thickness']
        1 / 0
        return db    