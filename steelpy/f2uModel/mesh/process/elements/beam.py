# 
# Copyright (c) 2019-2023 steelpy
# 

# Python stdlib imports
from __future__ import annotations
from array import array
from dataclasses import dataclass
from collections.abc import Mapping
#from typing import NamedTuple

# package imports
from steelpy.process.geometry.L3D import (DistancePointLine3D,
                                          LineLineIntersect3D)
from steelpy.f2uModel.mesh.process.elements.bstiffness import (beam3D_Klocal, trans_3d_beam,
                                                               Rmatrix, Rmatrix_new,
                                                               trans3Dbeam, beam3D_K)
from steelpy.f2uModel.mesh.process.elements.bstiffness2D import (beam2D_K, beam2D_Geom, Rmatrix2D)
#
import numpy as np
from numpy.linalg import inv

#
class BeamBasic(Mapping):

    
    def __init__(self) -> None:
        """
        Beam element 
        """
        #
        #
        self._labels: list[int|str] = []
        self._number: array = array('i', [])
        self._title: list[int|str] = []
        #
    #  
    #
    def __contains__(self, value) -> bool:
        return value in self._labels
    
    def __iter__(self):
        """
        """
        return iter(self._labels)    

    def __len__(self) -> float:
        return len(self._labels)
    #
    def __str__(self) -> str:
        """ """
        print('beam basic')
        1/0
#
#
#
@dataclass
class BeamItemBasic:
    """
    """
    #__slots__ = ['name', 'type']

    def __init__(self, element_name: int|str) -> None:
        """
        """
        self.name: int|str = element_name
        self._releases: list[bool] = [False] * 12
        #self._intersect = LineLineIntersect3D()
    #
    #
    #@property
    def T(self, m2D:bool = False):
        """
        Returns the transformation matrix for the member
        """
        if m2D:
            # Element length and orientation
            node1,  node2 = self.nodes
            Tlg = Rmatrix2D(node1, node2)
            return Tlg
        else:
            #if self.type in ['beam', 'truss']:
            return Rmatrix(*self.unit_vector, self.beta)
            #else:
            #    raise IOError("no yet included")    
    #
    #@property
    def K(self, m2D:bool = False):
        """
        m2d : Matrix 2D (False default --> 3D) 
        
        Return the stiffness matrix in global coordinates
        """
        if m2D:
            Kl = self.k2D()
            Tlg =  self.T(m2D=m2D)
            Kg = (np.transpose(Tlg).dot(Kl)).dot(Tlg)
            return Kg
        else:
            k = self.k3D()
            #
            #
            #self.beta = 30
            #disb = self.R
            #node1 = self._cls._f2u_nodes[self._cls._connectivity[self.index][0]]
            #node2 = self._cls._f2u_nodes[self._cls._connectivity[self.index][-1]]
            #dirc = Rmatrix_new(node1[:3], node2[:3], self.beta)
            #xxx = trans3Dbeam(K, node1[:3], node2[:3])
            #return trans3Dbeam(K, node1[:3], node2[:3])
            #return trans_3d_beam(K, dirc)
            return trans_3d_beam(k, self.T())
            #return self.trans3d(K)    
    #
    def _partition(self, unp_matrix):
        """
        Partitions a matrix into sub-matrices based on unreleased and released degree of freedom indices.
        """
        unp_matrix =  np.array(unp_matrix)
        # Create auxiliary lists of released/unreleased DOFs
        R1_indices, R2_indices = self._aux_list()
        # Partition the matrix by slicing
        if unp_matrix.shape[1] == 1:
            m1 = unp_matrix[R1_indices, :]
            m2 = unp_matrix[R2_indices, :]
            return m1, m2
        else:
            m11 = unp_matrix[R1_indices, :][:, R1_indices]
            m12 = unp_matrix[R1_indices, :][:, R2_indices]
            m21 = unp_matrix[R2_indices, :][:, R1_indices]
            m22 = unp_matrix[R2_indices, :][:, R2_indices]
            return  m11, m12, m21, m22
    #
    def _aux_list(self):
        """
        Builds lists of unreleased and released degree of freedom indices for the member.

        Returns
        -------
        R1_indices : list
            A list of the indices for the unreleased DOFs
        R2_indices : list
            A list of the indices for the released DOFs
        """
        R1i = [i for i, rel in enumerate(self._releases)
                      if not rel]
        R2i = [i for i, rel in enumerate(self._releases)
               if rel]
        return R1i, R2i
    #
    def _k_unc(self, k_unc):
        """Returns the uncondensed local stiffness matrix for the member"""
        # Partition the local stiffness matrix as 4 submatrices in
        # preparation for static condensation
        k11, k12, k21, k22 = self._partition(k_unc)
        #
        # Calculate the condensed local stiffness matrix
        k_condensed = np.subtract(k11, np.matmul(np.matmul(k12, inv(k22)), k21))
        #
        # Expand the condensed local stiffness matrix
        for i, DOF in enumerate(self._releases):
            if DOF:
                k_condensed = np.insert(k_condensed, i, 0, axis = 0)
                k_condensed = np.insert(k_condensed, i, 0, axis = 1)
        # Return the local stiffness matrix, with end releases applied
        return k_condensed
    #
    #
    def k3D(self):
        """Returns the condensed (and expanded) local stiffness matrix for the member"""
        # get section properties 
        section = self.section
        #section = self._cls._f2u_sections[section]
        section = section.properties()
        # get material properties
        material = self.material
        #material = self._cls._f2u_materials[material]        
        #
        #
        # solve K matrix
        #Kb = beam3D_Klocal(self.L,
        #                    self.area, self.J,
        #                    self.Iy, self.Iz,
        #                    self.E, self.G,
        #                    self.area, self.area)
        #        
        kb = beam3D_K(self.L,
                      section.area, section.J,
                      section.Iy, section.Iz,
                      material.E, material.G)
        #
        k_cond = self._k_unc(kb)
        return k_cond
        #return K
    #
    #
    def k2D(self, beta:int=0):
        """
        """
        # get section properties 
        section = self.section
        #section = self._cls._f2u_sections[section]
        section = section.properties()
        # get material properties
        material = self.material
        #
        k_cond =beam2D_K(L=self.L, 
                         A=section.area, I=section.Iz, 
                         E=material.E, beta=beta)
        return k_cond
    #
    # Mass
    #@property
    def Km(self):
        """
        """
        #FIXME: material 
        #section = self._cls._f2u_sections[self.section].properties
        #material = self._cls._f2u_materials[self.material]
        #
        if self.type in ['beam', 'truss']:
            em = beam_mass(self.length,
                           section, material,
                           ilump=2)
            return trans_3d_beam(em, self.R)
        else:
            raise IOError("no yet included")     
    #  
    #
    # Geometry
    #
    def intersect_point(self, point:list):
        """line intersection """
        p1, p2 = self.nodes
        point_line = DistancePointLine3D(p1[:3], p2[:3])
        return point_line.is_on_segment(point[:3])
    #