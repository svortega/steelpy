# 
# Copyright (c) 2009 steelpy
# 
# Python stdlib imports
from __future__ import annotations
#from array import array
from dataclasses import dataclass
from collections.abc import Mapping
#from typing import NamedTuple
import math

# package imports
from steelpy.ufo.mesh.process.bstiffness import beam3D_B3D2, B3D2_Ke
from steelpy.ufo.mesh.process.brotation import Tmatrix
from steelpy.ufo.mesh.process.bgeometry import B3D2_Kt, beam_KG
from steelpy.ufo.mesh.process.bmass import beam_mass
from steelpy.utils.geometry.L3D import DistancePointLine3D #, LineLineIntersect3D
from steelpy.utils.math.operations import remove_column_row
#
import numpy as np
from numpy.linalg import inv
#
#
#
#
class BeamBasic(Mapping):
    __slots__ = ['_name',]
    
    def __init__(self, name: int|str) -> None:
        """
        Beam element 
        """
        self._name = name
    #
    # ------------------------------------------------
    #
    def __contains__(self, value) -> bool:
        return value in self._labels
    
    def __iter__(self):
        """
        """
        labels = list(dict.fromkeys(self._labels))
        return iter(labels)

    def __len__(self) -> float:
        labels = list(dict.fromkeys(self._labels))
        return len(labels)
    #
    def __str__(self) -> str:
        """ """
        lenght = ' m'
        space = " "
        #
        output = "\n"
        output += "{:}\n".format(80*"_")
        output += "\n"
        output += f"{30*space}BEAM ELEMENTS\n"
        output += "\n"
        output += "\n"
        output += (f"Beam {7*space}Node1    Node2 {4*space}Material {5*space}Section")
        output += (f" {4*space}Beta {3*space}Len[{lenght}] {2*space}Title")
        output += "\n"
        output += "{:}\n".format(80*".")
        output += "\n"
        for beam_name in self._labels:
            beam = self.__getitem__(beam_name)
            output += beam.__str__()
        #
        return output
    #
    # ------------------------------------------------
    #
#
#
#
@dataclass
class BeamItemBasic:
    """
    """
    __slots__ = ['name', '_releases', 'type']

    def __init__(self, beam_name: int|str) -> None:
        """
        """
        self.name = beam_name
        self._releases: list[bool] = [False] * 12
        self.type: str = "beam"
    #
    # ------------------------------------------------
    # Transformation
    #
    #@property
    def T(self, plane2D: bool = False):
        """
        Returns Beam 2D/3D transformation matrix
        """
        Tlg = self.T3D()
        if plane2D:
            Tlg = self._M2D(Tlg)
        return Tlg
    #
    def T3D(self):
        """ """
        #nodei, nodej = self.nodes
        #return Rmatrix(*self.unit_vector, self.beta)
        #return Rmatrix2(nodei, nodej)
        unitvec = np.array(self.unit_vector)
        return Tmatrix(dirCos=unitvec)
    #
    # ------------------------------------------------
    # Stiffness
    #
    #@property
    #def Ke_global(self, plane2D: bool):
    #    """
    #    Return Beam stiffness matrix in global coordinates
    #    """
    #    return self.Ke(plane2D)
    #
    def Ke(self, plane2D: bool, item: None = None):
        """ Return Beam 2D/3D stiffness matrix in global coordinates """
        Tlg = self.T(plane2D)
        Kl = self.Ke_local(plane2D)
        return Tlg.T @ Kl @ Tlg        
    #
    #@property
    def Ke_local(self, plane2D: bool, item: None = None):
        """Return Beam 2D/3D stiffness matrix in local coordinates """
        ke = self.Ke3D()
        if plane2D:
            ke = self._M2D(ke)
        return ke
    #
    def Ke3D(self, item: None = None):
        """Returns the condensed (and expanded) local stiffness matrix for 3D beam"""
        # get material properties
        material = self.material        
        # get section properties 
        section = self.section.properties(poisson=material.poisson)     
        #
        # solve K matrix
        #
        kb = B3D2_Ke(Le=self.L,
                    Ax=section.area,
                    Asy=section.Asz,
                    Asz=section.Asy,
                    Jx=section.J,
                    Iy=section.Iz, Iz=section.Iy,
                    Emod=material.E, Gmod=material.G,
                    shear=True)
        #        
        # TODO : check if this element is reliable
        #kb2 = beam3D_B3D2(Le=self.L,
        #                 Ax=section.area,
        #                 Asy=section.Asy,
        #                 Asz=section.Asz,
        #                 Jx=section.J,
        #                 Iy=section.Iy, Iz=section.Iz,
        #                 Emod=material.E, Gmod=material.G,
        #                 shear=True)
        #
        k_cond = self._k_unc(kb)
        #1 / 0
        return k_cond   
    #
    # ------------------------------------------------
    # Geometry
    #
    def Kg(self, plane2D: bool, Un: list):
        """
        Return Beam geometrical 2D/3D stiffness matrix in global coordinates
        """
        Tb =  self.T
        Kg = self.Kg_local(Un=Un)
        return Tb.T @ Kg @ Tb
    #
    def Kg_local(self, plane2D: bool, Un:list):
        """
        Return Beam geometrical 2D/3D stiffness matrix in local coordinates
        """
        kg = self.Kg3D(Un=Un)
        if self._plane.plane2D:
            kg = self._M2D(kg)
        return kg
    #
    def Kg3D(self, Un: list):
        """
        Returns the condensed (and expanded) local geometrical stiffness matrix for 3D beam
        """
        material = self.material
        section = self.section.properties(poisson=material.poisson)
        #
        # ---------------------------------------------
        # convert global end-node disp in beam's local system
        nd_local = self.T @ Un # nd_global
        #
        #T = nd_local[6] - nd_local[0]
        #P = material.E * section.area / self.L * T
        #
        # ---------------------------------------------
        # convert beam end-node disp to force [F = Kd] in global system
        Fb = self.Ke_local @ nd_local
        #
        #kg = GeoStfBeamTimoshenko(Le=self.L,
        #                          Ax=section.area,
        #                          Asy=section.Asz,
        #                           Asz=section.Asy,
        #                           Jx=section.J,
        #                           Iy=section.Iz, Iz=section.Iy,
        #                           Emod=material.E, Gmod=material.G,
        #                           Fb=Fb,
        #                           shear=True,
        #                           order=4)
        #
        kg = beam_KG(Le=self.L,
                     Ax=section.area,
                     Asy=section.Asy,
                     Asz=section.Asz,
                     Jx=section.J,
                     Iy=section.Iy, Iz=section.Iz,
                     Emod=material.E, Gmod=material.G,
                     Fb=Fb,
                     shear=True)
        #
        k_cond = self._k_unc(kg)
        return k_cond
    #
    # ------------------------------------------------
    # Tangent
    #
    def Kt(self, plane2D: bool, Un: list):
        """ Return Beam tangent stiffness matrix in global coordinates"""
        Tb = self.T(plane2D)
        Kt = self.Kt_local(plane2D=plane2D, Un=Un)
        return Tb.T @ Kt @ Tb
    #
    def Kt_local(self, plane2D: bool, Un:list):
        """
        Un : node displacement global system

        Return :
        Kt : Beam geometrical stiffness matrix in local coordinates
        """
        # ---------------------------------------------
        # convert global end-node disp in beam's local system
        nd_local = self.T(plane2D) @ Un # nd_global
        # ---------------------------------------------
        # convert beam end-node disp to force [F = Kd] in local system
        Fb = self.Ke_local(plane2D) @ nd_local
        # ---------------------------------------------
        kt = self.Kt3D(Fb=Fb)
        if plane2D:
            kt = self._M2D(kt)
        return kt
    #
    def Kt3D(self, Fb: list):
        """
        Fb : Member force local system
        Returns :
        Kt : the condensed (and expanded) local tangent stiffness matrix for 3D beam
        """
        material = self.material
        section = self.section.properties(poisson=material.poisson)
        # ---------------------------------------------
        kg = B3D2_Kt(Le=self.L,
                     Ax=section.area,
                     Asy=section.Asz,
                     Asz=section.Asy,
                     Jx=section.J,
                     Iy=section.Iz, Iz=section.Iy,
                     Emod=material.E, Gmod=material.G,
                     Fb=Fb,
                     shear=True,
                     toler=0.01)
        #
        # ---------------------------------------------
        k_cond = self._k_unc(kg)
        return k_cond
    #
    # ------------------------------------------------
    # Mass
    #
    #@property
    #def Km_global(self, plane2D: bool):
    #    """
    #    Return Beam 2D/3D mass matrix in global coordinates
    #    """
    #    #Tlg =  self.T
    #    #Km = self.Km_local
    #    #return (np.transpose(Tlg).dot(Kl)).dot(Tlg)
    #    #return Tlg.T @ Km @ Tlg
    #    return self.Km(plane2D)
    #
    def Km(self, plane2D: bool, item: None = None):
        """
        Return Beam 2D/3D mass matrix in global coordinates
        """
        Tlg =  self.T(plane2D)
        Km = self.Km_local
        return Tlg.T @ Km @ Tlg        
    #
    #@property
    def Km_local(self, plane2D: bool):
        """
        Return Beam 2D/3D mass matrix in local coordinates
        """
        km = self.Km3D()
        if plane2D:
            km = self._M2D(km)
        return km
    #
    def Km3D(self, item: None = None):
        """ Returns the condensed (and expanded) local mass matrix for 3D beam"""
        # get section properties 
        material = self.material
        section = self.section.properties(poisson=material.poisson)        
        #
        em = beam_mass(self.length,
                       section, material,
                       ilump=2)
        #
        k_cond = self._k_unc(em)
        return k_cond        
    #
    # ------------------------------------------------
    # Matrix operations
    #
    def _partition(self, unp_matrix):
        """
        Partitions a matrix into sub-matrices based on
        unreleased and released degree of freedom indices.
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
            return m11, m12, m21, m22
    #
    def _aux_list(self):
        """
        Builds lists of unreleased and released degree of freedom
        indices for the member.

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
        """Returns the uncondensed local beam's stiffness matrix"""
        # Partition the local stiffness matrix as 4 submatrices in
        # preparation for static condensation
        k11, k12, k21, k22 = self._partition(k_unc)
        # Calculate the condensed local stiffness matrix
        k_condensed = np.subtract(k11, np.matmul(np.matmul(k12, inv(k22)), k21))
        # Expand the condensed local stiffness matrix
        for i, DOF in enumerate(self._releases):
            if DOF:
                k_condensed = np.insert(k_condensed, i, 0, axis = 0)
                k_condensed = np.insert(k_condensed, i, 0, axis = 1)
        # Return the local stiffness matrix, with end releases applied
        return k_condensed
    #
    def _M2D(self, ka, dof: list[str] = ['x', 'y', 'rz']):
        """Reduce matrix to 2D"""
        # removing z, Mx, My
        index_off = self.index_off(dof)
        for i in index_off:
            ka = remove_column_row(ka, i, i)
        return ka
    #
    # ------------------------------------------------
    # Operations
    #
    def intersect_point(self, point:list):
        """line intersection """
        p1, p2 = self.nodes
        point_line = DistancePointLine3D(p1[:3], p2[:3])
        return point_line.is_on_segment(point[:3])
    #
    #
    def find_coordinate(self, node_distance:float, node_end:int=0) -> tuple:
        """
        """
        node = self.nodes
        nodeNo3 = [0, 0, 0]
        #
        if math.isclose(node_distance, 0, rel_tol=0.01):
            nodeNo3[0] = node[node_end].x
            nodeNo3[1] = node[node_end].y
            nodeNo3[2] = node[node_end].z
        else:
            if node_end == 1:
                v1 = (node[0].x - node[1].x)
                v2 = (node[0].y - node[1].y)
                v3 = (node[0].z - node[1].z)
            else:
                v1 = (node[1].x - node[0].x)
                v2 = (node[1].y - node[0].y)
                v3 = (node[1].z - node[0].z)
            #
            norm = (v1**2 + v2**2 + v3**2)**0.50
            v1 /= norm
            v2 /= norm
            v3 /= norm
            #
            #L = math.dist(nodei[:3], nodej[:3])
            #uv = list(map(sub, nodej[:3], nodei[:3]))
            #l, m, n = [item / L for item in uv]
            #
            nodeNo3[0] = (node[node_end].x + v1 * node_distance)
            nodeNo3[1] = (node[node_end].y + v2 * node_distance)
            nodeNo3[2] = (node[node_end].z + v3 * node_distance)
        #
        return nodeNo3    
    #
    #
    #@property
    def index_off(self, dof: list[str]) -> list[int]:
        """return index"""
        index = {'x': 0, 'y': 1, 'z': 2, 'rx': 3, 'ry': 4, 'rz': 5} 
        idx3D = set(index.values())
        ndof = len(idx3D)
        #dofactive = set(self._index[item] for item in self.dof)
        getidx = [index[item] for item in dof]
        dofactive = set(getidx)
        indexes = list(idx3D - dofactive)
        indexes.extend([item + ndof for item in indexes])
        return list(reversed(indexes)) # [2, 3, 4, 8, 9, 10]
    #
    #