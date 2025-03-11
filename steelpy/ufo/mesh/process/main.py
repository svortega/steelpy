#
# Copyright (c) 2009 steelpy
#

# Python stdlib imports
from __future__ import annotations
from concurrent.futures import ThreadPoolExecutor
#import time
#from itertools import chain

# package imports
import numpy as np
#from scipy.sparse import coo_matrix
from scipy.sparse import dok_matrix

#
def Ke_matrix(elements,
              nodes,
              plane2D: bool):
              #sparse:bool = True):
    """
    elements: Class
    nodes: Class
    plane2D: True/False
    sparse : True/False

    Returns:
        The global Stiffness matrix [Ke]
    """
    matrix_type: str = "Ke"
    #start_time = time.time()
    Ka = assemble_KeKm_matrix(elements=elements,
                              nodes=nodes,
                              plane2D=plane2D,
                              matrix_type=matrix_type)
    #if sparse:
    #    Ka = coo_matrix(Ka)
    #    Ka = Ka.tolil()
    #
    #uptime = time.time() - start_time
    #print(f"** [{matrix_type}] assembly: {uptime:1.4e} sec")
    return Ka


#
def Km_matrix(elements,
              nodes,
              plane2D: bool):
              #sparse:bool = True):
    """
    elements: Class
    nodes: Class
    plane2D: True/False
    sparse : True/False

    Returns:
        The global Mass matrix [Km]
    """
    matrix_type: str = "Km"
    #start_time = time.time()
    Km = assemble_KeKm_matrix(elements=elements,
                              nodes=nodes,
                              plane2D=plane2D,
                              matrix_type=matrix_type)
    #if sparse:
    #    Km = coo_matrix(Km)
    #    Km = Km.tolil()
    #
    #uptime = time.time() - start_time
    #print(f"** [{matrix_type}] assembly: {uptime:1.4e} sec")
    return Km
#
# ------------------------------------------------------
#
def Kg_matrix(elements,
              nodes,
              #Un,
              plane2D: bool):
              #sparse:bool = True):
    """
    elements: Class
    nodes: Class
    Un: 
    plane2D: True/False
    sparse : True/False

    Returns:
        The Geometric stiffness matrix [Kg]
    """
    matrix_type: str = "Kg"
    #start_time = time.time()
    Kg = assemble_KtKg_matrix(elements=elements,
                              nodes=nodes,
                              plane2D=plane2D,
                              matrix_type=matrix_type)
                              #Un=Un)
    #if sparse:
    #    Kg = coo_matrix(Kg)
    #    Kg = Kg.tolil()
    #
    #uptime = time.time() - start_time
    #print(f"** [{matrix_type}] assembly: {uptime:1.4e} sec")
    return Kg
#
def Kt_matrix(elements,
              nodes,
              Un,
              plane2D: bool):
              #sparse:bool = True):
    """
    elements: Class
    nodes: Class
    Un: global node displacement
    plane2D: True/False
    sparse : True/False

    Returns:
        The Tangent stiffness matrix [Kt]
    """
    matrix_type: str = "Kt"
    #start_time = time.time()
    Kt = assemble_KtKg_matrix(elements=elements,
                              nodes=nodes,
                              plane2D=plane2D,
                              matrix_type=matrix_type,
                              Un=Un)
    #if sparse:
    #    Kt = coo_matrix(Kt)
    #    Kt = Kt.tolil()
    #
    #uptime = time.time() - start_time
    #print(f"** [{matrix_type}] assembly: {uptime:1.4e} sec")
    return Kt
#
def Kt_matrix_R(elements,
                nodes,
                Fb,
                plane2D: bool):
                #sparse:bool = True):
    """
    elements: Class
    nodes: Class
    Un: global node displacement
    plane2D: True/False
    sparse : True/False

    Returns:
        The Tangent stiffness matrix [Kt]
    """
    matrix_type: str = "Kt"
    #start_time = time.time()
    Kt = assemble_KtKg_R_matrix(elements=elements,
                                nodes=nodes,
                                plane2D=plane2D,
                                matrix_type=matrix_type, 
                                Fb=Fb)
    #if sparse:
    #    Kt = coo_matrix(Kt)
    #    Kt = Kt.tolil()
    #
    #uptime = time.time() - start_time
    #print(f"** [{matrix_type}] assembly: {uptime:1.4e} sec")
    return Kt
#
# ------------------------------------------------------
#
def assemblyXX(element, Ka: np.array,
             ndof:int, matrix_type:str,
             plane2D:bool,
             Fb: list|tuple|np.array|None):
    """
    element : Class
    Ka : The stiffness matrix's container
    ndof : node's degree of freedom
    matrix_type: Ke, Km, Kg, Kt
    plane2D: True/False
    Fb: Beam force

    Returns:
        The update stiffness matrix's container [Ka]
    """
    # TODO : check applicable to all element type
    keg = getattr(element, matrix_type)(plane2D, Fb)
    idof, jdof = element.DoF
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
def assembly(element, Ka: np.array,
             ndof:int, matrix_type:str,
             plane2D:bool,
             Fb: list|tuple|np.array|None):
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
    keg = getattr(element, matrix_type)(plane2D, Fb)
    idof, jdof = element.DoF
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
# ------------------------------------------------------
#
def assemble_KeKm_matrix(elements, nodes,
                         plane2D: bool,
                         matrix_type: str):
                         #Fb: None):
    """
    Assemble the global stiffness matrix [Ka] based on elements' Stiffness or Mass  matrix
    -------------------------------------------------
    elements: class
    nodes: class
    matrix_type : Ke, Km
    Fb : None

    Returns:
        Ka : Global stiffness matrix
    """
    Fb = None
    nn = len(nodes.keys())
    #Ka = form_matrix(elements=elements,
    #                 nn=nn, plane2D=plane2D,
    #                 matrix_type=matrix_type,
    #                 item=Fb)
    #
    ndof: int = 6
    if plane2D:
        ndof: int = 3
    # Initialize a K matrix of zeros
    Ka = dok_matrix((nn * ndof, nn * ndof), dtype=np.float32)
    #Ka = np.zeros((nn * ndof, nn * ndof), dtype=np.float64)
    with ThreadPoolExecutor() as executor:
        for key, element in elements.items():
            executor.submit(assembly, element, Ka, ndof, matrix_type, plane2D, Fb)
    #
    return Ka
#
#
def form_matrixXX(elements, nn: int, plane2D: bool,
                matrix_type: str, item: list | tuple | None):
    """
    Form global stiffness matrix [K] based on elements' Stiffness or Mass  matrix

    elements :
    nn  : node number
    ndof : node degree of freedom
    matrix_type: Ke/Km
    item : None

    :Returns
        Ka : global stiffness matrix
    """
    ndof: int = 6
    if plane2D:
        ndof: int = 3
    # Initialize a K matrix of zeros
    Ka = np.zeros((nn * ndof, nn * ndof), dtype=np.float64)
    with ThreadPoolExecutor() as executor:
        for key, element in elements.items():
            executor.submit(assembly, element, Ka, ndof, matrix_type, plane2D, item)
    return Ka
#
#
# ------------------------------------------------------
#
def form_KtKg_matrixX(elements, nn: int, plane2D: bool,
                     matrix_type: str, item: list | tuple | None):
    """
    Form global stiffness matrix [K] based on elements' Tangent or Geometric  matrix

    elements :
    nn  : node number
    ndof : node degree of freedom
    matrix_type: Kt/Kg
    item: Beam force

    :Returns
        Ka : global stiffness matrix
    """
    ndof = 6
    if plane2D:
        ndof = 3
    #
    # Initialize a K matrix of zeros
    Ka = np.zeros((nn * ndof, nn * ndof), dtype=np.float64)
    with ThreadPoolExecutor() as executor:
        for key, element in elements.items():
            nodes = element.connectivity
            # ---------------------------------------------
            # Beam's end forces
            nd_global = np.concatenate((item.loc[nodes[0]],
                                        item.loc[nodes[1]]), axis=None)
            #
            executor.submit(assembly, element, Ka, ndof, matrix_type, plane2D, nd_global)
            # ---------------------------------------------
    return Ka


#
#
def assemble_KtKg_matrix(elements, nodes,
                         plane2D: bool,
                         matrix_type: str):
                         #Un: list | tuple):
    """
    Assemble the stiffness matrix [Ka] based on elements' Tangent or Geometric  matrix
    -------------------------------------------------
    elements: class
    nodes: class
    matrix_type : Kt, Kg
    Un : node global displacement

    Returns:
        Ka : Global stiffness matrix
    """
    nn = len(nodes.keys())
    #Ka = form_KtKg_matrix(elements=elements,
    #                      nn=nn, plane2D=plane2D,
    #                      matrix_type=matrix_type, item=Fb)
    #
    #
    ndof = 6
    if plane2D:
        ndof = 3
    #
    Ka = dok_matrix((nn * ndof, nn * ndof), dtype=np.float32)
    for element in elements:
        #
        assembly(element, Ka, ndof, matrix_type,
                 plane2D=plane2D, Un=element.du)    
    #
    #beams = elements._beams
    #
    # Initialize a K matrix of zeros
    #Ka = np.zeros((nn * ndof, nn * ndof), dtype=np.float64)
    #Ka = dok_matrix((nn * ndof, nn * ndof), dtype=np.float32)
    #with ThreadPoolExecutor() as executor:
    #    # ---------------------------------------------
    #    # Beam section    
    #    for key, beam in beams.items():
    #        nodes = beam.connectivity
    #        Un_global = np.concatenate((Un.loc[nodes[0]],
    #                                    Un.loc[nodes[1]]), axis=None)
    #        #
    #        #beam_mod = beam.deformed(du=Un_global, R=R)
    #        #
    #        executor.submit(assembly, beam, Ka, ndof,
    #                        matrix_type, plane2D, Un_global)
    #    # ---------------------------------------------
    #    # TODO : implement plates
    #
    return Ka
#
def assemble_KtKg_R_matrix(elements, nodes,
                           plane2D: bool,
                           matrix_type: str, 
                           Fb: list | tuple):
    """
    Assemble the stiffness matrix [Ka] based on elements' Tangent or Geometric  matrix
    -------------------------------------------------
    elements: class
    nodes: class
    matrix_type : Kt, Kg
    Fb : Element internal force

    Returns:
        Ka : Global stiffness matrix
    """
    nn = len(nodes.keys())
    #Ka = form_KtKg_matrix(elements=elements,
    #                      nn=nn, plane2D=plane2D,
    #                      matrix_type=matrix_type, item=Fb)
    #
    #
    ndof = 6
    if plane2D:
        ndof = 3
    #
    #
    Ka = dok_matrix((nn * ndof, nn * ndof), dtype=np.float32)
    for element in elements:
        #
        assembly(element, Ka, ndof, matrix_type,
                 plane2D=plane2D, Fb=Fb[element.name])
    #
    #
    # Initialize a K matrix of zeros
    #Ka = np.zeros((nn * ndof, nn * ndof), dtype=np.float64)
    #Ka = dok_matrix((nn * ndof, nn * ndof), dtype=np.float32)
    #with ThreadPoolExecutor() as executor:
    #    # ---------------------------------------------
    #    # Beam section    
    #    for beam in elements:
    #        #nodes = beam.connectivity
    #        #Un_global = np.concatenate((Un.loc[nodes[0]],
    #        #                            Un.loc[nodes[1]]), axis=None)
    #        #
    #        #beam_mod = beam.deformed(du=Un_global, R=R)
    #        #
    #        executor.submit(assembly, beam, Ka, ndof,
    #                        matrix_type, plane2D, beam.du)
    #    # ---------------------------------------------
    #    # TODO : implement plates
    #
    return Ka
