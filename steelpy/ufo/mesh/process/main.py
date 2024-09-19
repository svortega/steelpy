#
# Copyright (c) 2009 steelpy
#

# Python stdlib imports
from __future__ import annotations
from concurrent.futures import ThreadPoolExecutor
import time
#from itertools import chain

# package imports
import numpy as np
from scipy.sparse import coo_matrix

#
def Ke_matrix(elements,
              nodes,
              plane2D: bool,
              sparse: bool = True,
              mitem:str = "Ke"):
    """
    Stiffness matrix

    elements:
    nodes:
    dof : ['x', 'y', 'z', 'rx', 'ry', 'rz']
    solver :
    condensed : Matrix with dof = 0 removed
    """
    start_time = time.time()
    Ka = assemble_matrix(elements=elements,
                         nodes=nodes,
                         plane2D=plane2D,
                         mitem=mitem, item=None)
    if sparse:
        Ka = coo_matrix(Ka)
    #
    uptime = time.time() - start_time
    print(f"** [{mitem}] assembly: {uptime:1.4e} sec")
    return Ka


#
def Km_matrix(elements,
              nodes,
              plane2D: bool,
              sparse:bool = True,
              mitem:str = "Km"):
    """ Mass matrix"""
    start_time = time.time()
    Ka = assemble_matrix(elements=elements,
                         nodes=nodes,
                         plane2D=plane2D,
                         mitem=mitem, item=None)
    if sparse:
        Ka = coo_matrix(Ka)
    #
    uptime = time.time() - start_time
    print(f"** [{mitem}] assembly: {uptime:1.4e} sec")
    return Ka


#
def Kg_matrix(elements,
              nodes, D,
              plane2D: bool,
              sparse:bool = True,
              mitem:str="Kg"):
    """ Geometric stiffness matrix """
    start_time = time.time()
    Kg = assemble_Gmatrix(elements=elements,
                          nodes=nodes,
                          plane2D=plane2D,
                          mitem=mitem, item=D)
    if sparse:
        Kg = coo_matrix(Kg)
    #
    uptime = time.time() - start_time
    print(f"** [{mitem}] assembly: {uptime:1.4e} sec")
    return Kg
#
def Kt_matrix(elements,
              nodes, D,
              plane2D: bool,
              sparse:bool = True,
              mitem:str="Kt"):
    """ Tangent stiffness matrix """
    start_time = time.time()
    Kt = assemble_Gmatrix(elements=elements,
                          nodes=nodes,
                          plane2D=plane2D,
                          mitem=mitem, item=D)
    if sparse:
        Kt = coo_matrix(Kt)
    #
    uptime = time.time() - start_time
    print(f"** [{mitem}] assembly: {uptime:1.4e} sec")
    return Kt
#
#
# ------------------------------------------------------
#
def assembly(element, Ka: np.array,
             ndof:int, mitem:str,
             plane2D:bool, item: list|tuple|np.array|None):
    """
    """
    # TODO : check applicable to all element type
    keg = getattr(element, mitem)(plane2D, item)
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
def form_matrix(elements, nn: int, plane2D: bool,
                 mitem: str, item: list|tuple|None):
    """
    Global system stiffness matrix

    elements :
    nn  : node number
    ndof : node degree of freedom
    :return
    Ka : global stiffness matrix
    """
    ndof: int = 6
    if plane2D:
        ndof: int = 3
    # Initialize a K matrix of zeros
    Ka = np.zeros((nn * ndof, nn * ndof), dtype=np.float64)
    with ThreadPoolExecutor() as executor:
        for key, element in elements.items():
            executor.submit(assembly, element, Ka, ndof, mitem, plane2D, item)
    #
    return Ka
#
#
def form_Gmatrix(elements, nn:int, plane2D: bool,
                 mitem: str, item:list|tuple|None):
    """
    Global system stiffness matrix 
    
    elements : 
    nn  : node number
    ndof : node degree of freedom
    :return
    Ka : global stiffness matrix
    """
    ndof = 6
    if plane2D:
        ndof = 3
    #
    # Initialize a K matrix of zeros
    Ka = np.zeros((nn*ndof, nn*ndof), dtype=np.float64)
    with ThreadPoolExecutor() as executor:
        for key, element in elements.items():
            nodes = element.connectivity
            # ---------------------------------------------
            # displacement
            nd_global = np.concatenate((item.loc[nodes[0]],
                                        item.loc[nodes[1]]), axis=None)
            #
            executor.submit(assembly, element, Ka, ndof, mitem, plane2D, nd_global)
            # ---------------------------------------------
    # 
    return Ka
#
#
# ------------------------------------------------------
#
def assemble_matrix(elements, nodes,
                    plane2D: bool,
                    mitem: str, item:list|tuple|None):
    """
    Asseable the element matrices
    -------------------------------------------------
    Ka : stiffness matrix
    jbc : nodes freedom
    ndof : node degree of freedom
    mitem : matrix item
    item : None
    """
    nn = len(nodes.keys())
    Ka = form_matrix(elements=elements,
                     nn=nn, plane2D=plane2D,
                     mitem=mitem, item=item)
    return Ka
#
def assemble_Gmatrix(elements, nodes,
                     plane2D: bool,
                     mitem: str, item:list|tuple|None):
    """
    Asseable the element matrices
    -------------------------------------------------
    Ka : stiffness matrix
    jbc : nodes freedom
    ndof : node degree of freedom
    mitem : matrix item
    item :
    """
    nn = len(nodes.keys())
    Ka = form_Gmatrix(elements=elements,
                      nn=nn, plane2D=plane2D,
                      mitem=mitem, item=item)
    return Ka
#
