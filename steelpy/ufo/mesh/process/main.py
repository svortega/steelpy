#
# Copyright (c) 2009 steelpy
#

# Python stdlib imports
from __future__ import annotations
import time
from itertools import chain

# package imports
#steelpy.f2uModel.mesh
#from steelpy.f2uModel.mesh.process.matrix.Kassemble import assemble_banded_Kmatrix
#
#from steelpy.utils.math.operations import remove_column_row
#from steelpy.ufo.mesh.process.operations import assemble_matrix, assemble_Gmatrix #, assemble_Lmatrix

import numpy as np

#
def Ke_matrix(elements,
              nodes,
              ndof: int = 6,
              #condensed: bool = True,
              sparse: bool = True,
              mitem:str = "Ke"):
    """
    Stiffness matrix

    elements:
    nodes:
    ndof :
    solver :
    condensed : Matrix with dof = 0 removed
    """
    start_time = time.time()
    #mitem = "Ke"
    Ka = assemble_matrix(elements=elements,
                         nodes=nodes,
                         ndof=ndof,
                         mitem=mitem, item=None)
    if sparse:
        from scipy.sparse import coo_matrix
        Ka = coo_matrix(Ka)
    #
    uptime = time.time() - start_time
    print(f"** [{mitem}] assembly: {uptime:1.4e} sec")
    return Ka


#
def Km_matrix(elements,
              nodes,
              ndof:int = 6,
              sparse:bool = True,
              mitem:str = "Km"):
    """ Mass matrix"""
    start_time = time.time()
    Ka = assemble_matrix(elements=elements,
                         nodes=nodes,
                         ndof=ndof,
                         mitem=mitem, item=None)
    if sparse:
        from scipy.sparse import coo_matrix
        Ka = coo_matrix(Ka)
    #
    uptime = time.time() - start_time
    print(f"** [{mitem}] assembly: {uptime:1.4e} sec")
    return Ka


#
def Kg_matrix(elements,
              nodes, D,
              ndof:int = 6,
              sparse:bool = True,
              mitem:str="Kg"):
    """ Geometric stiffness matrix """
    start_time = time.time()
    Kg = assemble_Gmatrix(elements=elements,
                          nodes=nodes,
                          ndof=ndof,
                          mitem=mitem, item=D)
    if sparse:
        from scipy.sparse import coo_matrix
        Kg = coo_matrix(Kg)
    #
    uptime = time.time() - start_time
    print(f"** [{mitem}] assembly: {uptime:1.4e} sec")
    return Kg
#
def Kt_matrix(elements,
              nodes, D,
              ndof:int = 6,
              sparse:bool = True,
              mitem:str="Kt"):
    """ Tangent stiffness matrix """
    start_time = time.time()
    Kt = assemble_Gmatrix(elements=elements,
                          nodes=nodes,
                          ndof=ndof,
                          mitem=mitem, item=D)
    if sparse:
        from scipy.sparse import coo_matrix
        Kt = coo_matrix(Kt)
    #
    uptime = time.time() - start_time
    print(f"** [{mitem}] assembly: {uptime:1.4e} sec")
    return Kt
#
#
# ------------------------------------------------------
#
#
def form_matrix(elements, nn:int, ndof: int,
                mitem: str, item):
    """
    Global system stiffness matrix 
    
    elements : 
    nn  : node number
    ndof : node degree of freedom
    :return
    Ka : global stiffness matrix
    """
    # Initialize a K matrix of zeros
    Ka = np.zeros((nn*ndof, nn*ndof), dtype=np.float64)
    for key, element in elements.items():
        # TODO : check applicable to all element type
        keg = getattr(element, mitem)(item)
        idof, jdof = element.DoF
        # node and corresponding dof (start, end)
        niqi, niqj = idof*ndof, idof*ndof + ndof
        njqi, njqj = jdof*ndof, jdof*ndof + ndof
        # assemble global stiffness matrix, quadrant 1 to 4
        Ka[niqi:niqj, niqi:niqj] += keg[:ndof, :ndof]             # 2nd
        Ka[niqi:niqj, njqi:njqj] += keg[:ndof, ndof:2*ndof]       # 1st
        Ka[njqi:njqj, niqi:niqj] += keg[ndof:2*ndof, :ndof]       # 3rd
        Ka[njqi:njqj, njqi:njqj] += keg[ndof:2*ndof, ndof:2*ndof] # 4th
    return Ka
#
#
def form_Gmatrix(elements, nn:int, ndof: int,
                 mitem: str, item):
    """
    Global system stiffness matrix 
    
    elements : 
    nn  : node number
    ndof : node degree of freedom
    :return
    Ka : global stiffness matrix
    """
    # Initialize a K matrix of zeros
    Ka = np.zeros((nn*ndof, nn*ndof), dtype=np.float64)
    for key, element in elements.items():
        nodes = element.connectivity
        # ---------------------------------------------
        # displacement
        nd_global = np.concatenate((item.loc[nodes[0]],
                                    item.loc[nodes[1]]), axis=None)
        # ---------------------------------------------
        # convert global end-node disp in beam's local system
        # [x,y,z,rx,ry,rz]
        #nd_local = Tb @ nd_global
        #P = nd_local[6]-nd_local[0]
        #
        # ---------------------------------------------
        # get matrix
        keg = getattr(element, mitem)(nd_global)
        #
        # ---------------------------------------------
        # Assemble matrix
        idof, jdof = element.DoF
        # node and corresponding dof (start, end)
        niqi, niqj = idof*ndof, idof*ndof + ndof
        njqi, njqj = jdof*ndof, jdof*ndof + ndof
        # assemble global stiffness matrix, quadrant 1 to 4
        Ka[niqi:niqj, niqi:niqj] += keg[:ndof, :ndof]             # 2nd
        Ka[niqi:niqj, njqi:njqj] += keg[:ndof, ndof:2*ndof]       # 1st
        Ka[njqi:njqj, niqi:niqj] += keg[ndof:2*ndof, :ndof]       # 3rd
        Ka[njqi:njqj, njqi:njqj] += keg[ndof:2*ndof, ndof:2*ndof] # 4th
    # 
    return Ka
#
#
def form_Lmatrix(elements, nn: int, ndof: int,
                 mitem: str, item):
    """
    Global system stiffness matrix

    elements :
    nn  : node number
    ndof : node degree of freedom
    miten : str
    item:
    :return
    Ka : global stiffness matrix
    """
    # Initialize a K matrix of zeros
    Ka = np.zeros((nn * ndof, nn * ndof), dtype=np.float64)
    for key, element in elements.items():
        # nodes = element.nodes
        nodes = element.connectivity
        #Tb = element.T
        # ---------------------------------------------
        # displacement
        nd_global = np.concatenate((item.loc[nodes[0]],
                                    item.loc[nodes[1]]), axis=None)
        # ---------------------------------------------
        # get matrix
        keg = getattr(element, mitem)(nd_global)
        #
        #
        idof, jdof = element.DoF
        # node and corresponding dof (start, end)
        niqi, niqj = idof * ndof, idof * ndof + ndof
        njqi, njqj = jdof * ndof, jdof * ndof + ndof
        # assemble global stiffness matrix, quadrant 1 to 4
        Ka[niqi:niqj, niqi:niqj] += keg[:ndof, :ndof]  # 2nd
        Ka[niqi:niqj, njqi:njqj] += keg[:ndof, ndof:2 * ndof]  # 1st
        Ka[njqi:njqj, niqi:niqj] += keg[ndof:2 * ndof, :ndof]  # 3rd
        Ka[njqi:njqj, njqi:njqj] += keg[ndof:2 * ndof, ndof:2 * ndof]  # 4th
    #
    return Ka
#
# ------------------------------------------------------
#
def assemble_matrix(elements, nodes,
                    ndof: int,
                    mitem: str, item):
    """
    Asseable the element matrices
    -------------------------------------------------
    Ka : stiffness matrix
    jbc : nodes freedom
    ndof : node degree of freedom
    mitem : matrix item
    """
    #print(f"** Processing Global [{mitem}] Matrix")
    #start_time = time.time()
    nn = len(nodes.keys())
    Ka = form_matrix(elements=elements,
                     nn=nn, ndof=ndof,
                     mitem=mitem, item=item)
    #
    #uptime = time.time() - start_time
    #print(f"** [{mitem}] assembly: {uptime:1.4e} sec")
    return Ka
#
def assemble_Gmatrix(elements, nodes,
                    ndof: int,
                    mitem: str, item):
    """
    Asseable the element matrices
    -------------------------------------------------
    Ka : stiffness matrix
    jbc : nodes freedom
    ndof : node degree of freedom
    mitem : matrix item
    """
    #print(f"** Processing Global [{mitem}] Matrix")
    #start_time = time.time()
    nn = len(nodes.keys())
    Ka = form_Gmatrix(elements=elements,
                      nn=nn, ndof=ndof,
                      mitem=mitem, item=item)
    #uptime = time.time() - start_time
    #print(f"** [{mitem}] assembly: {uptime:1.4e} sec")
    return Ka
#
