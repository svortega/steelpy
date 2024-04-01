# 
# Copyright (c) 2009 steelpy
#
from __future__ import annotations
#
# Python stdlib imports
import time
from itertools import chain
#from multiprocessing import Process, Manager

# package imports
#from steelpy.utils.math.operations import remove_column_row
import numpy as np

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
    # 
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
        #nodes = element.nodes
        nodes = element.connectivity
        Tb = element.T
        # ---------------------------------------------
        # displacement
        nd_global = np.concatenate((item.loc[nodes[0]],
                                    item.loc[nodes[1]]), axis=None)
        # ---------------------------------------------
        # convert global end-node disp in beam's local system
        # [x,y,z,rx,ry,rz]
        nd_local = Tb @ nd_global        
        P = nd_local[6]-nd_local[0]
        #
        keg = getattr(element, mitem)(P)
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
def sparse_matrix(elements, nn:int, ndof: int, mitem: str):
    """ """
    from scipy.sparse import coo_matrix
    #
    neq = nn*ndof
    # Initialize a K matrix of zeros
    Ka = np.zeros((neq, neq), dtype=np.float64)    
    #
    ki = []
    row = []
    col = []    
    for key, element in elements.items():
        # TODO : check applicable to all element type
        keg = getattr(element, mitem)
        idof, jdof = element.DoF #
        # node and corresponding dof (start, end)
        # Start i
        niqi = idof*ndof
        niqj = niqi + ndof
        # End j
        njqi = jdof*ndof
        njqj = njqi + ndof
        #
        q2nd = [item for item in range(niqi, niqj, 1)]
        q1st = [item for item in range(njqi, njqj, 1)]
        #
        #2nd
        #ki.append(keg[:ndof, :ndof])
        #row.append(q2nd)
        #col.append(q2nd)
        # 1st
        #ki.append(keg[:ndof, ndof:2*ndof])
        #row.append(q2nd)
        #col.append(q1st)
        # 3rd
        #ki.append(keg[ndof:2*ndof, :ndof])
        #row.append(q1st)
        #col.append(q2nd)
        # 4th
        #ki.append(keg[ndof:2*ndof, ndof:2*ndof])
        #row.append(q1st)
        #col.append(q1st)
        print('-->')
        #ki.append(np.concatenate([keg[:ndof, :ndof],
        #                          keg[:ndof, ndof:2*ndof],
        #                          keg[ndof:2*ndof, :ndof],
        #                          keg[ndof:2*ndof, ndof:2*ndof]]))
        #
        # assemble global stiffness matrix, quadrant 1 to 4
        Ka[niqi:niqj, niqi:niqj] += keg[:ndof, :ndof]             # 2nd
        Ka[niqi:niqj, njqi:njqj] += keg[:ndof, ndof:2*ndof]       # 1st
        Ka[njqi:njqj, niqi:niqj] += keg[ndof:2*ndof, :ndof]       # 3rd
        Ka[njqi:njqj, njqi:njqj] += keg[ndof:2*ndof, ndof:2*ndof] # 4th        
        #
        row.append(list(chain(q2nd, q2nd, q1st, q1st)))
        col.append(list(chain(q2nd, q1st, q2nd, q1st)))
    #
    K = coo_matrix(Ka)
    #K = coo_matrix((ki, (row, col)), shape=(neq, neq))
    return K
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
    print(f"** Processing Global [{mitem}] Matrix")
    start_time = time.time()
    #
    nn = len(nodes.keys())
    #nn = len(jbc)
    #neq = nn*ndof
    Ka = form_matrix(elements=elements,
                     nn=nn,
                     ndof=ndof,
                     mitem=mitem, item=item)
    #
    #Kas = sparse_matrix(elements, nn, ndof, mitem)
    #
    #jbcc = jbc.stack().values
    #index = list(reversed([i for i, jbitem in enumerate(jbcc)
    #                       if jbitem == 0]))
    #
    #for i in index:
    #    Ka = remove_column_row(Ka, i, i)
    #
    uptime = time.time() - start_time
    print(f"** [{mitem}] assembly Finish Time: {uptime:1.4e} sec")
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
    print(f"** Processing Global [{mitem}] Matrix")
    start_time = time.time()
    #
    nn = len(nodes.keys())
    Ka = form_Gmatrix(elements=elements,
                      nn=nn,
                      ndof=ndof,
                      mitem=mitem, item=item)
    #
    #Kas = sparse_matrix(elements, nn, ndof, mitem)
    #
    #jbcc = jbc.stack().values
    #index = list(reversed([i for i, jbitem in enumerate(jbcc)
    #                       if jbitem == 0]))
    #
    #for i in index:
    #    Ka = remove_column_row(Ka, i, i)
    #
    uptime = time.time() - start_time
    print(f"** [{mitem}] assembly Finish Time: {uptime:1.4e} sec")
    return Ka
#
#