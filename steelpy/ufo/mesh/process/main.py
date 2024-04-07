#
# Copyright (c) 2009 steelpy
#

# Python stdlib imports
from __future__ import annotations

# package imports
#steelpy.f2uModel.mesh
#from steelpy.f2uModel.mesh.process.matrix.Kassemble import assemble_banded_Kmatrix
#
#from steelpy.utils.math.operations import remove_column_row
from steelpy.ufo.mesh.process.operations import assemble_matrix, assemble_Gmatrix #, assemble_Lmatrix


#

def Ke_matrix(elements,
              nodes,
              ndof: int = 6,
              #condensed: bool = True,
              sparse: bool = True):
    """
    Stiffness matrix

    elements:
    nodes:
    ndof :
    solver :
    condensed : Matrix with dof = 0 removed
    """
    #jbc = nodes.jbc()
    #neq = nodes.neq()
    #nn = len(nodes.keys())
    #
    # Banded matrix
    #if solver == 'banded':
    #    iband = elements.max_bandwidth(jbc=jbc)
    #    aa =  assemble_banded_Kmatrix(elements=elements , jbc=jbc,
    #                                  neq=neq, iband=iband)
    #else:
    # numpy matrix
    #
    Ka = assemble_matrix(elements=elements,
                         nodes=nodes,
                         ndof=ndof,
                         mitem="Ke", item=None)
    #
    #
    if sparse:
        from scipy.sparse import coo_matrix
        Ka = coo_matrix(Ka)

    #elif condensed:
    #    #dof_index = nodes.DOF_unreleased()
    #    jbcc = jbc.stack().values
    #    index = list(reversed([i for i, item in enumerate(jbcc)
    #                           if item == 0]))
    #    for i in index:
    #        Ka = remove_column_row(Ka, i, i)
    #
    return Ka


#
def Km_matrix(elements,
              nodes,
              ndof: int = 6,
              sparse: bool = True):
    """ Mass matrix"""
    Ka = assemble_matrix(elements=elements,
                         nodes=nodes,
                         ndof=ndof,
                         mitem="Km", item=None)
    #
    #
    if sparse:
        from scipy.sparse import coo_matrix
        Ka = coo_matrix(Ka)
    return Ka


#
def Kg_matrix(elements,
              nodes, D,
              ndof: int = 6,
              sparse: bool = True):
    """ Geometric stiffness matrix """
    Kg = assemble_Gmatrix(elements=elements,
                          nodes=nodes,
                          ndof=ndof,
                          mitem="Kg", item=D)
    #
    if sparse:
        from scipy.sparse import coo_matrix
        Kg = coo_matrix(Kg)
    return Kg
#
def Kt_matrix(elements,
              nodes, D,
              ndof: int = 6,
              sparse: bool = True):
    """ Tangent stiffness matrix """
    Kt = assemble_Gmatrix(elements=elements,
                          nodes=nodes,
                          ndof=ndof,
                          mitem="Kt", item=D)
    #
    if sparse:
        from scipy.sparse import coo_matrix
        Kt = coo_matrix(Kt)
    return Kt
#
