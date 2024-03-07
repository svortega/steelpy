#
# Copyright (c) 2009 steelpy
#

# Python stdlib imports
from __future__ import annotations

# package imports
#steelpy.f2uModel.mesh
#from steelpy.f2uModel.mesh.process.matrix.Kassemble import assemble_banded_Kmatrix
#
from steelpy.utils.math.operations import remove_column_row
from steelpy.ufo.mesh.process.operations import assemble_matrix
#

def Kmatrix(elements,
            nodes, 
            ndof: int,
            condensed: bool = True,
            sparse: bool = False):
    """
    elements:
    nodes:
    ndof :
    solver :
    condensed : Matrix with dof = 0 removed
    """
    jbc = nodes.jbc()
    #neq = nodes.neq()    
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
                         jbc=jbc, #neq=neq,
                         ndof=ndof, mitem="K")
    #
    #
    if sparse:
        from scipy.sparse import coo_matrix
        Ka = coo_matrix(Ka)
    
    elif condensed:
        #dof_index = nodes.DOF_unreleased()
        jbcc = jbc.stack().values
        index = list(reversed([i for i, item in enumerate(jbcc)
                               if item == 0]))
        for i in index:
            Ka = remove_column_row(Ka, i, i)
    #
    return Ka
#
def Mmatrix(elements, jbc, neq, ndof: int, solver: str|None = None):
    """ """
    aa = assemble_matrix(elements=elements, jbc=jbc, neq=neq, ndof=ndof, mitem="M")
    #
    return aa    
#
#
#
def Gmatrix(elements, jbc, neq, plane, solver: str|None = None):
    """ """
    aa = assemble_matrix(elements=elements, jbc=jbc, neq=neq, plane=plane, mitem="Kg")
    #
    return aa    
#
#