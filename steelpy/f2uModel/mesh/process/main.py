#
# Copyright (c) 2009-2023 fem2ufo
#

# Python stdlib imports
from __future__ import annotations

# package imports
#steelpy.f2uModel.mesh
from steelpy.f2uModel.mesh.process.matrix.Kassemble import assemble_Kmatrix_np, assemble_banded_Kmatrix
#
#

def Kmatrix(elements, jbc, neq, solver: str|None = None,
            m2D:bool = False):
    """ """
    
    # Banded matrix
    if solver == 'banded':
        iband = elements.max_bandwidth(jbc=jbc)
        aa =  assemble_banded_Kmatrix(elements=elements , jbc=jbc,
                                      neq=neq, iband=iband)
    else:
        # numpy matrix
        aa = assemble_Kmatrix_np(elements=elements, jbc=jbc, neq=neq, m2D=m2D)
    #
    return aa
#
def Mmatrix():
    """ """
    pass
#
#
#
def Gmatrix():
    """ """
    pass
#
#