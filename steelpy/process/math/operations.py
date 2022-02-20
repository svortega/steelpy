# 
# Copyright (c) 2009-2022 steelpy
#

# Python stdlib imports
from array import array
from math import fsum

# package imports
from steelpy.process.math.vector import Vector

# --------------------
# Matrix Operations       
# --------------------
#
#
def zeros(m, n=None, code:str = 'd'):
    """
    Create zero matrix
    """
    if n: 
        new_matrix = [array(code, [0 for row in range(n)]) 
                      for col in range(m)]
        #new_matrix = [[0 for row in range(n)] for col in range(m)]
        #new_matrix = np.zeros((m, n), dtype=np.float64, order='F')
    else: 
        new_matrix = array(code, [0 for row in range(m)])
    
    return new_matrix
#
def zeros_vector(m, n, code:str = 'd'):
    """
    Create zero matrix
    """
    new_matrix = [Vector([0 for row in range(n)]) for col in range(m)]
    return new_matrix
#
def ones(m, n=None, code:str = 'd'):
    """
    Create zero matrix
    """
    if n: 
        new_matrix = [array(code, [1 for row in range(n)]) 
                      for col in range(m)]        
        #new_matrix = [[1 for row in range(n)] for col in range(m)]
        #new_matrix = np.zeros((m, n), dtype=np.float64, order='F')
    else: 
        new_matrix = array(code, [1 for row in range(m)])
    
    return new_matrix
#
def to_matrix(l, n):
    return [l[i:i+n] for i in range(0, len(l), n)]
#
#
def matAbd(a, b, code:str = 'd'):
    """
    Multiply [a]{b} = {c}
    """
    n = min(len(a), len(b))
    c = array(code, [fsum([a[i][j] * b[j] for j in range(n)])
                     for i in range(n)])
    return c
#
def transposeM(matrix):
    """
    """
    # rt = list( map( list, zip( *r_matrix ) ) )
    return list(zip(*matrix))
#
#
def trns_3Dv(gloads, r_matrix):
    """
    Makes 3-d vector transformations.
    """
    lloads = zeros(12)    
    for i in range(0, 12, 3):
        lloads[i:i+3] = matAbd(r_matrix, gloads[i:i+3])
    #yy = [matAbd(r_matrix, gloads[i:i+3]) 
    #      for i in range(0, 12, 3)]
    return lloads
#
#
#
def Tmatrix(Rmat):
    """
    """
    tmx = zeros(12, 12)
    #tmx = Matrix(12, 12)
    for j1 in range ( 0, 12, 3 ):
        for k in range( 3 ):
            for l in range( 3 ):
                tmx[ j1 + k ][ j1 + l ]  = Rmat[ k ][ l ]
    return tmx
#