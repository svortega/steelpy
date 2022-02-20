#
# Copyright (c) 2009-2021 steelpy
#
# Python stdlib imports
from array import array
from copy import copy
from math import fsum
import pickle
from typing import List #, NamedTuple, Union
from itertools import chain
import time
#
# package imports
from steelpy.trave3D.processor.operations import get_deflection
from steelpy.process.math.operations import to_matrix, zeros_vector, matAbd, trns_3Dv
from steelpy.process.math.vector import Vector
import numpy as np
#
#
# ---------------------------------------------
# solver using numpy
#
def solve_displacement_numpy(elements, nodes, boundaries,
                             basic_load, load_combination):
    """
    """
    print("** Calculating Joint Displacements")
    print("** reloaded [k] & {p}")    
    #
    file = open( "stfmx.f2u", "rb" )
    jbc = pickle.load( file )
    stf = pickle.load( file )
    file.close()
    #
    #
    ndof = 6    # nodal degrees of freedom
    #en = 2      # number of nodes per element    
    # global system stiffness matrix
    nn = len(nodes)
    #
    # global system force matrix consisting of summation of loads
    #Fs = np.zeros((nn*ndof, 1), dtype=np.float64)    
    #
    jbcc = list(chain.from_iterable(jbc))
    #
    print("** Calculating Joint Displacements")
    print("** reloaded [k] & {p}")
    #
    #fmemb = {}
    memf_basic = {}
    basic_res = {}
    for key, item in basic_load.items():
        load_res = basic_load.get_basic_load(key, elements, nodes, boundaries)
        # member end load in local system
        memf_basic[item.title] = load_res.member_load
        # node load in global system
        #nload = np.array(list(chain.from_iterable(load_res.nodal_load)))
        #
        # node load in global system
        nodal_load = load_res.nodal_load        
        # reduce degree of fredooms
        nloads = [nodal_load[i][j] 
                  for i in range(len(jbc)) for j in range(6) 
                  if jbc[i][j] != 0]
        nload = np.array(nloads)
        #
        #Fs[:,0] += nload
        basic = np.linalg.solve(stf, nload)
        #
        #if not np.allclose(np.dot(stf, basic), nload): #, rtol=1e-05, atol=1e-08
        #    raise RuntimeError('Solution fail')
        #
        basic = [basic[ieqnum - 1] if ieqnum != 0 else ieqnum
                 for ieqnum in jbcc]
        basic_res[item.title] = Vector(basic)        
    #
    #
    comb_res, memf_comb = load_combination.solve_combinations(basic_res, memf_basic)
    #
    with open("elemfout.f2u", "wb") as f:
        pickle.dump(memf_basic, f)
        pickle.dump(memf_comb, f)
    #
    # all supports are initially assumed linear
    #Xs = np.linalg.solve(Ks, Fs)
    #
    #flag = np.allclose(np.dot(Ks, Xs), Fs)
    #if flag:
    #    print('** Solution ok')
    #else:
    #    raise RuntimeError('Solution fail')
    #
    #Xb = Xs.reshape(len(nodes), 6)
    #
    #return Xs
    #
    deflection = get_deflection(nodes, basic_res, comb_res)
    print("** Finished Calculating Joint Displacements")
    return deflection
#
#