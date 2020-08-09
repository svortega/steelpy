# 
# Copyright (c) 2009-2020 fem2ufo
#
# Python stdlib imports
#from array import array
#from itertools import chain
#from copy import  copy
import math as math
import pickle
#
# package imports
#from steelpy.frame3D.processor.operations import zeros, to_matrix
#
#
#
# --------------------
# output 
# --------------------
#
def print_member_forces(elements, element_force):
    """
    """
    print("***************************************")
    print("")
    print("** Reloaded Elements Forces")
    print("")
    #file = open( "elemout.f2u", "rb" )
    #basic_load = pickle.load(file)
    for key, item in element_force["basic"].items():
        print("-- Basic Load {:} {:}".format(item.name, item.title))
        print_member_end_forces(elements, item.items["global"], item.items["local"])
        print("")
    #
    try:
        #comb_load = pickle.load(file)
        print("***************************************")
        for key, item in element_force["combination"].items():
            print("-- Load Combination {:} {:}".format(item.name, item.title))
            print_member_end_forces(elements, item.items["global"], item.items["local"])
            print("")
    except EOFError:
        pass
    #file.close()
#
#                   
def print_member_end_forces(elements, load_global, load_local):
    """
    """
    for index in load_global.keys():
        element = elements[index]
        print ( "\n" )
        print ( " Element number: {:}   Type: {:}"
                .format ( element.name, element.type ) )
        #
        nodes = element.connectivity
        #
        print ( "{:>7}|{:>9} Member Load {:>10}|{:>7} Member Moments {:>9}|"
                .format ( "", "", "", "", "" ) )
        print ( "Node   |    axial     shearY     shearZ |   torque    momentY    momentZ |" )
        print ( "{:6d} |".format( nodes[0] ), end='' )
        for i in range ( 6 ):
            print ( "{: 1.3e} ".format ( load_local[index][i] ), end='' )
        print ( "" )
        print ( "{:6d} |".format( nodes[1] ), end='' )
        for i in range ( 6, 12 ):
            print ( "{: 1.3e} ".format ( load_local[index][i] ), end='' )
        #
        print ( "\n" )
        print ( "Loads in global coords" )
        print ( "{:6d} |".format( nodes[0] ), end='' )
        for i in range ( 6 ):
            print ( "{: 1.3e} ".format ( load_global[index][i] ), end='' )
        print ( "" )
        print ( "{:6d} |".format( nodes[1] ), end='' )
        for i in range ( 6, 12 ):
            print ( "{: 1.3e} ".format ( load_global[index][i] ), end='' )
    #
    print("\n")
#          
#
#
# -------------------------
#
#
def print_deflections(nodes, disp_result):
    """
    """
    print("")
    print("***************************************")
    print("** Reloaded Joint Displacements")
    #file = open( "nodout.f2u", "rb" )
    #basic = pickle.load(file)
    #
    print("")
    print("** Static Displacements")
    print("")
    for wk in disp_result["basic"].values():
        print("-- Basic Load {:} {:}".format(wk.name, wk.title))
        out_dis(wk.items, nodes)
        print("")
    #
    if disp_result["combination"]:
    #try:
        #comb = pickle.load(file)
        print("***************************************")
        print("** Static Displacements")
        print("")        
        for key, wk in disp_result["combination"].items():
            print("-- Load Combination {:} {:}".format(wk.name, wk.title))
            out_dis(wk.items, nodes)
    #except EOFError:
    #    print("** No Load Combinations")
    #
    #file.close()
#
def get_max_displacement(dispp):
    """ """
    columns = list(zip(*dispp))
    #
    maxval = []
    nodeitem = []
    print("")
    print("Maximum displacements")    
    for column in columns:
        nmax = [max(column), min(column)]
        maxval.append(nmax[0] if nmax[0] > abs(nmax[1]) else nmax[1])
        nodeitem.append(column.index(maxval[-1]))
    #
    print("node {:>10}".format(""), end='')
    for _node in nodeitem:
        print("{:}{:>10}".format(_node, ""), end='')
    #
    print("")
    print("value ", end='')
    for _value in maxval:
        print("{: 1.3e} ".format(_value), end='')
    #    
#
def out_dis(disp, nodes, str1=""):
    """
    Report the relative displacement of the structure 
    """
    # TODO: send this to mesh node?
    print("    node   x-disp     y-disp     z-disp     x-rot      y-rot      z-rot")
    #
    #dispp = to_matrix(disp, 6)
    #node_user = {key: node.number for key, node in nodes.items()}
    node_user = sorted([key for key in nodes.keys()])
    # calculate maximum values for easy reference
    #for key, node in nodes.items():
    for key in node_user:
        node = nodes[key]
        #index = node.number - 1
        print("{:8.0f}  ".format(node.name), end='')
        for _disp in disp[node.number]:
            print("{: 1.3e} ".format(_disp), end='')
        print("")
    #
    #get_max_displacement(dispp)
    #print("-->")
#
#
#
# -------------------------
# 
def eig_out(a, eigv, jbc, ivib, nnp, neq):
    """
    Report the mode shapes for eigenvalue problems
    """
    #
    # L82:
    #
    # open report.txt for app# end as #1
    print("CHOOSE nodal storage: ")
    print("     0 = return, 1 = mode shapes, 2 = nodal vectors, 3 = full <binary>")
    #istore = input(" --> ")
    istore = 1
    try: 
        istore = int(istore)
    except : 
        istore = 0
    #
    #iout = open('stadyn.out','a')
    #ilog = open('stadyn.log','a')
    #ilog.write(" \n")
    #ilog.write(" {:} :: istore  1 = shape, 2 = vector, 3 = Matrix <binary> \n".format(istore))
    #
    print("{:} :: istore  1 = shape, 2 = vector, 3 = Matrix <binary>".format(istore))
    
    if istore == 0 : 
        #iout.close()
        #ilog.close()        
        return
    #
    elif istore == 1 :
        print("")
        print("INPUT: No of modes to report")
        #neigen = input(" --> ")
        neigen = 5
        try: 
            neigen = int(neigen)
        except : 
            neigen = 1
        #
        if neigen  >  neq : 
            neigen = neq
        
        print(" {:} :: No of nodes".format(neigen))
        #ilog.write(" {:} :: No of nodes \n".format(neigen))
        
        #disp = np.zeros(nnp*6 +1, dtype = np.float64, order = 'F')
        disp = zeros(nnp*6)
        #dispp = rddisp(neq)
        # assign displacements to each node 
        #for idof in range(1, nnp*6 + 1):
        #    ieqnum = jbc[idof]
        #    if ieqnum > 0 : disp[idof] = dispp[ieqnum]
        #    else : disp[idof] = 0.0e0
        #L70:        
        #iout.write(" ")
        for nn in range(neigen):
            # fill in the displacement vector
            for idof in range(nnp*6):
                ieqnum = jbc[idof]
                if ieqnum > 0 : 
                    disp[idof] = a[ieqnum-1][nn]
                else : 
                    disp[idof] = 0.0
            # L70:
            
            # print the displacements 
            if ivib == 2 :
                freq = abs(eigv[nn])
                freq =  math.sqrt(freq)
                print("resonant freq: {: 1.4e} rad/s".format(freq))
                #iout.write(" \n")
                #iout.write("resonant freq: {: 1.4e} rad/s \n".format(freq))
                #
                str1 = "Mode shapes"
                out_dis(nnp, disp, str1)
            elif ivib == 1 :
                freq = eigv[nn]
                print("BUCKLING L0AD: {: 1.4e}".format(freq))
                #iout.write("BUCKLING L0AD: {: 1.4e} \n".format(freq))
                str1 = "Mode shapes"
                out_dis(nnp, disp, str1)
            # 23  format(1x, a, 3x, 1g13.5, 1x, a)
        # L20:
    #
    elif istore == 2 :
        #
        print("INPUT: No of vectors to report")
        neigen = raw_input(" --> ")
        try: neigen = int(neigen)
        except : neigen = 1  
        #
        if neigen > neq : neigen = neq
        
        ilog.write(" {:} :: No of vectors \n".format(neigen))
        print(" {:} :: No of vectors \n".format(neigen))
        #iout.write(" \n")
        iout.write("MODAL VECTORS \n")
        
        for nn in range(neigen):
            iout.write("{: 1.4e} ".format(float(eigv[nn])))
            for i in range(neq):
                iout.write("{: 1.4e} ".format(float(a[i, nn])))
            iout.write(" \n")
        # L80:
    #
    elif istore == 3 :
        for nn in range(neq):
            iout.write("{: 1.4e} ".format(float(eigv[nn])))
            for i in range(neq):
                iout.write("{: 1.4e} ".format(float(a[i, nn])))
            iout.write(" \n")
        # L84:
    # end if
    #
    #iout.close()
    #ilog.close()
    #
    print('-->')
#
#
def nodeout(_time,forc,disp,vel,acc,neq,nout,inod,idyn):
    """
    """
    # Dim inod(111,2) as integer
    # Dim  disp(neq)  as double
    # Dim  vel(neq)  as double
    # Dim  acc(neq) as double
    # Dim velout(111) as double
    #
    for n in range(nout):
        iprnode = inod[n][0]

        if inod[n][1] == 0 :
            velout[n] = disp[iprnode]
        elif inod[n][1] == 1 :
            velout[n] = vel[iprnode]
        elif inod[n][1] == 2 :
            velout[n] = acc[iprnode]
    #L72:
    # open report.txt for app# end as #1
    print ([_time, forc, velout[n]]  for n in range(1, nout))
    # 122   format(1x,8(g12.5, 1x) )
    # close #1:exit def
    # close #1:exit def
    #
    # close #1
    # end def
#
