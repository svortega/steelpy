# Copyright (c) 2009-2021 fem2ufo

# Python stdlib imports
from typing import List, ClassVar, Dict, NamedTuple, Union
#
# package imports
from steelpy.f2uModel.load.operations.operations import duplicates, indices


#
class Results(NamedTuple):
    """ Basic load transfer"""
    name: int
    number:int
    title: str
    load_type: str
    items: List[List[float]]
# -------------------------
#
#
def print_deflections(disp_result):
    """
    """
    print('{:}'.format(52 * '-'))
    print("** Reloaded Joint Displacements")
    print('{:}'.format(52 * '-'))
    #print("")
    #print("** Static Displacements")
    print("")
    for wk in disp_result["basic"].values():
        print("-- Basic Load  Name: {:}  Number: {:}  Title: {:}  System: {:}"
              .format(wk.name, wk.number, wk.title, "global"))
        out_dis(wk.items)
        print("")
    #
    if disp_result["combination"]:
    #try:
        #comb = pickle.load(file)
        print('{:}'.format(52 * '-'))
        #print("** Static Displacements")
        print("")
        for key, wk in disp_result["combination"].items():
            print("-- Load Combination  Name: {:}  Number: {:}  Title: {:}  System: {:}"
                  .format(wk.name, wk.number, wk.title, "global"))
            out_dis(wk.items)
            print("")
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
def out_dis(disp, str1=""):
    """
    Report the relative displacement of the structure
    """
    # TODO: send this to mesh node?
    print("     node  x-disp     y-disp     z-disp     x-rot      y-rot      z-rot")
    for item in disp:
        print("{:9d} ".format(item[0]), end='')
        for _disp in item[1:]:
            print("{: 1.3e} ".format(_disp), end='')
        print("")
    #
    #get_max_displacement(dispp)
    #print("-->")
#
#
#
# --------------------
# output
# --------------------
#
def print_member_forces(element_force,
                        system:str="local"):
    """
    """
    print('{:}'.format(52 * '-'))
    print("** Reloaded Elements Forces")
    print('{:}'.format(52 * '-'))
    print("")
    #
    for key, item in element_force["basic"].items():
        print("-- Basic Load  Name: {:}  Number: {:}  Title: {:}  System: {:}"
              .format(item.name, item.number, item.title, system))
        mforces = [memb for memb in item.items if memb[3] == system]
        print_member_end_forces(mforces)
        print("")
    #
    try:
        print('{:}'.format(52 * '-'))
        for key, item in element_force["combination"].items():
            print("-- Load Combination  Name: {:}  Number: {:}  Title: {:}  System: {:}"
                  .format(item.name, item.number, item.title, system))
            mforces = [memb for memb in item.items if memb[3] == system]
            print_member_end_forces(mforces)
            print("")
    except EOFError:
        pass
#
#
def print_member_end_forces(element_forces):
    """
    """
    #dispp = list(chain.from_iterable(element_forces))
    forces = list(zip(*element_forces))
    dup = duplicates(forces[0])
    elements = indices(forces[0], dup)
    #print("{:>17}|{:>9} Member Load {:>10}|{:>7} Member Moments {:>7}|"
    #        .format ( "", "", "", "", "" ))
    print("  Element       Node  Fx{:} Fy{:} Fz{:} Mx{:} My{:} Mz"
          .format(8*" ", 8*" ", 8*" ", 8*" ", 8*" "))
    for key, items in elements.items():
        print("{:9d} ".format(key), end='')
        flag = 1
        for index in items:
            try:
                1/flag
                flag = 0
            except ZeroDivisionError:
                print("{:} ".format(9*" "), end='')
            print("{:10d} ".format(element_forces[index][1]), end='')
            for i in range(6):
                print("{: 1.3e} ".format(element_forces[index][i+4]), end='')
            #print("{:}".format(8*" "))
            print("")
    #print("\n")
#
def print_node_reactions(reactions):
    """
    :return:
    """
    print('{:}'.format(52 * '-'))
    print("** Reloaded Joint Reactions")
    print('{:}'.format(52 * '-'))
    print("")
    for wk in reactions["basic"].values():
        print("-- Basic Load  Name: {:}  Number: {:}  Title: {:}  System: {:}"
              .format(wk.name, wk.number, wk.title, "global"))
        out_rac(wk.items)
        print("")
    if reactions["combination"]:
    #try:
        #comb = pickle.load(file)
        print('{:}'.format(52 * '-'))
        #print("** Static Displacements")
        print("")
        for key, wk in reactions["combination"].items():
            print("-- Load Combination  Name: {:}  Number: {:}  Title: {:}  System: {:}"
                  .format(wk.name, wk.number, wk.title, "global"))
            out_rac(wk.items)
            print("")
    #print('-->')
#
def out_rac(disp, str1=""):
    """
    Report the relative displacement of the structure
    """
    # TODO: send this to mesh node?
    print("     node  x-force    y-force    z-force    x-moment   y-moment   z-moment")
    for item in disp:
        print("{:9d} ".format(item[0]), end='')
        for _disp in item[1:]:
            print("{: 1.3e} ".format(_disp), end='')
        print("")