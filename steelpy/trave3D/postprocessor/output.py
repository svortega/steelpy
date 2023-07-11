# Copyright (c) 2009-2023 fem2ufo
#
from __future__ import annotations
# Python stdlib imports
from typing import NamedTuple
#
# package imports
#from steelpy.f2uModel.load.operations.operations import duplicates, indices
#
#
#
class Results(NamedTuple):
    """ Basic load transfer"""
    name: int
    number: int
    title: str
    load_type: str
    items: list[list[float]]
# -------------------------
#
#
def print_node_deflections(disp_result, m2D: bool):
    """
    """
    disphead = ['x', 'y', 'z', 'rx', 'ry', 'rz']
    if m2D:
        disphead = ['x', 'y', 'rz']
    #
    #ndgrp = disp_result.groupby(['load_type', 'load_name'])
    ndgrp = disp_result.groupby('load_type')
    #
    print('{:}'.format(52 * '-'))
    print("** Reloaded Joint Displacements")
    print('{:}'.format(52 * '-'))
    print("")
    #
    ndtype = ndgrp.get_group('basic')
    nditems = ndtype.groupby('load_name')
    # Basic
    for key, wk in nditems:
        header =  wk[['load_name', 'load_number', 'load_title', 'load_system']].values
        print("-- Basic Load  Name: {:}  Number: {:}  Title: {:}  System: {:}"
              .format(*header[0]))
        if m2D:
            print("     node     x [ m]     y [ m]   rz [rad]")
        else:
            print("     node     x [ m]     y [ m]     z [ m]   rx [rad]   ry [rad]   rz [rad]")
        printout(wk.node_name, wk[disphead].values)
        print("")
    #
    # Combination
    try:
        ndtype = ndgrp.get_group('combination')
        nditems = ndtype.groupby('load_name')
        #
        #nditems2 = (ndtype.groupby(['load_name', 'load_number', 'load_title', 'system', 'node_name'])
        #           [['x', 'y', 'z', 'rx', 'ry', 'rz']].sum())
        #
        #nditems2 = nditems2.groupby('load_name')
        #
        #if nditems.columns:
        print('{:}'.format(52 * '-'))
        print("")
        for key, wk in nditems:
            header =  wk[['load_name', 'load_number', 'load_title', 'load_system']].values
            print("-- Load Combination  Name: {:}  Number: {:}  Title: {:}  System: {:}"
                  .format(*header[0]))
            #
            ndsum =  wk.groupby('node_name')[disphead].sum()
            if m2D:
                print("     node     x [ m]     y [ m]   rz [rad]")
            else:
                print("     node     x [ m]     y [ m]     z [ m]   rx [rad]   ry [rad]   rz [rad]")
            printout(ndsum.index, ndsum.values)
            print("")
    except KeyError:
        pass
    #
    #print('-->')
#
def print_node_reactions(reactions, m2D: bool):
    """
    :return:
    """
    forcehead = ['Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']
    if m2D:
        forcehead = ['Fx', 'Fy', 'Mz']
    #
    nreacgrp =  reactions.groupby('load_type')
    #1 / 0
    #
    print('{:}'.format(52 * '-'))
    print("** Reloaded Joint Reactions")
    print('{:}'.format(52 * '-'))
    print("")
    # Basic
    nrtype = nreacgrp.get_group('basic')
    nritems = nrtype.groupby('load_name')
    for key, wk in nritems:
        header =  wk[['load_name', 'load_number', 'load_title', 'load_system']].values
        print("-- Basic Load  Name: {:}  Number: {:}  Title: {:}  System: {:}"
              .format(*header[0]))
        #
        if m2D:
            print("     node    Fx [ N]    Fy [ N]  Mz[ N-m ]")
        else:
            print("     node    Fx [ N]    Fy [ N]    Fz [ N]  Mx[ N-m ]  My[ N-m ]  Mz[ N-m ]")
        printout(wk.node_name, wk[forcehead].values)
        print("")
    #
    # Combinations
    try:
        nrtype = nreacgrp.get_group('combination')
        nritems = nrtype.groupby('load_name')
        #if nritems.columns:
        print('{:}'.format(52 * '-'))
        print("")
        for key, wk in nritems:
            header =  wk[['load_name', 'load_number', 'load_title', 'load_system']].values
            print("-- Load Combination  Name: {:}  Number: {:}  Title: {:}  System: {:}"
                  .format(*header[0]))
            #
            if m2D:
                print("     node    Fx [ N]    Fy [ N]  Mz[ N-m ]")
            else:
                print("     node    Fx [ N]    Fy [ N]    Fz [ N]  Mx[ N-m ]  My[ N-m ]  Mz[ N-m ]")
            printout(wk.node_name, wk[forcehead].values)
            print("")
    except KeyError:
        pass
    #print('-->')
#
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
def printout(node_id:list, disp:list, str1=""):
    """
    Report the relative displacement of the structure
    """
    # TODO: send this to mesh node?
    for x, item in enumerate(node_id):
        print("{:9d} ".format(item), end='')
        for _disp in disp[x]:
            print("{: 1.3e} ".format(_disp), end='')
        print("")
    #
    #get_max_displacement(dispp)
    #print("-->")
#
#
# --------------------
# beam print out
# --------------------
#
def print_beam_forces(df_membforce, m2D: bool):
    """
    """
    blgrp = df_membforce.groupby('load_type')
    #
    print('{:}'.format(52 * '-'))
    print("** Reloaded Elements Forces")
    print('{:}'.format(52 * '-'))
    print("")
    #   
    # basic
    bltype = blgrp.get_group('basic')
    blitems = bltype.groupby('load_name')     
    for key, wk in blitems:
        header =  wk[['load_name', 'load_number', 'load_title', 'load_system']].values
        print("-- Basic Load  Name: {:}  Number: {:}  Title: {:}  System: {:}"
              .format(*header[0]))
        #
        print_beam_int_forces(wk, m2D)
        print("")
    #
    # combinations    
    try:
        bltype = blgrp.get_group('combination')
        blitems = bltype.groupby('load_name')        
        print('{:}'.format(52 * '-'))
        for key, wk in blitems:
            header =  wk[['load_name', 'load_number', 'load_title', 'load_system']].values
            print("-- Load Combination  Name: {:}  Number: {:}  Title: {:}  System: {:}"
                  .format(*header[0]))
            #
            print_beam_int_forces(wk, m2D)
            print("")
    except KeyError:
        pass
#
#
def print_beam_int_forces(eforces, m2D: bool):
    """
    """
    mgroup = eforces.groupby("element_name") #.groups .get_group([0,1])
    #
    if m2D:
        ndof = 3
        forcehead = ['node_end', 'F_Vx', 'F_Vy', 'F_Mz']
        print("     Beam  Len[ m]    Fx [ N]    Fy [ N]  Mz[ N-m ]")
    else:
        ndof = 6
        forcehead = ['node_end', 'F_Vx', 'F_Vy', 'F_Vz', 'F_Mx', 'F_My', 'F_Mz']
        print("     Beam  Len[ m]    Fx [ N]    Fy [ N]    Fz [ N]  Mx[ N-m ]  My[ N-m ]  Mz[ N-m ]")
    #
    for key, mgroup in mgroup:
        print("{:9d} ".format(key), end='')
        items = mgroup[forcehead].values
        
        for x, member in enumerate(items):
            #
            try:
                1 / x
                print("{:} ".format(9 * " "), end='')
            except ZeroDivisionError:
                pass
            
            print("{:1.2e} ".format(member[0]), end='')
            for i in range(ndof):
                print("{: 1.3e} ".format(member[i+1]), end='')            
            #
            #member = eforces.iloc[index]
            #try:
            #    1/member.end
            #    print("{:} ".format(9 * " "), end='')
            #except ZeroDivisionError:
            #    #print('---=======')
            #    pass
            #print("{:10d} ".format(member.node), end='')
            #for i in range(6):
            #    print("{: 1.3e} ".format(member[i+4]), end='')
            #print("{:}".format(8*" "))
            print("")
        #print("\n")
    #print("\n")
#
#
def print_beam_disp(df_membforce, m2D: bool):
    """ """
    blgrp = df_membforce.groupby('load_type')
    #
    print('{:}'.format(52 * '-'))
    print("** Reloaded Elements Displacement")
    print('{:}'.format(52 * '-'))
    print("")
    #   
    # basic
    bltype = blgrp.get_group('basic')
    blitems = bltype.groupby('load_name')     
    for key, wk in blitems:
        header =  wk[['load_name', 'load_number', 'load_title', 'load_system']].values
        print("-- Basic Load  Name: {:}  Number: {:}  Title: {:}  System: {:}"
              .format(*header[0]))
        #
        print_beam_int_disp(wk, m2D)
        print("")
    #
    # combinations    
    try:
        bltype = blgrp.get_group('combination')
        blitems = bltype.groupby('load_name')        
        print('{:}'.format(52 * '-'))
        for key, wk in blitems:
            header =  wk[['load_name', 'load_number', 'load_title', 'load_system']].values
            print("-- Load Combination  Name: {:}  Number: {:}  Title: {:}  System: {:}"
                  .format(*header[0]))
            #
            print_beam_int_disp(wk, m2D)
            print("")
    except KeyError:
        pass    
#
#
def print_beam_int_disp(eforces, m2D: bool):
    """
    """
    mgroup = eforces.groupby("element_name") #.groups .get_group([0,1])
    #
    if m2D:
        ndof = 3
        forcehead = ['node_end', 'F_wx', 'F_wy', 'F_thetaz']
        print("     Beam  Len[ m]     x [ m]     y [ m]   rz [rad]")
    else:
        ndof = 6
        forcehead = ['node_end', 'F_wx', 'F_wy', 'F_wz', 'F_phix', 'F_thetay', 'F_thetaz']
        print("     Beam  Len[ m]     x [ m]     y [ m]     rx [rad]   ry [rad]   rz [rad]")
    #
    for key, mgroup in mgroup:
        print("{:9d} ".format(key), end='')
        items = mgroup[forcehead].values
        
        for x, member in enumerate(items):
            #
            try:
                1 / x
                print("{:} ".format(9 * " "), end='')
            except ZeroDivisionError:
                pass
            
            print("{:1.2e} ".format(member[0]), end='')
            for i in range(ndof):
                print("{: 1.3e} ".format(member[i+1]), end='')            
            #
            #member = eforces.iloc[index]
            #try:
            #    1/member.end
            #    print("{:} ".format(9 * " "), end='')
            #except ZeroDivisionError:
            #    #print('---=======')
            #    pass
            #print("{:10d} ".format(member.node), end='')
            #for i in range(6):
            #    print("{: 1.3e} ".format(member[i+4]), end='')
            #print("{:}".format(8*" "))
            print("")
        #print("\n")
    #print("\n")
#
#
# --------------------
# output
# --------------------
#
#
class ResultInmemory:
    
    __slots__ = ['_displacement', '_reaction', 
                 '_node_force', '_beam_force', '_m2D']

    def __init__(self, m2D: bool):
        """
        """
        #self._basic_loads = basic_loads
        #self._displacement = {}
        #self._beam_force = BeamForce()
        #self._reaction = BoundaryNodes()
        self._m2D = m2D
    #
    #
    @property
    def node_displacement(self):
        """
        [node, load_name, system, x, y, z, rx, ry, rz]
        """
        return self._displacement
    #
    #@node_displacement.setter
    #def node_displacement(self, values:List):
    #    """
    #    values = [node, load_name, system, x, y, z, rx, ry, rz]
    #    """
    #    self._displacement = values
    #    #print("-->")
    #
    def print_node_displacement(self):
        """
        """
        print_node_deflections(self._displacement, self._m2D)
    #
    @property
    def beam_force(self):
        """
        [element, load_name, node_name, position, system, fx, fy, fz, mx, my, mz]
        """
        return self._beam_force

    #@element_force.setter
    #def element_force(self, values:List):
    #    """
    #    values = [element, load_name, node_name, position, system, fx, fy, fz, mx, my, mz]
    #    """
    #    self._beam_force = values
    #    #print("-->")
    #
    def print_element_forces(self):
        """ """
        print_beam_forces(self.beam_force, self._m2D)
    #
    def print_element_disp(self):
        """ """
        print_beam_disp(self.beam_force, self._m2D)
    #
    #
    @property
    def node_reaction(self):
        """
        [node, load_name, system, x, y, z, rx, ry, rz]
        """
        return self._reaction
    #
    #@node_reaction.setter
    #def node_reaction(self, values:List):
    #    """
    #    values = [node, load_name, system, x, y, z, rx, ry, rz]
    #    """
    #    self._reaction = values
    #    #print("-->")
    #
    def print_node_reactions(self):
        """
        """
        print_node_reactions(self._reaction, self._m2D)
    #
#