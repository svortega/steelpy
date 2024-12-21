# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
#from dataclasses import dataclass
from collections.abc import Mapping
from typing import NamedTuple
import re
#
# package imports
from steelpy.utils.dataframe.main import DBframework


#
#
class BeamResBasic(Mapping):
    __slots__ = ['_mesh', '_beams']
    
    def __init__(self, mesh) -> None:
        """
        Beam element 
        """
        self._mesh = mesh
        #self._labels = list(self._mesh._elements._beams.keys())
        #self.plane = mesh._plane # self._plane()
        self._beams = mesh._elements._beams
    #
    #
    def __contains__(self, value) -> bool:
        return value in self._labels
    
    def __iter__(self):
        """
        """
        return iter(self._labels)

    def __len__(self) -> float:
        return len(self._labels)
    
    def __str__(self) -> str:
        """ """
        output = "\n"
        output += self.force().__str__()
        output += self.displacement().__str__()
        output += self.stress().__str__()
        return output
#
#
# --------------------
# Beam Ops
#
class BeamForce(NamedTuple):
    """ Basic load transfer"""
    df: DBframework
    units: str
    #
    def __str__(self):
        """ """
        unitsout = "Units : SI [N/N-m]"
        if re.match(r"\b(us|imperial)\b", self.units, re.IGNORECASE):
            unitsout = "units : US [lb/lb-ft]"
        #
        output = "\n"
        output += '{:}\n'.format(52 * '-')
        output += "** Beam Forces | "
        output += unitsout
        output += "\n"
        output += '{:}\n'.format(52 * '-')
        output += print_beam_items(self.df)
        return output
#
class BeamDeflection(NamedTuple):
    """ Basic load transfer"""
    df: DBframework
    units: str
    #
    def __str__(self):
        """ """
        unitsout = "Units : SI [m/radians]"
        if re.match(r"\b(us|imperial)\b", self.units, re.IGNORECASE):
            unitsout = "units : US [ft/radians]"
        #        
        output = "\n"
        output += '{:}\n'.format(52 * '-')
        output += "** Beam Displacement | "
        output += unitsout
        output += "\n"
        output += '{:}\n'.format(52 * '-')
        output += print_beam_items(self.df)
        return output
#
class BeamStress(NamedTuple):
    """ Basic load transfer"""
    df: DBframework
    units: str
    #
    def __str__(self):
        """ """
        unitsout = "Units : SI [Pa]"
        if re.match(r"\b(us|imperial)\b", self.units, re.IGNORECASE):
            unitsout = "units : US [psi]"
        #         
        output = "\n"
        output += '{:}\n'.format(52 * '-')
        output += "** Beam Stress | "
        output += unitsout
        output += "\n"        
        output += '{:}\n'.format(52 * '-')
        output += print_beam_items(self.df)
        return output
#
# --------------------
# Printing
#
def print_beam_items(items,
                 cols: list|None = None):
    """
    """
    if not cols:
        cols = ['number', 'load_name', 'component_name',
                'load_level', 'system']
    #
    items.rename(columns={'element_name': 'beam',
                          'length': 'len',
                          'stress_points': 's_points',},
                 inplace=True)
    blgrp = items.groupby('load_level')
    #
    output = ''
    try:
        bltype = blgrp.get_group('basic')
        blitems = bltype.groupby(['load_name', 'component_name', 'system'])
        for key, wk in blitems:
            output += "-- Basic Load  Name: {:}  Component: {:} System: {:}\n".format(*key)
            output += "\n"
            #
            vals = wk.drop(cols, axis=1)
            header2 = list(vals)
            header2 = get_gap(header2)
            #
            output += header2
            output += "\n"        
            output += printout(vals)
            output += "\n"
    except KeyError:
        pass    
    #
    # combination
    try:
        bltype = blgrp.get_group('combination')
        blitems = bltype.groupby(['load_name', 'component_name', 'system'])
        #
        output += '{:}\n'.format(52 * '-')
        for key, wk in blitems:
            output += "-- Load Combination  Name: {:} Component: {:} System: {:}\n".format(*key)
            output += "\n"
            #
            vals = wk.drop(cols, axis=1)
            header2 = list(vals)
            header2 = get_gap(header2)
            #
            output += header2
            output += "\n"               
            output += printout(vals)
            output += "\n"
    except KeyError:
        pass
    #
    return output
#
def printout(bforces):
    """
    """
    output = ""
    mgroup = bforces.groupby("beam")
    #
    for key, mgroup in mgroup:
        output += "{:9d} ".format(key)
        items = mgroup.set_index('beam')
        for x, item in enumerate(items.itertuples()):
            try:
                1 / x
                output += "{:} ".format(9 * " ")
            except ZeroDivisionError:
                pass
            
            output += "{: 1.3e} ".format(item.len)
            for val in item[2:]:
                output += "{: 1.3e} ".format(val)
            output += "\n"
    #
    return output
#
def get_gap(header, step:int=9):
    """ """
    new = []
    for item in header:
        gap = step - len(item) 
        new.extend([item, " " * gap])
    #
    new = " ".join(new)
    return new
#
#