# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
from dataclasses import dataclass
from collections.abc import Mapping
from typing import NamedTuple
import re
#
# package imports
from steelpy.ufo.mesh.main import MeshPlane
from steelpy.utils.dataframe.main import DBframework
from steelpy.utils.sqlite.utils import create_connection #, create_table


#
#
class BeamResBasic(Mapping):
    __slots__ = ['_labels', '_mesh', 'plane']
    
    def __init__(self, mesh) -> None:
        """
        Beam element 
        """
        self._mesh = mesh
        self._labels = list(self._mesh._elements._beams.keys())
        self.plane = mesh._plane # self._plane()
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
    #
    def __str__(self) -> str:
        """ """
        output = "\n"
        output += self.force().__str__()
        output += self.displacement().__str__()
        output += self.stress().__str__()
        return output
    #    
#
#
#
class Beam(BeamResBasic):
    __slots__ = ['_labels', '_mesh', 'plane', 'db_file']
    
    def __init__(self, mesh, db_file: str):
        """ """
        self.db_file = db_file
        super().__init__(mesh)
    #
    #
    def __getitem__(self, beam_name: int|str):
        """ """
        try:
            self._labels.index(beam_name)
            return BeamResItem(beam=self._mesh._elements._beams[beam_name],
                               plane=self.plane, 
                               db_file=self.db_file)
        except ValueError:
            raise IndexError(' ** element {:} does not exist'.format(beam_name))    
    #
    # -----------------------------
    #
    #@property
    def _plane(self):
        """ """
        conn = create_connection(self.db_file)
        with conn:
            sql = 'SELECT * FROM Component'
            cur = conn.cursor()
            cur.execute(sql)
            data = cur.fetchone()
        #
        if data[5] == '2D':
            return MeshPlane(plane2D=True)
        return MeshPlane(plane2D=False)
    #
    # -----------------------------
    #
    #@property
    def force(self, units:str='si'):
        """beam integration forces"""
        #
        conn = create_connection(self.db_file)
        with conn:        
            df = get_force(conn)
        #
        if self.plane.plane2D:
            df.drop(['Fz', 'Mx', 'My'], axis=1, inplace=True)
        #
        return BeamForce(df, units=units)
    #
    #@property
    def displacement(self, units:str='si'):
        """beam integration forces"""
        #
        conn = create_connection(self.db_file)
        with conn:        
            df = get_displacement(conn)
        #
        if self.plane.plane2D:
            df.drop(['z', 'rx', 'ry'], axis=1, inplace=True)        
        #
        return BeamDeflection(df, units=units)  
    #
    def stress(self, units:str='si'):
        """ beam stress """
        #
        conn = create_connection(self.db_file)
        with conn:        
            df = get_stress(conn)
        #
        #if self.plane.plane2D:
        #    df.drop(['tau_z', 'sigma_x', 'sigma_y'], axis=1, inplace=True)        
        #
        return BeamStress(df, units=units)         
    #
    def code_check(self):
        """ beam code check"""
        1 / 0
    #
    #def pull_data(self, query: str, cols: list):
    #    """ """
    #    conn = create_connection(self.db_file)
    #    with conn:
    #        cur = conn.cursor()
    #        cur.execute(query)
    #        data = cur.fetchall()
    #    #       
    #    db = DBframework()
    #    return db.DataFrame(data=data, columns=cols)    
    #
#
#
# --------------------
# Beam Ops
#
#
@dataclass
class BeamResItem:
    __slots__ = ['_beam', '_db_file', '_plane']
    
    def __init__(self, beam, plane, db_file: str):
        """ """
        self._beam = beam
        self._plane = plane
        self._db_file = db_file
    #
    def force(self, units:str='si'):
        """node force"""
        beam_name = self._beam.name
        conn = create_connection(self._db_file)
        with conn:        
            df = get_force(conn, beam_name)
        #
        if self._plane.plane2D:
            df.drop(['Fz', 'Mx', 'My'], axis=1, inplace=True)
        #
        return BeamForce(df, units=units)
    #
    def displacement(self, units:str='si'):
        """ """
        beam_name = self._beam.name
        conn = create_connection(self._db_file)
        with conn:        
            df = get_displacement(conn, beam_name)
        #
        if self._plane.plane2D:
            df.drop(['z', 'rx', 'ry'], axis=1, inplace=True)        
        #        
        return BeamDeflection(df, units=units)
    #
    def stress(self, units:str='si'):
        """node force"""
        beam_name = self._beam.name
        conn = create_connection(self._db_file)
        with conn:        
            df = get_stress(conn, beam_name)
        #
        #if self._plane.plane2D:
        #    df.drop(['Fz', 'Mx', 'My'], axis=1, inplace=True)
        #
        return BeamStress(df, units=units)
    #
#
#
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
        output += print_bitems(self.df)
        return output
#
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
        output += print_bitems(self.df)
        return output
#
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
        output += print_bitems(self.df)
        return output
#
#
#
#
# --------------------
# sql operations
#
def get_force(conn, element_name: int|str|None=None,
              item:str='*'):
    """ """
    if element_name:
        #project = (node_name,)
        query = f'SELECT {item} FROM ResultBeamForce \
                                WHERE element_name = {element_name}'
    else:
        query = f'SELECT {item} FROM ResultBeamForce'
    
    cols = ['number', 'component_name', 
            'load_name', 'load_level', 'load_system', 
            'element_name', 'node_end',
            'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz', 'B', 'Tw']
    df = pull_data(conn, query, cols)
    df.drop(columns=['B', 'Tw'], inplace=True)
    return df
#
#
def get_displacement(conn, element_name: int|str|None=None,
                     item:str='*'):
    """ """
    if element_name:
        #project = (node_name,)
        query = f'SELECT {item} FROM ResultBeamU \
                  WHERE element_name = {element_name}'
    else:
        query = f'SELECT {item} FROM ResultBeamU'
    
    cols = ['number', 'load_name', 'component_name', 
            'load_level', 'load_system', 
            'element_name', 'node_end',
            'x', 'y', 'z', 'rx', 'ry', 'rz']
    
    return pull_data(conn, query, cols)
#
#
def get_stress(conn, element_name: int|str|None=None,
                     item:str='*'):
    """ """
    if element_name:
        #project = (node_name,)
        query = f'SELECT {item} FROM ResultBeamStress \
                  WHERE element_name = {element_name}'
    else:
        query = f'SELECT {item} FROM ResultBeamStress'
    
    cols = ['number', 'load_name', 'component_name', 
            'load_level', 'load_system', 
            'element_name', 'node_end', 'stress_points', 
            'y', 'z', 'tau_x', 'tau_y', 'tau_z','sigma_x', 'sigma_y', 'sigma_z']
    
    return pull_data(conn, query, cols)
#
#
def pull_data(conn, query: str, cols: list):
    """ """
    cur = conn.cursor()
    cur.execute(query)
    data = cur.fetchall()
    #
    db = DBframework()
    return db.DataFrame(data=data, columns=cols)
#
#
# --------------------
# Printing
#
def print_bitems(items,
                 cols: list = ['number', 'load_name', 'component_name',
                               'load_level', 'load_system']):
    """
    """
    items.rename(columns={'element_name': 'beam',
                          'node_end': 'len',
                          'stress_points': 's_points',},
                 inplace=True)
    #
    blgrp = items.groupby('load_level')
    #  
    # basic
    bltype = blgrp.get_group('basic')
    blitems = bltype.groupby(['load_name', 'component_name', 'load_system'])
    #
    output = ''
    # basic
    for key, wk in blitems:
        #header =  wk[['load_name', 'load_system']].values
        #output += "-- Basic Load  Name: {:}  System: {:}".format(*header[0])
        output += "-- Basic Load  Name: {:}  Component: {:} System: {:}\n".format(*key)
        output += "\n"
        #
        vals = wk.drop(cols, axis=1)
        header2 = list(vals)
        #header2 = '       '.join(header2)
        header2 = get_gap(header2)
        #
        output += header2
        output += "\n"        
        output += printout(vals)
        output += "\n"
    #
    # combination
    try:
        bltype = blgrp.get_group('combination')
        blitems = bltype.groupby(['load_name', 'component_name', 'load_system'])
        #
        output += '{:}\n'.format(52 * '-')
        for key, wk in blitems:
            #header =  wk[['load_name', 'load_system']].values
            #output += "-- Load Combination  Name: {:}  System: {:}".format(*header[0])
            output += "-- Load Combination  Name: {:} Component: {:} System: {:}\n".format(*key)
            output += "\n"
            #
            vals = wk.drop(cols, axis=1)
            header2 = list(vals)
            #header2 = '       '.join(header2)
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
#
def printout(bforces):
    """
    """
    #
    output = ""
    #
    mgroup = bforces.groupby("beam")
    #bforces.set_index('beam', inplace=True)
    #
    for key, mgroup in mgroup:
    #for item in bforces.itertuples():
        output += "{:9d} ".format(key)
        items = mgroup.set_index('beam')
        #for x, member in enumerate(items):
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