# 
# Copyright (c) 2009 steelpy
# 

# Python stdlib imports
from __future__ import annotations
#from array import array
#from dataclasses import dataclass
#from collections.abc import Mapping
#from typing import NamedTuple
#from collections import defaultdict
#import re

# package imports
#
from steelpy.ufo.load.concept.beam import BeamLoadItemIM
from steelpy.ufo.load.concept.combination import LoadCombConcept
from steelpy.sections.main import SectionIM
from steelpy.sections.utils.operations import get_section
from steelpy.material.main import MaterialIM
from steelpy.trave.beam.operation import BeamBasic
from steelpy.utils.dataframe.main import DBframework
#
from steelpy.material.utils.operations import get_mat_properties
from steelpy.ufo.load.process.combination import (get_comb_dict,
                                                  get_comb_list)
# 
#
class Beam:
    __slots__ = ['_support', '_basic', '_sections', '_materials',
                 'steps', '_length', '_response', '_Pdelta', 
                 '_combination', '_design', 'name', '_db']

    def __init__(self, name:int|str, steps: int = 10):
        """
        beam_length : total length of the beam 
        section_type:str
        """
        print("-- module : Beam Version 0.25dev")
        print('{:}'.format(52*'-'))
        #
        self.name = name
        self.steps: int = steps
        #
        self._materials = MaterialIM()
        self._sections = SectionIM()
        #
        self._basic = BeamLoadItemIM(load_name="beam", 
                                     load_title="beam",
                                     component=self.name)
        self._basic._load.local_system()
        #
        self._combination = LoadCombConcept(basic_load=self._basic,
                                            component=self.name)
        #
        self._db = DBframework()
        self._Pdelta = True
    #
    #
    #
    @property
    def L(self):
        """ """
        return self._length
    
    @L.setter
    def L(self, length):
        """ """
        try:
            self._length = length.value
        except AttributeError: # unit must be in metres
            self._length = length
        #
        
    #
    # -----------------------------------------------------
    #    
    @property
    def section(self):
        """
        """
        return self._sections['beam_sect']

    @section.setter
    def section(self, properties:list|tuple|dict):
        """
        section_type : select beam's cross-section shape (tubular/rectangular,I_section)
        """
        properties = get_section(properties)
        self._sections["beam_sect"] = properties
    #
    @property
    def material(self):
        """
        """
        return self._materials['beam_mat']

    @material.setter
    def material(self, properties:list|tuple|dict):
        """
        """
        # [elastic, Fy, Fu, E, G, poisson, density, alpha]
        properties = get_mat_properties(properties)
        self._materials['beam_mat'] = properties
    #
    # -----------------------------------------------------
    #
    @property
    def support(self):
        """
        """
        return self._support

    @support.setter
    def support(self, value:list):
        """
        [end1, end2, torsion1, torsion2, k1, k2]
        """
        support = [item for item in value[:2]]
        try:
            support.append(value[2]) # t1
        except:
            support.append(value[0])
        
        try:
            support.append(value[3]) # t2
        except:
            support.append(value[1])
        
        try:
            support.append(value[4]) # k1
        except:
            support.append(None)

        try:
            support.append(value[5]) # k1
        except:
            support.append(None)

        self._support = support
    #
    # -----------------------------------------------------
    # Loading
    #
    @property
    def selfweight(self):
        """ """
        return self._basic._load._selfweight
    
    @selfweight.setter
    def selfweight(self, value:list|dict):
        """ """
        self._basic._load._selfweight = value
    #
    #
    @property
    def P(self):
        """
        """
        return self._basic[self.name].point

    @P.setter
    def P(self, value:list|dict):
        """
        """
        self._basic[self.name].point = value
    #
    @property
    def q(self):
        """line load"""
        return self._basic[self.name].line

    @q.setter
    def q(self, value:list|tuple|dict):
        """line load"""
        self._basic[self.name].line = value
    #
    @property
    def load_combination(self):
        """load combination"""
        return self._combination

    @load_combination.setter
    def load_combination(self, values:list|dict):
        """load combination"""
        if isinstance(values, dict):
            comb = get_comb_dict(values)
            sname = comb[2]  # lname
            if isinstance(sname, (list | tuple)):
                for x, lname in enumerate(sname):
                    cname = self._check_type(comb[0], step=x)
                    ltype = self._check_type(comb[1], step=x)
                    factor = self._check_type(comb[3], step=x)
                    title = self._check_type(comb[4], step=x)
                    self._comb_setup(values=[cname, ltype, lname, factor, title])
            else:
                self._comb_setup(values=comb)
        else:
            comb = get_comb_list(values)
            self._comb_setup(values=comb)
    #
    # ----------------------------
    #
    def _check_type(self, item:list|str|float|int, step:int):
        """ """
        if isinstance(item, list):
            return item[step]
        else:
            return item
    #
    def _comb_setup(self, values: list):
        """" values : [name, type, load_id, factor, title] """
        cname = values[0]
        ltype = values[1]
        lname = values[2]
        factor = float(values[3])
        try:
            title = values[4]
        except IndexError:
            title = cname
        #
        try:
            self._combination[cname]
        except KeyError:
            self._combination[cname] = title
        #
        match ltype:
            case 'basic':
                self._combination[cname]._basic[lname] = factor
            case 'combination':
                self._combination[cname]._combination[lname] = factor
            case _:
                raise IOError(f"Combination load type {ltype} not available")
    #
    # -----------------------------------------------------
    # TODO: beam with multiple supports
    #
    def join(self, other, join_type):
        """
        b1 = Beam()
        b2 = Beam()
        b = b1.join(b2, 'fixed')
        """
        pass
    #
    #
    # -----------------------------------------------------
    # Calculations
    #
    #    
    #
    def _beam(self, Pdelta: bool):
        """get beam""" 
        material = self.material
        section = self.section.properties(poisson=material.poisson)
        #
        return BeamBasic(L=self._length, area=section.area, 
                         Iy=section.Iy, Iz=section.Iz,
                         J=section.J, Cw=section.Cw, 
                         E=material.E, G=material.G,
                         Asy=section.Asy,
                         Asz=section.Asz,
                         Pdelta=Pdelta)
    #
    def _get_loads(self)-> list:
        """get beam loading"""
        bloads:list = []
        # basic load
        for item in self._basic._load._line.values():
            bloads.extend(item)
        #
        for item in self._basic._load._point.values():
            bloads.extend(item)
        # combination
        if self._combination:
            comb = self._combination.to_basic()
            #comb
            #for items in self._combination.values():
            #    # basic
            #    for item in items.basic.values():
            #        bloads.extend(item)
            #    # combination
            #    for item in items.combination.values():
            #        bloads.extend(item)

        return bloads
    #
    #    
    #
    # -----------------------------------------------------
    # Results    
    #
    #
    #@property
    #def shear(self):
    #    """ """
    #    return self._response.shear(self.steps)
    #
    #@property
    #def bending_moment(self):
    #    """ """
    #    return self._response.bending(self.steps)
    #
    #@property
    #def slope(self):
    #    """ """
    #    return self._response.slope(self.steps)
    #
    #@property
    #def deflection(self):
    #    """ """
    #    return self._response.deflection(self.steps)
    #
    def reactions(self):
        """
        Return: 
        dataframe : [load_name, load_title, load_type, system,
                     element_name, support_type, support_end,
                     Fx, Fy, Fz, Mx, My, Mz,
                     delta_x, delta_y, delta_z,
                     phi_x, theta_y, theta_z]
        """
        bloads = self._get_loads()
        #
        beam = self._beam(Pdelta=self._Pdelta)
        beam.supports(*self.support)
        #
        # load_name : [load_title, load_type,
        #              load_level, load_system, 
        #              beam_name,
        #              R0[Fa, Tx, Ry, Rz],
        #              R1[Fa, Tx, Ry, Rz]]
        #
        # Tx = [FT, FB, Fpsi, Fphi, Tw]
        # Fa, Ry & Rz = [V, M, theta, w]        
        reactions = beam.reactions(bloads)
        #
        # Organize data for reaction output
        msupport = []
        for key, item in reactions.items():
            for reac in item:
                msupport.extend([[*reac[:7],            # load details & beam name
                                  self._support[x], x+1, # support type & beam end
                                  lbf[0][0],  # Fx [axial]
                                  lbf[1][0],  # Mx [torsion]
                                  lbf[1][3],  # phi_x [torsion]
                                  lbf[0][3],  # delta_x [axial]
                                  *lbf[2],    # Bending in_plane
                                  *lbf[3]]    # Bending out_plane
                                 for x, lbf in enumerate(reac[7:])]) # R0 & R1
        #
        # Dataframe setup 
        #
        header = ['load_name', 'mesh_name',
                  'load_title', 'load_type',
                  'load_level', 'load_system',
                  'element_name', 
                  'support_type', 'support_end',
                  'Fx', 'Mx', 'theta_x', 'delta_x',
                  'Fy', 'Mz', 'theta_z', 'delta_y',
                  'Fz', 'My', 'theta_y', 'delta_z']
        #
        df_mload = self._db.DataFrame(data=msupport, columns=header, index=None)
        df_mload = df_mload[['load_name', 'mesh_name',
                             'load_title', 'load_type',
                             'load_level', 'load_system',
                             'element_name', 'support_type', 'support_end',
                             'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz',
                             'delta_x', 'delta_y', 'delta_z',
                             'theta_x', 'theta_y', 'theta_z']]
        #
        return df_mload    
    #
    #
    #@property
    def response(self):
        """
        Results:
        Beam force response 
        Dataframe [load_name, load_title, load_type, load_system,
                   element_name, node_end,
                   Fx, Fy, Fz, Mx, My, Mz,
                   delta_x, delta_y, delta_z,
                   phi_x, theta_y, theta_z,
                   psi_t (phi_x'), B(bimoment), Tw(Warping torque)]
        """
        bloads = self._get_loads()
        #        
        # -----------------------------------------------------
        #        
        beam = self._beam(Pdelta=self._Pdelta)
        beam.supports(*self.support)
        df_mload = beam.forces(bloads, steps=self.steps)
        return df_mload
    #
    #
    def stress(self):
        """calculate beam stress"""
        mat = self.material
        #
        actions_df = self.response(steps=self.steps)
        members = actions_df.groupby(['element_name',
                                      'load_name', 'component_name',
                                      'load_title', 
                                      'load_level', 'load_system'])
        #
        frames = []
        for key, item in members:
            frames.append(self.section.stress(df=item,
                                              E=mat.E, G=mat.G))
        #
        try:
            stress = DBframework.concat(frames, ignore_index=True)
        except TypeError:
            stress = frames[0]
        #1 / 0
        return stress
    #
    @property
    def design(self, material=None):
        """
        """
        self._design.material = material
        if not material:
            self._design.material = self._material
        #
        self._design.section = self._section
        self._design.L = self.length
        self._design._stress = self._response.stress()
        return self._design
#
#
#
#
