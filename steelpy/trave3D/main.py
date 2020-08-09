# 
# Copyright (c) 2009-2020 fem2ufo
#

# Python stdlib imports
from array import array
#import datetime
from typing import List, ClassVar, Dict, NamedTuple, Union
from itertools import chain
import pickle
import time

#import multiprocessing

# package imports
#from scipy.linalg import solve_banded
#from steelpy.process.units.units import Units
#from steelpy.trave3D.preprocessor.assemble import assemble_matrix
from steelpy.trave3D.processor.static_solver import (UDUt, memload,
                                                     solve_basic_load, solve_combinations)
from steelpy.trave3D.postprocessor.output import  print_deflections, print_member_forces
from steelpy.trave3D.postprocessor.operations import Results
#
from steelpy.trave3D.preprocessor.assemble import (form_Kmatrix, shape_cond,
                                                   max_bandwidth)
from steelpy.trave3D.processor.operations import to_matrix





#
#
class Trave:
    """
    A program for static & dynamic analysis
    of 3-d framed structures
    """
    __slots__ = ['_mesh', '_load', '_sql']
    
    def __init__(self) -> None:
        """
        """       
        #self._units = Units()
        #self._mesh = Mesh()
        print ( "-- module : trave3D version 2.10" )
        print ( '{:}'.format ( 52 * '-' ) )
    
    
    @property
    def f2u_model(self):
        """
        """
        return self._sql
    
    @f2u_model.setter
    def f2u_model(self, model):
        """
        """
        self._sql = model.sql
        #self._mesh = model.mesh
        self._load = model.load
    
    #
    @property
    def load(self):
        """
        """
        return self._load
    
    @load.setter
    def load(self, value):
        """
        """
        self._load = value
    
    #@property
    #def mesh(self):
    #    """
    #    """
    #    return self._mesh
    
    #@mesh.setter
    #def mesh(self, value):
    #    """
    #    """
    #    self._mesh = value
    #
    #
    def solve_static(self):
        """
        Solves the static system 
        """
        #print("** 3Dframe version 2.10")
        #print ( '{:}'.format ( 52 * '-' ) )
        #print("** Date : {:}".format(datetime.date.today()))
        #
        #self.mesh.load_sql_model()
        #
        #self._banded_stiffness_matrix()
        self.assemble_banded_matrix()
        #
        self.solve_joint_displacement()
        #
        #print_deflections(self.mesh._nodes)
        #
        self.solve_member_forces()
        #
        #
        #print("-->")
    #
    #
    def assemble_banded_matrix(self):
        """
        Asseable the element matrices in upper band form;
        call separatly from formstif, formmass, formgeom 
        -------------------------------------------------
        aa : stiffness matrix
        a  : members rotation matrix
        jbc : nodes freedom
        i1  : node end 1
        j1  : node end 2
        
        aa, a, jbc, i1, j1
        """
        #
        nodes = self._sql.nodes
        boundaries = self._sql.boundaries
        materials = self._sql.materials
        sections = self._sql.sections
        elements = self._sql.elements
        free_nodes = self._sql.get_free_nodes
        #
        print("** Processing Global [K] Matrix")
        jbcc, neq = shape_cond(elements=elements, 
                               nodes=nodes, 
                               boundaries=boundaries, 
                               free_nodes=free_nodes)
        #
        iband = max_bandwidth(elements=elements, 
                              nodes=nodes, jbc=jbcc )
        print("** From datack: half band = {:}".format(iband))
        #
        #
        # ---------------------------
        # multiprocessing
        #
        #call(["python", "steelpy//frame3D//preprocessor//assemblyMatrix.py"])    
        #
        #with Manager() as manager:
        #    d = manager.dict(elements)
        #    l = manager.list(jbc)
        #    st = manager.list(aa)
        #    p = Process(target=loop_members, args=(d, l, st))
        #    p.start()
        #    p.join()
        #
        #
        # ---------------------------
        # Normal 
        #
        jbc = list(chain.from_iterable( jbcc ))
        stf = form_Kmatrix(elements= elements, 
                           nodes=nodes,
                           materials=materials, 
                           sections=sections,
                           jbc=jbc, neq=neq, iband=iband)
        #print("** Finished Processing Global [K] Matrix")
        #
        #print("** Processing Banded [K] Matrix")
        aa = UDUt(stf, neq, iband)
        #
        #
        with open("stfmx.f2u", "wb") as f:
            pickle.dump(neq, f)
            pickle.dump(iband, f)
            pickle.dump(jbcc, f)
            pickle.dump(aa, f)             
        #
        #print("** End Processing Banded [K]")
    #    
    def solve_joint_displacement(self):
        """
        """
        nodes = self._sql.nodes
        materials = self._sql.materials
        sections = self._sql.sections
        elements = self._sql.elements
        #free_nodes = self._mesh.get_free_nodes
        #
        #basic_load = self._mesh.sql.basic_load
        #
        print("** Calculating Joint Displacements")
        print("** reloaded [k] & {p}")
        basic_load = self._load._basic
        basic_res, memf_basic = solve_basic_load(elements=elements, 
                                                 nodes=nodes, materials=materials,
                                                 sections=sections, basic_load=basic_load)
        #
        load_combination = self._load._combination
        comb_res, memf_comb = solve_combinations(basic_res, memf_basic,
                                                 load_combination)
        #
        with open("elemfout.f2u", "wb") as f:
            pickle.dump(memf_basic, f)
            pickle.dump(memf_comb, f)
        #
        #file = open( "elemfout.f2u", "wb" )
        #pickle.dump( memf_basic, file)
        #pickle.dump( memf_comb, file)
        #file.close()
        #
        #self._mesh.sql.element_forces = memf_basic
        #
        # expand displacement to full vector
        print("** reassigning node displacement")
        #
        file = open( "stfmx.f2u", "rb" )
        neq = pickle.load(file )
        iband = pickle.load( file )
        jbc = pickle.load( file )
        #stff = pickle.load( file )
        file.close()    
        #
        jbcc = list(chain.from_iterable(jbc))
        disp_ext = {}
        for key, basic in basic_res.items():
            disp = [basic[ieqnum - 1] if ieqnum != 0 else ieqnum for ieqnum in jbcc]
            disp = to_matrix(disp, 6)
            disp_ext[key] = Results(name=key, title=basic_load[key].name,
                                    load_type="basic", items=disp)
        ## dispp = to_matrix(disp, 6)
        self._sql.displacements = disp_ext
        #
        comb_exp = {}
        for key, comb in comb_res.items():
            disp = [comb[ieqnum - 1] if ieqnum != 0 else ieqnum for ieqnum in jbcc]
            disp = to_matrix(disp, 6)
            comb_exp[key] = Results(name=key, title=load_combination[key].name,
                                    load_type="combination", items=disp)
        if comb_exp:
            self._sql.displacements = comb_exp
        #
        print("** Finished Calculating Joint Displacements")    
    #
    def solve_member_forces(self):
        """
        """
        start_time = time.time()
        #print(" ")
        print("** Calculating Member Forces")
        print("** Reloaded Joint Displacements")
        #
        basic_load = self._load._basic
        load_combination = self._load._combination
        #
        nodes = self._sql.nodes
        materials = self._sql.materials
        sections = self._sql.sections
        elements = self._sql.elements      
        #
        filem = open( "elemfout.f2u", "rb" )
        #
        displacement = self._sql.displacements
        #basic_res = pickle.load(filed)
        mbload = pickle.load(filem)
        memf_basic = {}
        for key, disp in displacement["basic"].items():
            m_load = mbload[key]
            memf_basic[key] = Results(name=key, title=basic_load[key].name,
                                      load_type="basic",
                                      items=memload(elements, nodes, materials, sections,
                                                     disp.items, m_load))
        #
        self._sql.element_forces = memf_basic
        #
        memf_comb = {}
        if displacement["combination"]:
            #comb_res = pickle.load(filed)
            mcload = pickle.load(filem)
            for key, disp in displacement["combination"].items():
                m_load = mcload[key]
                memf_comb[key] = Results(name=key, title=load_combination[key].name,
                                         load_type="combination",
                                         items=memload(elements, nodes, materials, sections,
                                                        disp.items, m_load))
            #
            #if memf_comb:
            self._sql.element_forces = memf_comb
        #
        #
        filem.close()
        #
        end_time = time.time()
        uptime = end_time - start_time
        print("** Calculating Member Forces Process Time: {:1.4e} sec".format(uptime))       
        print("** End Calculating Member Forces")
    #    
    #
    def print_results(self):
        """
        """
        #
        disp_result = self._sql.displacements
        #
        nodes = self._sql.nodes
        print_deflections(nodes, disp_result)
        #
        elements = self._sql.elements 
        member_forces = self._sql.element_forces
        print_member_forces(elements, member_forces)
        #
#    
    
    
    
    