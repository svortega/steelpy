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
#
from steelpy.trave3D.processor.static_solver import UDUt, beam_force, solve_basic_load
#
from steelpy.trave3D.preprocessor.assemble import form_Kmatrix, get_bandwidth
from steelpy.trave3D.processor.operations import to_matrix





#
#
class Trave3D:
    """
    A program for static & dynamic analysis
    of 3-d framed structures
    """
    __slots__ = ['_mesh', '_load', '_f2u']
    
    def __init__(self) -> None:
        """
        """       
        #self._units = Units()
        #self._mesh = Mesh()
        print ("-- module : trave3D version 2.10")
        print ('{:}'.format(52 * '-'))
    
    
    @property
    def f2u_model(self):
        """
        """
        return self._f2u
    
    @f2u_model.setter
    def f2u_model(self, model):
        """
        """
        self._f2u = model
        #self._mesh = model.mesh
        #self._load = model.load
    
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
    def run_static(self):
        """
        Solves the static system 
        """
        #self._banded_stiffness_matrix()
        self.assemble_banded_matrix()
        self.solve_joint_displacement()
        self.solve_member_forces()
        #print("-->")
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
        materials = self._f2u.materials
        sections = self._f2u.sections
        elements = self._f2u.mesh.elements
        nodes = self._f2u.mesh.nodes
        boundaries = self._f2u.mesh.boundaries
        #
        free_nodes = elements.get_free_nodes
        #
        print("** Processing Global [K] Matrix")
        #
        jbcc, neq, iband = get_bandwidth(elements=elements, 
                                         nodes=nodes, 
                                         boundaries=boundaries, 
                                         free_nodes=free_nodes)
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
        jbc = list(chain.from_iterable(jbcc))
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
    #
    def solve_joint_displacement(self):
        """
        """
        materials = self._f2u.materials
        sections = self._f2u.sections
        elements = self._f2u.mesh.elements
        nodes = self._f2u.mesh.nodes
        #
        print("** Calculating Joint Displacements")
        print("** reloaded [k] & {p}")
        basic_load = self._f2u.load.basic
        basic_res, memf_basic = solve_basic_load(elements=elements, 
                                                 nodes=nodes, materials=materials,
                                                 sections=sections, basic_load=basic_load)
        #
        load_combination = self._f2u.load.combination
        comb_res, memf_comb = load_combination.solve_combinations(basic_res, memf_basic)
        #
        with open("elemfout.f2u", "wb") as f:
            pickle.dump(memf_basic, f)
            pickle.dump(memf_comb, f)
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
        # basic load
        jbcc = list(chain.from_iterable(jbc))
        # TODO node name instead number
        node_index = {item.index:key for key, item in nodes.items()}
        for load_tile, basic in basic_res.items():
            disp = [basic[ieqnum - 1] if ieqnum != 0 else ieqnum
                    for ieqnum in jbcc]
            disp = to_matrix(disp, 6)
            disp = [[node_index[number], load_tile, "global"] + item
                    for number, item in enumerate(disp)]
            self._f2u._results.node.deflection = disp
        # load combination
        for load_tile, comb in comb_res.items():
            disp = [comb[ieqnum - 1] if ieqnum != 0 else ieqnum
                    for ieqnum in jbcc]
            disp = to_matrix(disp, 6)
            disp = [[node_index[number], load_tile, "global"] + item
                    for number, item in enumerate(disp)]
            self._f2u._results.node.deflection = disp
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
        elements = self._f2u.mesh.elements
        basic_load = self._f2u.load.basic
        load_combination = self._f2u.load.combination
        displacement = self._f2u._results.node.deflection
        #
        filem = open("elemfout.f2u", "rb")
        mbload = pickle.load(filem)
        memf_basic = []
        for key, disp in displacement["basic"].items():
            m_load = mbload[key]
            # select [x, y, z, rx, ry, rz] with list in node index order
            ndisp = [item[1:] for item in disp.items]
            force = beam_force(elements, ndisp, m_load)
            memf_basic.extend([[item[0], key, item[2], 0, item[1], *item[3:9]]
                               for item in force])
            memf_basic.extend([[item[0], key, item[9], 1, item[1], *item[10:]]
                               for item in force])
        self._f2u._results.element.force = memf_basic
        #
        if displacement["combination"]:
            memf_comb = []
            mcload = pickle.load(filem)
            for key, disp in displacement["combination"].items():
                m_load = mcload[key]
                ndisp = [item[1:] for item in disp.items]
                force = beam_force(elements, ndisp, m_load)
                memf_comb.extend([[item[0], key, item[2], 0, item[1], *item[3:9]]
                                  for item in force])
                memf_comb.extend([[item[0], key, item[9], 1, item[1], *item[10:]]
                                  for item in force])
            self._f2u._results.element.force = memf_comb
        #
        filem.close()
        end_time = time.time()
        uptime = end_time - start_time
        print("** Calculating Member Forces Process Time: {:1.4e} sec".format(uptime))       
        print("** End Calculating Member Forces")
    #
    def print_results(self):
        """
        """
        self._f2u._results.node.print_deflections()
        self._f2u._results.element.print_forces()
    #
    #
#
#    
    
    
    
    