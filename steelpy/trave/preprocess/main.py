#
# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
from datetime import datetime as dt

# package imports
#from steelpy.trave.process.dynamic import eigen, trnsient
#from steelpy.trave.utils.solution import UnSolver
from steelpy.trave.process.static.main import StaticSolver
#from steelpy.trave.postprocess.main import PostProcess
#
from steelpy.trave.preprocess.sql.main import TraveItemSQL
#from steelpy.utils.dataframe.main import DBframework
#from steelpy.utils.sqlite.utils import create_connection, create_table
#
class TraveItemBasic(TraveItemSQL):
    """
    A program for static & dynamic analysis
    of 3-d framed structures
    """
    __slots__ = ['_name', '_mesh', 'db_file',
                 '_log', '_plane2D'] # 'Un'

    def __init__(self, mesh,
                 #name: str |None = None,
                 sql_file: str |None = None,
                 log: bool = False) -> None:
        """
        """
        # set system
        plane = "3D"
        if self._plane2D:
            plane = "2D"
        #
        # -----------------------------
        print('{:}'.format(52 * '-'))
        print (f"-- module : Trave{plane} version 2.60")
        print ('{:}'.format(52 * '-'))
        #
        # -----------------------------
        self._mesh = mesh
        self._mesh.plane(self._plane2D)
        self._log = log
        #
        if not  sql_file:
            sql_file = self._mesh.db_file
        self.db_file = sql_file
        #
        # -----------------------------
        super().__init__(mesh=self._mesh,
                         #name=self._name,
                         sql_file=self.db_file,
                         log=self._log)
    #
    #
    # --------------------------------------------
    #
    def static(self, name: str |None = None,
               second_order: bool = False,
               inelastic: bool = False):
        """
        Solves the static system by the Direct Stiffness Method (DSM)

        method : banded, frontal
        second_order : Second order (True/False)
        """
        #
        if not name:
            name = self._mesh._name
        self._name = name
        # TODO : check if name is more suitable than id
        result_id = self._push_analysis(name=self._name,
                                        analysis_type='static',
                                        plane=self._plane2D,
                                        Pdelta=second_order,
                                        ineleastic=inelastic)
        #
        if  self._mesh:
            return  StaticSolver(mesh=self._mesh,
                                 result_name=self._name,
                                 db_file=self.db_file,
                                 log = self._log,
                                 second_order=second_order,
                                 nonlinear=inelastic)
        else:
            raise IOError('** error: mesh missing')
    #
    #
    def modal(self, name: str |None = None,
              second_order: bool = False,
              ineleastic: bool = False):
        """ Natural Period"""
        if not name:
            name = self._mesh._name
        self._name = name
        #
        result_id = self._push_analysis(name=self._name,
                                        analysis_type='static',
                                        plane=self._plane2D,
                                        Pdelta=second_order,
                                        ineleastic=ineleastic)
        pass
    #
    #
    def buckling(self, name: str |None = None,
                 second_order: bool = False,
                 ineleastic: bool = False):
        """ Eigen Buckling"""
        if not name:
            name = self._mesh._name
        self._name = name
        #
        result_id = self._push_analysis(name=self._name,
                                        analysis_type='static',
                                        plane=self._plane2D,
                                        Pdelta=second_order,
                                        ineleastic=ineleastic)
        pass
    #
    #
    def time_history(self, name: str |None = None,
                     second_order: bool = False,
                     ineleastic: bool = False):
        """ """
        if not name:
            name = self._mesh._name
        self._name = name
        #
        result_id = self._push_analysis(name=self._name,
                                        analysis_type='static',
                                        plane=self._plane2D,
                                        Pdelta=second_order,
                                        ineleastic=ineleastic)
        pass
    #
    #
    def dynamic(self, end_time: float, delta_t: float,
                name: str | None = None,
                type: str | None = None,
                second_order: bool = False,
                ineleastic: bool = False
                ):
        """
        Solves the dynamic system

        end_time: simulation time [seconds]
        deltat : time increment [seconds]
        type : modal, time history
        mass : lumped, consistent
        damping :
        """
        if not name:
            name = self._mesh._name
        self._name = name
        #
        result_id = self._push_analysis(name=self._name,
                                        analysis_type='static',
                                        plane=self._plane2D,
                                        Pdelta=second_order,
                                        ineleastic=ineleastic)
        # file = open( "stfmx.f2u", "rb" )
        # jbc = pickle.load( file )
        # stf = pickle.load( file )
        # mass =  pickle.load( file )
        # file.close()
        #
        ibandm = 1
        # load = []
        npt = int(end_time // delta_t)
        #
        trnsient(stf, mass, jbc, npt,
                 # load,
                 # disp, vel, acc, fmag, olddis,
                 # wk, damp, loadin, maxnode,
                 ibandm)
        #
        print('--')
        #

    def nfreq(self, name: str |None = None,
              second_order: bool = False,
              ineleastic: bool = False
              ):
        """ """
        if not name:
            name = self._mesh._name
        self._name = name
        #
        result_id = self._push_analysis(name=self._name,
                                        analysis_type='static',
                                        plane=self._plane2D,
                                        Pdelta=second_order,
                                        ineleastic=ineleastic)
        # geometry
        elements = self._mesh.elements()
        nodes = self._mesh.nodes()
        boundaries = self._mesh.boundaries()
        # pure python solution
        # assemble_banded_Kmatrix(elements, nodes, boundaries)
        #
        # mss, ibandm = form_mass(elements, nodes, boundaries)
        #eigen(ibandm=ibandm, ivib=2)
        print('--')

    #
    #
#