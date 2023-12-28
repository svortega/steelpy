# 
# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations

#
# package imports
#from .marine_growth import MarineGrowth
#from .morison import CdCmCoefficients
from .hydro_diametre import HydroDiametre
from .flooding import Flooding
#
from steelpy.metocean.hydrodynamic.marine_growth import MarineGrowth
from steelpy.metocean.hydrodynamic.morison import CdCmCoefficients
from steelpy.metocean.hydrodynamic.wkf import WaveKinFactor
from steelpy.metocean.hydrodynamic.cbf import CurrentBlockFactor
#
from steelpy.utils.sqlite.utils import create_connection, create_table

#
#@dataclass
class Hydrodynamic:
    """
    """
    __slots__ = ['flooding', 
                 'conductor_shielding',
                 'element_refining',
                 '_buoyancy_area',
                 '_cdcm', '_wkf', '_cbf', 
                 '_air_drag',
                 '_marine_growth',
                 '_hydro_diametre',
                 '_non_hydro',
                 'rho_w', '_db_file']
    
    def __init__(self, rho_w:float,
                 db_file: str) -> None:
        """
        """
        #
        self._db_file = db_file
        #
        self.flooding = Flooding()
        self._hydro_diametre = HydroDiametre()
        #self._marine_growth = MarineGrowth()
        # TODO: 
        self._non_hydro = {}
        self._buoyancy_area = {}
        self._air_drag = {}
        self.rho_w: float = rho_w  # 1032.0  * kg / m^3
        #
        self._marine_growth = MarineGrowth(rho_w=self.rho_w,
                                           db_file=self._db_file)        
        #
        self._cdcm = CdCmCoefficients(db_file=self._db_file)
        self._wkf = WaveKinFactor(db_file=self._db_file)
        self._cbf = CurrentBlockFactor(db_file=self._db_file)
        #
        # create table
        conn = create_connection(self._db_file)
        with conn:
            self._create_table(conn)        
    #
    # ------------------
    # SQL ops
    # ------------------
    #
    def _create_table(self, conn) -> None:
        """ """
        #
        table = "CREATE TABLE IF NOT EXISTS tb_Properties (\
                    number INTEGER PRIMARY KEY NOT NULL,\
                    condition_number INTEGER NOT NULL REFERENCES tb_Condition(number),\
                    mg_number INTEGER REFERENCES tb_MarineGrowth(number),\
                    cdcm_number INTEGER REFERENCES tb_CdCm(number), \
                    flooding_number INTEGER, \
                    cshielding_number INTEGER, \
                    element_refinament INTEGER, \
                    title TEXT);"
        create_table(conn, table)    
    #
    # ---------------------------------------
    #
    def MG(self, values:None|list=None,
           df=None):
        """Marine Growth"""
        if values:
            print('-->')
            1/0
        else:
            try:
                df.columns            
                self._marine_growth.df(df)
            except AttributeError:
                pass           
        return self._marine_growth
    #
    #
    def WKF(self, values:None|list=None,
            df=None):
        """ Wave kinematic factor"""
        if values:
            print('-->')
            1/0
        else:
            try:
                df.columns            
                self._wkf.df(df)
            except AttributeError:
                pass
        #
        return self._wkf
    #
    #
    def CBF(self, values:None|list=None,
            df=None):
        """ Current Blockage factor"""
        if values:
            print('-->')
            1/0
        else:
            try:
                df.columns            
                self._cbf.df(df)
            except AttributeError:
                pass
        #
        return self._cbf    
    #
    #
    def WIP(self, values:None|list=None,
            df=None):
        """ Wave integration points"""
        if values:
            print('-->')
            1/0
        else:
            try:
                df.columns            
                #self._cbf.df(df)
            except AttributeError:
                pass
        #
        #return self._cbf
        1 / 0
    #
    # ---------------------------------------
    #
    #@property
    def CdCm(self):
        """
        """
        return self._cdcm
    #
    #@property
    #def Morison(self):
    #    """
    #    """
    #    return self._cdcm
    #
    @property
    def flooded(self):
        """
        """
        return self.flooding
    #
    #@property
    #def marine_growth(self):
    #    """
    #    """
    #    return self._marine_growth
    #
    @property
    def diameter(self):
        """
        """
        return self._hydro_diametre
    #
    @property
    def buoyancy_area(self):
        """
        """
        return self._buoyancy_area
    #
    @property
    def non_hydro(self):
        """
        """
        return self._non_hydro
    #
    #
    @property
    def air_drag(self):
        """
        """
        return self._air_drag
#
#
