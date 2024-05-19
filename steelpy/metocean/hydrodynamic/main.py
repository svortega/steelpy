# 
# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations

#
# package imports
#
from steelpy.metocean.hydrodynamic.marine_growth import MarineGrowth
from steelpy.metocean.hydrodynamic.morison import CdCmCoefficients
from steelpy.metocean.hydrodynamic.wkf import WaveKinFactor
from steelpy.metocean.hydrodynamic.cbf import CurrentBlockFactor
from steelpy.metocean.hydrodynamic.element_segment import ElementSegmentation
#
from steelpy.utils.sqlite.utils import create_connection, create_table

#
#@dataclass
class HydroProperty:
    """
    """
    __slots__ = ['conductor_shielding',
                 '_wip',
                 #'_buoyancy_area',
                 '_cdcm', '_wkf', '_cbf', 
                 #'_air_drag',
                 '_marine_growth',
                 # 'flooding', 
                 #'_non_hydro',
                 'db_file'] # 'rho_w', 
    
    def __init__(self, db_file: str) -> None:
        """
        """
        #
        self.db_file = db_file
        # TODO : class to SQL
        #self.flooding = Flooding()
        #self._hydro_diametre = HydroDiametre()
        #self._marine_growth = MarineGrowth()
        # TODO: define classes
        #self._non_hydro = {}
        #self._buoyancy_area = {}
        #self._air_drag = {}
        #self.rho_w: float = rho_w  # 1032.0  * kg / m^3
        #
        self._marine_growth = MarineGrowth(db_file=self.db_file)        
        #
        self._cdcm = CdCmCoefficients(db_file=self.db_file)
        self._wkf = WaveKinFactor(db_file=self.db_file)
        self._cbf = CurrentBlockFactor(db_file=self.db_file)
        self._wip = ElementSegmentation(db_file=self.db_file)
        #
        # create table
        conn = create_connection(self.db_file)
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
        table = "CREATE TABLE IF NOT EXISTS Property (\
                    number INTEGER PRIMARY KEY NOT NULL,\
                    condition_id INTEGER NOT NULL REFERENCES Condition(number),\
                    mg_id INTEGER REFERENCES MarineGrowth(number),\
                    cdcm_id INTEGER REFERENCES CdCm(number), \
                    element_segment_id INTEGER REFERENCES ElementSegmentation(number), \
                    flooding_id INTEGER, \
                    cshielding_id INTEGER REFERENCES ConductorShielding(number), \
                    title TEXT);"
        create_table(conn, table)    
    #
    # ---------------------------------------
    #
    @property
    def rho(self) -> float:
        """Water density [kg / m^3]"""
        data = ('metocean', )
        table = "SELECT name FROM Criteria \
                 WHERE type = ?"
        conn = create_connection(self.db_file)
        with conn:
            cur = conn.cursor()
            cur.execute(table, data)
        1 / 0
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
                self._wip.df(df)
            except AttributeError:
                pass
        #
        return self._wip
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
    #
    #
    #@property
    #def marine_growth(self):
    #    """
    #    """
    #    return self._marine_growth
    #
    #@property
    #def diameter(self):
    #    """
    #    """
    #    return self._hydro_diametre
    #
    #@property
    #def buoyancy_area(self):
    #    """
    #    """
    #    return self._buoyancy_area
    #
    #@property
    #def non_hydro(self):
    #    """
    #    """
    #    return self._non_hydro
    #
    #
    #@property
    #def air_drag(self):
    #    """
    #    """
    #    return self._air_drag
#
#
