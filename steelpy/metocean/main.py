# 
# Copyright (c) 2019 steelpy
# 
# Python stdlib imports
from __future__ import annotations
#from collections.abc import Mapping
from datetime import datetime as dt
#import os
import re


# package imports
from steelpy.utils.sqlite.main import ClassMainSQL
#from steelpy.utils.io_module.inout import check_input_file
#from steelpy.metocean.irregular.main import WaveIrregular
from steelpy.metocean.wave.regular.process.bsotm import BSOTM
from steelpy.metocean.wave.main import WaveMain
#from steelpy.metocean.irregular.spectrum import Sprectrum
from steelpy.utils.units.main import Units
#from steelpy.metocean.interface.wave import RegularWaves, IregularWaves
from steelpy.metocean.wind.main import Wind
from steelpy.metocean.current.main import Current
#from steelpy.metocean.hydrodynamic.marine_growth import MarineGrowth
#from steelpy.metocean.hydrodynamic.morison import CdCmCoefficients
#from steelpy.metocean.hydrodynamic.wkf import WaveKinFactor
#from steelpy.metocean.hydrodynamic.cbf import CurrentBlockFactor
from steelpy.metocean.hydrodynamic.main import Hydrodynamic
from steelpy.metocean.process.conditions import MetoceanCondition
from steelpy.utils.sqlite.utils import create_table
#
#
class Metocean(ClassMainSQL):
    """
    FE Metocean Class
    
    Metocean
        |_ name
        |_ number
        |_ data
        |_ type
        |_ wave
        |_ current
        |_ wind
        |_ cdcm
        |_ non_hydro
        |_ elevation
        |_ hydro_diameter
        |_ buoyancy
        |_ flooded
        |_ seastate
    
    **Parameters**:  
      :number:  integer internal number 
      :name:  string node external name
      
      : seastate : metocean combination
    """
    #
    __slots__ = ['_name', '_wave', '_current', '_wind',
                 '_cdcm','_combination', '_rho_w',
                 '_properties',
                 #'_marine_growth','_wkf', '_cbf', 
                 '_build', 'db_file']
    #
    def __init__(self, name:str|None = None,
                 sql_file:str|None = None):
        """
        """
        self._rho_w: float = 1032.0  # kg / m^
        #self._build = True        
        #
        super().__init__(name=name, sql_file=sql_file)
        #
        #
        #if sql_file:
        #    sql_file = check_input_file(file=sql_file,
        #                                file_extension="db")
        #    self.db_file = sql_file
        #    self._build = False
        #else:
        #    self.db_file = self._get_file(name=name)
        #    #self.data_type = mesh_type
        #    self._name = name
        #    conn = create_connection(self.db_file)
        #    with conn:
        #        self._create_table(conn)
        #
        #
        self._wave = WaveMain(rho_w=self._rho_w,
                              db_file=self.db_file)
        
        self._wind = Wind(db_file=self.db_file)
        
        self._current = Current(db_file=self.db_file)
        #
        #self._marine_growth = MarineGrowth(rho_w=self._rho_w,
        #                                   db_file=self.db_file)
        #
        #self._cdcm = CdCmCoefficients(db_file=self.db_file)
        #self._wkf = WaveKinFactor(db_file=self.db_file)
        #self._cbf = CurrentBlockFactor(db_file=self.db_file)
        #
        self._properties = Hydrodynamic(rho_w=self._rho_w,
                                        db_file=self.db_file)
        #
        self._combination = MetoceanCondition(db_file=self.db_file,
                                              properties=self._properties)
    #
    #def _get_file(self, name: str):
    #    """ """
    #    filename = name + ".db"
    #    path = os.path.abspath(filename)
    #    try: # remove file if exist
    #        os.remove(path)
    #    except FileNotFoundError:
    #        pass
    #    #
    #    return path
    #
    #
    # ------------------
    # SQL ops
    # ------------------
    #
    def _create_table(self, conn) -> None:
        """ """
        # conn = create_connection(self.db_file)
        table = "CREATE TABLE IF NOT EXISTS tb_Main (\
                    number INTEGER PRIMARY KEY NOT NULL,\
                    name NOT NULL,\
                    type TEXT NOT NULL,\
                    units TEXT NOT NULL,\
                    rho_water DECIMAL NOT NULL,\
                    date TEXT);"
        #
        create_table(conn, table)
        #
        #
        table = 'INSERT INTO tb_Main(name, type, units,\
                                    rho_water, date)\
                                    VALUES(?,?,?,?,?)'
        #
        time=dt.now().strftime('%Y-%m-%d')
        data = (self._name, 'metocean', 'si', self._rho_w, time)
        # push
        cur = conn.cursor()
        cur.execute(table, data)
    #
    # ------------------------------------------
    #
    @property
    def rho_w(self):
        """ """
        units = Units()
        return self._rho_w * units.kg / units.m**3

    @rho_w.setter
    def rho_w(self, value:Units):
        """ """
        self._rho_w = value.convert('kilogram/metre^3').value
    #
    # ------------------------------------------
    #
    def wave(self, values:None|list=None,
             df=None):
        """ """
        if values:
            print('-->')
            1/0
        else:
            try:
                columns = list(df.columns)
                grpwave = df.groupby('type')
                for wtype, item in grpwave:
                    if re.match(r"\b(regular(_wave)?)\b", wtype, re.IGNORECASE):
                        self._wave._regular.df(df=item)
                    elif re.match(r"\b(iregular(_wave)?)\b", wtype, re.IGNORECASE):
                        # self._wave._iregular.df(df=item)
                        raise NotImplementedError
                    else:
                        raise ValueError("wave type invalid")
            except AttributeError:
                pass
        #
        return self._wave
    #
    #@property
    def wind(self, values:None|list=None,
             df=None):
        """
        """
        if values:
            print('-->')
            1/0
        else:
            try:
                df.columns            
                self._wind.df(df)
            except AttributeError:
                pass
        #
        return self._wind
    #
    #@property
    def current(self, values:None|list=None,
                df=None):
        """
        """
        if values:
            print('-->')
            1/0
        else:
            try:
                df.columns            
                self._current.df(df)
            except AttributeError:
                pass
        return self._current
    #
    #
    # ------------------------------------------
    #@property
    def properties(self):
        """
        """
        return self._properties
    #
    # ------------------------------------------
    #@property
    def condition(self, values:None|list=None,
                 df=None):
        """
        """
        if values:
            print('-->')
            1/0
        else:
            try:
                df.columns            
                self._combination.df(df)
            except AttributeError:
                pass        
        return self._combination
    #
    #
    # ------------------------------------------
    #
    #
    def solve(self, surface_points:int = 36,
              depth_points:int = 100):
        """ """
        # wave
        self._wave.solve(surface_points=surface_points,
                         depth_points=depth_points)
        #
        # current
        #
        # wind
        #print('-->')
        #1 / 0
    #
    #
    # ------------------------------------------
    #
    def get_load(self, mesh, kinematic,
                 condition:int|None = None,
                 rho:float = 1025):
        """
        condition :
            1 - Linear wave (dafault)
            2 - Non-linear wave
        """
        #if not condition:
            
        wforce = BSOTM(kinematic, condition, rho)
        wforce.wave_force(mesh=mesh)
        print('-->')
    #
    #
    def pile_response(self, D:float|Units,
                      L:float|Units,
                      kinematic,
                      condition: int = 1,
                      rho:float = 1025):
        """
        D : Pile diametre
        L : Pile length
        kinematics : kinematic dataframe
        condition :
            1 - Linear wave (dafault)
            2 - Non-linear wave
        """
        Dp = D.convert('metre').value
        Lp = L.convert('metre').value
        #
        #if kinematic._type == 'regular':
        #bs, ovtm = bsotm_reg(kinematic, D_pile, condition)
        #else:
        #    #bs, ovtm = bsvtm(kinematic, D_pile, condition)
        #    raise NotImplemented
        #
        bsotm = BSOTM(kinematic, condition, rho)
        bs, otm = bsotm.solveBSOTM(D=Dp, L=Lp)
        #
        return bs, otm 
        #return surface
#
#

