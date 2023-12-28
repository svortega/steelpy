#
# Copyright (c) 2009 steelpy
#
from __future__ import annotations
#
# Python stdlib imports
#from typing import NamedTuple
from dataclasses import dataclass
import re
from math import tau
#import sqlite3

# package imports
from steelpy.metocean.wave.utils.main import WaveBasic
from steelpy.metocean.wave.regular.process.waveops import get_wave_data #, WaveItem
from steelpy.metocean.wave.regular.stokes.Stokes import  WaveStokes
from steelpy.metocean.wave.regular.fourier.Fourier import WaveFourier
from steelpy.metocean.wave.regular.cnoidal.Cnoidal import WaveCnoidal
#
from steelpy.utils.sqlite.utils import create_connection, create_table
from steelpy.metocean.wave.regular.process.surface import SurfaceResults
from steelpy.metocean.wave.regular.process.kinematic import KinematicResults
from steelpy.utils.dataframe.main import DBframework

from numpy.matlib import repmat
import numpy as np

#
#
#
class RegularWaves(WaveBasic):
    __slots__ = ['db_file', '_condition', 'nstep', 'niter', 'accuracy']

    def __init__(self, db_file: str,
                 nstep:int=2, niter:int=40, accuracy:float=1e-6)-> None:
        """
        db_file  : sqlite database name
        nsteps   : Number of height steps to reach H/d (2)
        niter    : Maximum number of iterations for each step (20)
        accuracy : Criterion for convergence (1e-6)
        """
        super().__init__(db_file)
        self._condition = 2
        self.nstep = nstep
        self.niter = niter
        self.accuracy = accuracy
    #
    #
    def __setitem__(self, name: int|str,
                    case_data: list[float]|dict[str, float]) -> None:
        """
        case_name : Wave name
         H : Wave height [unit length]
         T : Wave Period [second]
         d : Water depth LTH + Tide and Surge [unit length]
         Lw : Wave Length [unit length]
         Phase : Wave phase [degree]
         order :  Number of Fourier components or Stokes' order or cnoidal theory.
         crest_elevation: [unit lenght]
        """
        #
        #self._cls._input(name, 'regular')
        #
        data = get_wave_data(case_data)
        wave_type = data.pop()
        Lw = data.pop()
        #
        #if re.match(r"\b(stoke(s)?)\b", wave_type, re.IGNORECASE):
        #    self._wave[name] = 'stokes'
        #    self._stokes[name] = data
        #
        #elif re.match(r"\b(cnoidal)\b", wave_type, re.IGNORECASE):
        #    self._wave[name] = 'cnoidal'
        #    self._cnoidal[name] = data
        #else:
        #    self._wave[name] = 'fourier'
        #    self._fourier[name] = data
        #
        #
        wave_data = (name, 'regular', *data, wave_type, None, None, Lw, None)
        conn = create_connection(self.db_file)
        with conn:
            self._push_data(conn, wave_data)
    #
    #
    def __getitem__(self, name: int|str) -> tuple:
        """
        case_name : Wave name
        """
        #
        conn = create_connection(self.db_file)
        with conn:
            data = self._pull_data(conn, wave_name=name)
        #
        #
        return RegWaveItem(number=data[0], name=data[1],
                           Hw=data[3], Tw=data[4], d=data[5], 
                           theory=data[6], db_file=self.db_file)
    #
    @property
    def _label(self):
        """ """
        conn = create_connection(self.db_file)
        with conn:        
            table = "SELECT * FROM tb_Wave \
                     WHERE type = 'regular'"
            cur = conn.cursor()
            cur.execute(table)
            rows = cur.fetchall()
        #
        labels = set([item[1] for item in rows])
        return list(labels)
    #
    # ------------------
    # SQL ops
    # ------------------
    #
    def _create_table(self, conn) -> None:
        """ """
        # Wave main
        table = "CREATE TABLE IF NOT EXISTS tb_WaveRegular (\
                    number INTEGER PRIMARY KEY NOT NULL,\
                    wave_number INTEGER NOT NULL REFERENCES tb_Wave(number),\
                    height  DECIMAL NOT NULL,\
                    period DECIMAL NOT NULL,\
                    water_depth DECIMAL NOT NULL,\
                    wave_theory TEXT,\
                    wave_order DECIMAL,\
                    phase DECIMAL, \
                    wave_length DECIMAL,\
                    crest_elevation DECIMAL);"
        create_table(conn, table)
        #
        # Wave surface
        table = "CREATE TABLE IF NOT EXISTS tb_WaveSurface (\
                    number INTEGER PRIMARY KEY NOT NULL,\
                    wave_number INTEGER NOT NULL REFERENCES tb_Wave(number), \
                    type TEXT NOT NULL, \
                    x INTEGER NOT NULL, \
                    eta INTEGER NOT NULL, \
                    phase INTEGER NOT NULL);"
        create_table(conn, table)
        #
        # Wave Kinematics
        table = "CREATE TABLE IF NOT EXISTS tb_WaveKinematics (\
                    number INTEGER PRIMARY KEY NOT NULL,\
                    wave_number INTEGER NOT NULL REFERENCES tb_Wave(number), \
                    surface_number INTEGER NOT NULL REFERENCES tb_WaveSurface(number), \
                    type TEXT NOT NULL, \
                    z INTEGER NOT NULL, \
                    u INTEGER NOT NULL, \
                    v INTEGER NOT NULL, \
                    dphidt INTEGER NOT NULL, \
                    ut INTEGER NOT NULL, \
                    vt INTEGER NOT NULL, \
                    ux INTEGER NOT NULL, \
                    uz INTEGER NOT NULL, \
                    pressure INTEGER NOT NULL, \
                    Bernoulli_check INTEGER NOT NULL);"
        create_table(conn, table)
        #
    #
    #
    def _push_data(self, conn, wave_data):
        """
        Create a new project into the projects table
        """
        project = (wave_data[0], wave_data[1])
        table = 'INSERT INTO tb_Wave(name, type) \
                 VALUES(?,?)'
        # push
        cur = conn.cursor()
        cur.execute(table, project)
        wave_number = cur.lastrowid
        #
        #
        project = (wave_number, *wave_data[2:])
        cur = conn.cursor()
        table = 'INSERT INTO tb_WaveRegular(wave_number, \
                                            height, period, water_depth,  \
                                            wave_theory, wave_order, phase, \
                                            wave_length, crest_elevation) \
                                            VALUES(?,?,?,?,?,?,?,?,?)'
        # push
        cur = conn.cursor()
        cur.execute(table, project)
    #
    def _pull_data(self, conn, wave_name: str|int,
                  item: str = "*"):
        """ """
        project = (wave_name, )
        table = "SELECT * FROM tb_Wave \
                 WHERE name = ? AND type = 'regular'"
        cur = conn.cursor()
        cur.execute(table, project)
        wave = cur.fetchone()
        #
        project = (wave[0],)
        sql = 'SELECT {:} FROM tb_WaveRegular WHERE wave_number = ?'.format(item)
        cur = conn.cursor()
        cur.execute(sql, project)
        data = cur.fetchone()
        return [*wave[:3], *data[2:]]
    #
    # ------------------------------
    #
    def df(self, df):
        """ """
        grpname = df.groupby(['name'])
        for key, item in grpname:
            values = item[['Hw', 'Tw', 'd', 'theory']].values
            for idx, value in enumerate(values):
                # TODO : define input better
                #name = f'{key[0]}_{idx + 1}'
                name = key[0]
                self.__setitem__(name=name,
                                 case_data=value.tolist())
        #print('--')
    #
    # ------------------------------
    #
    def solve(self, surface_points:int = 36,
              depth_points:int = 100):
        """ """
        #print('-->')
        for name in self._label:
            wave = self.__getitem__(name)
            wave(surface_points=surface_points,
                 depth_points=depth_points)
        #print('-->')
        #1 / 0
    #
    # ------------------------------
    #
    @property
    def Stokes(self):
        """ """
        return WaveStokes
    #
    @property
    def Fourier(self):
        """ """
        return WaveFourier
    #
    @property
    def Cnoidal(self):
        """ """
    #    print('stokes')
    #    1 / 0
        return WaveCnoidal
    #
    #
#
#
@dataclass
class RegWaveBasic:
    #1 / 0
    number:int
    name:str|int
    Hw:float
    Tw:float
    d:float
    theory:str 
    #title:str
    #Lw:float|None=None   
    #direction:float = 0
    #kinematic_factor:float = 1.0
    #crest_elevation:float|bool = None
    #title:str
    #
    #@property
    #def Stokes(self):
    #    """ """
    #    return WaveStokes(H=self.Hw, T=self.Tw, d=self.d,
    #                      title=self.name)
    #
    #@property
    #def Cnoidal(self):
    #    """ """
    #    return CnoidalModule
    #
    #@property
    #def Fourier(self):
    #    """ """
    #    return FourierModule
    #
    #@property
    #def Lw(self):
    #    """ """
    #    1 / 0
    #
    #
    #def kinematics(self):
    #    """vave kinemactis"""
    #    print('here')    
#
#
@dataclass
class RegWaveItem: # (WaveItem)
    
    __slots__ = ['number', 'name', 'db_file', '_wave', 'finite_depth']
   
    
    def __init__(self, number:int, name:str|int, 
                 Hw:float, Tw:float, d:float, db_file:str, 
                 Lw:float|None=None, theory: str = 'fourier', 
                 infinite_depth:bool=False,
                 current:float = 0.0, c_type:int = 1,
                 order:int=5, nstep:int=2,
                 niter:int=40, accuracy:float=1e-6) -> None:
        """
        """
        #super().__init__(H=Hw, Tw=Tw, d=d, title=name,
        #                 order=order, nstep=nstep, niter=niter,
        #                 accuracy=accuracy, 
        #                 current=current, c_type=c_type,
        #                 infinite_depth=infinite_depth)
        #
        self.number = number
        self.name = name
        self.db_file = db_file
        #
        if theory.lower() == 'stokes':
            self._wave = WaveStokes(H=Hw, T=Tw, d=d, title=name,
                                    order=order, nstep=nstep, niter=niter,
                                    accuracy=accuracy, 
                                    current=current, c_type=c_type,
                                    infinite_depth=infinite_depth)
        #elif theory.lower() == 'cnoidal':
        #    WaveCnoidal.__init__(self)
        else:
            self._wave = WaveFourier(H=Hw, T=Tw, d=d, title=name,
                                     order=order, nstep=nstep, niter=niter,
                                     accuracy=accuracy, 
                                     current=current, c_type=c_type,
                                     infinite_depth=infinite_depth)
        #
        self.finite_depth = True
        if infinite_depth:
            self.finite_depth = False
        #
    #
    def __call__(self, surface_points:int,
                 depth_points:int):
        """ Solver """
        self._wave()
        #wave_length = (tau / self._z[ 1 ]) * self.d
        #
        conn = create_connection(self.db_file)
        #
        # Update wave lenght
        with conn:
            items = (self.number, )
            table = f"UPDATE tb_WaveRegular \
                     SET wave_length = {self._wave._Lw} \
                     WHERE wave_number = ?"
            cur = conn.cursor()
            cur.execute(table, items)
        #
        # Calculate Surface
        #
        surface = self._wave.get_surface(surface_points)
        #surface['wave_name'] = self.name
        surface['wave_number'] = self.number        
        #surface =  surface[['wave_number', 'type', 'x', 'eta', 'phase']]
        #
        conn = create_connection(self.db_file)
        with conn:
            # push surface data
            surface.to_sql('tb_WaveSurface', conn,
                         index_label=['wave_number', 'type', 'x', 'eta', 'phase'], 
                         if_exists='append', index=False)        
            #
            # get surface data from sql 
            surface = self._pull_surface(conn)
        #
        #
        kindf = self._wave.get_kinematics(depth_points=depth_points)
        #
        #kpoints = depth_points
        #kindf = get_kinematic(self.order, self._z, self._B, self._Tanh,
        #                      self.d, surface, kpoints,
        #                      self.finite_depth)
        #
        kindf['wave_number'] = self.number
        #kindf['type'] = 'order_1'
        #
        depth_steps = depth_steps = np.arange(depth_points + 1) / depth_points
        kindf['surface_number'] = repmat(surface['number'].to_numpy(),
                                         depth_steps.size, 1).flatten('F')
        #
        kindf.drop(columns=['phase', 'x'], inplace=True, axis=1)
        #
        conn = create_connection(self.db_file)
        with conn:
            kindf.to_sql('tb_WaveKinematics', conn,
                         index_label=['wave_number', 'surface_number','type', 
                                      'z', 'u', 'v', 'dphidt',
                                      'ut', 'vt', 'ux', 'uz',
                                      'pressure', 'Benoulli_check'], 
                         if_exists='append', index=False)        
        #
        #print('---> here')
    #
    @property
    def Lw(self):
        """ Wave Length"""
        conn = create_connection(self.db_file)
        with conn:        
            data = self._pull_wave(conn)
            wave_length = data[9]
        #
        if not wave_length:
            #self.solver()
            #wave_length = (tau / self._z[ 1 ]) * self.d
            #conn = create_connection(self.db_file)
            # Update wave lenght
            #with conn:
            #    items = (self.name, )
            #    table = f"UPDATE tb_Wave \
            #             SET wave_length = {wave_length} \
            #             WHERE name = ?"
            #    cur = conn.cursor()
            #    cur.execute(table, items)            
            self.__call__()
            data = self._pull_wave(conn)
            wave_length = data[9]            
        #
        return wave_length
    #
    @property
    def wave(self):
        """ """
        conn = create_connection(self.db_file)
        with conn:        
            data = self._pull_wave(conn)
        #
        return RegWaveBasic(number=data[0], 
                            name=data[1],
                            Hw=data[3], 
                            Tw=data[4], 
                            d=data[5],
                            theory=data[6])
    #
    def _pull_wave(self, conn, item: str = "*"):
        """ """
        items = (self.name,)
        table = "SELECT * FROM tb_Wave \
                 WHERE name = ? AND type = 'regular'"
        cur = conn.cursor()
        cur.execute(table, items)
        wave = data = cur.fetchone()
        #return data
        project = (wave[0],)
        sql = 'SELECT {:} FROM tb_WaveRegular WHERE wave_number = ?'.format(item)
        cur = conn.cursor()
        cur.execute(sql, project)
        data = cur.fetchone()
        return [*wave[:3], *data[2:]]
    #
    #
    def surface(self, surface_points:int = 36,
                step_size:float|None=None):
        """Free surface (wave) profile
        surface_points : Number of points on free surface (the program clusters them near crest)
        step_size: deg
        """
        conn = create_connection(self.db_file)
            #
            #df = DBframework()
            #header = ['wave_name', 'number', 'wave_number', 'type', 'x', 'eta', 'phase']
            #surface = df.DataFrame(data=data, columns=header)
            #surface = surface[['wave_name', 'type', 'x', 'eta', 'phase']]            
        #
        try:
            with conn:
                surface = self._pull_surface(conn)            
            surface.columns
        except AttributeError:
            #surface = self.get_surface(surface_points)
            #surface =  surface[['wave_number', 'type', 'x', 'eta', 'phase']]
            #
            #conn = create_connection(self.db_file)
            #with conn:
            #    surface.to_sql('tb_WaveSurface', conn,
            #                 index_label=['wave_number', 'type', 'x', 'eta', 'phase'], 
            #                 if_exists='append', index=False)
            #    surface['wave_name'] = self.name
            #
            self.__call__(surface_points=surface_points)
            with conn:
                surface = self._pull_surface(conn)            
        #
        #if sdata:
        #    print('df here')
        #    1 / 0
        #    surface 
        #else:
        #n = self.order
        #try:
        #    kd = self._z[1]
        #except :
        #    self.__call__()
        #    kd = self._z[1]
        #
        #surface = get_surface(n, kd, self._Y, self.d, 
        #                      surface_points, self.finite_depth)
        #
        #surface['wave_number'] = self.number
        #surface['type'] = 'order_1'
        #
        #conn = create_connection(self.db_file)
        #with conn:
        #    surface.to_sql('tb_WaveSurface', conn,
        #                 index_label=['wave_number', 'type', 'x', 'eta', 'phase'], 
        #                 if_exists='append', index=False)
        #
        #surface['wave_name'] = self.name
        #return surface[['wave_name', 'type', 'x', 'eta', 'phase']]
        #
        #return SurfaceResults(surface, self.H, self.Tw, self.d, self.finite_depth)
        #
        #
        #with conn:        
        #    data = self._pull_wave(conn)
        #    data
        #    1 / 0
        wave = self.wave
        #
        return SurfaceResults(surface=surface,
                              Hw=wave.Hw, Tw=wave.Tw, d=wave.d,
                              finite_depth=self.finite_depth)
    #
    def _pull_surface(self, conn) -> list | None:
        """get wave surface"""
        #
        #wave = self._get_wave(self, conn)
        #
        items = (self.number,)
        table = 'SELECT tb_Wave.name, tb_WaveSurface.* \
                 FROM tb_Wave, tb_WaveSurface \
                 WHERE tb_WaveSurface.wave_number = ? \
                 AND tb_Wave.number = tb_WaveSurface.wave_number'
        #try:
        cur = conn.cursor()
        cur.execute(table, items)
        data = cur.fetchall()
        #
        df = DBframework()
        header = ['wave_name', 'number', 'wave_number', 'type', 'x', 'eta', 'phase']
        surface = df.DataFrame(data=data, columns=header)
        return surface[['number', 'wave_name', 'wave_number', 'type', 'x', 'eta', 'phase']]
        #except sqlite3.OperationalError:
        #    return None
    #
    def kinematics(self, depth_points:int = 100,
                   surface_points:int = 36):
        """ """
        #
        conn = create_connection(self.db_file)
        #
        #self.__call__()
        try:
            with conn:
                kindf = self._pull_kinematics(conn)            
            kindf.columns
        except AttributeError:
            #surface = self.surface(surface_points)
            #
            #n = self.order
            #kd = self._z[1]
            #
            #crest = surface.eta
            #
            #depth = np.arange(depth_points + 1) / depth_points
            #
            #zdepth =  crestmax * depth  # --> fix sign
            #zdepth = surfacePoints(d=self.d, points=depth_points, eta=crest)
            #
            #get_kinematicX(n, self._z, self._Y, self._B, self._Tanh,
            #               surface_points, depth_points, self.finite_depth)
            #
            #1 / 0
            #kpoints = depth_points
            #kindf = get_kinematic(self.order, self._z, self._B, self._Tanh,
            #                      self.d, surface, kpoints,
            #                      self.finite_depth)
            #
            self.__call__(surface_points=surface_points,
                          depth_points=depth_points)
            with conn:
                kindf = self._pull_kinematics(conn)            
        #
        with conn:
            surface = self.surface(conn)
        #
        #return KinematicResults(surface, kindf, depth_points+1)
        return KinematicResults(surface=surface,
                                kindata=kindf,
                                depth_points=depth_points)
        #
        #return dataframe
    #
    def _pull_kinematics(self, conn):
        """read kin from sql"""
        items = (self.number,)
        table = 'SELECT tb_Wave.name, tb_WaveSurface.x, tb_WaveSurface.phase, \
                 tb_WaveKinematics.* \
                 FROM tb_Wave, tb_WaveSurface, tb_WaveKinematics \
                 WHERE tb_WaveSurface.wave_number = ? \
                 AND tb_Wave.number = tb_WaveKinematics.wave_number \
                 AND tb_WaveSurface.number = tb_WaveKinematics.surface_number'
        #try:
        cur = conn.cursor()
        cur.execute(table, items)
        data = cur.fetchall()
        #
        df = DBframework()
        header = ['wave_name', 'x', 'phase', 'index', 
                  'wave_number', 'surface_number','type', 
                  'z', 'u', 'v', 'dphidt',
                  'ut', 'vt', 'ux', 'uz',
                  'pressure', 'Benoulli_check']
        surface = df.DataFrame(data=data, columns=header)
        return surface