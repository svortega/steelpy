#
# Copyright (c) 2009 steelpy
#
from __future__ import annotations
#
# Python stdlib imports
#from typing import NamedTuple
from dataclasses import dataclass
#import re
#from math import tau

# package imports
from steelpy.metocean.wave.utils.main import WaveBasic
from steelpy.metocean.wave.regular.process.waveops import get_wave_data
from steelpy.metocean.wave.regular.stokes.Stokes import WaveStokes
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
class RegularWave(WaveBasic):
    __slots__ = ['db_file', '_criteria', '_wave',
                 '_hpoints', '_vpoints',
                 '_condition', 'nstep', 'niter', 'accuracy']

    def __init__(self, criteria: str | int,
                 db_file: str,
                 nstep: int = 2, niter: int = 40,
                 accuracy: float = 1e-6,
                 depth_points: int = 100,
                 surface_points: int = 36
                 ) -> None:
        """
        db_file  : sqlite database name
        nsteps   : Number of height steps to reach H/d (2)
        niter    : Maximum number of iterations for each step (20)
        accuracy : Criterion for convergence (1e-6)
        """
        super().__init__(db_file)
        self._criteria = criteria
        self._condition = 2
        self.nstep = nstep
        self.niter = niter
        self.accuracy = accuracy
        #
        self._hpoints = surface_points
        self._vpoints = depth_points
        #self._wave = {}

    #
    #
    def __setitem__(self, name: int | str,
                    case_data: list[float] | dict[str, float]) -> None:
        """
        case_name : Wave name
         H : Wave height [unit length]
         T : Wave Period [second]
         d : Water depth LTH + Tide and Surge [unit length]
         Lw : Wave Length [unit length]
         order :  Number of Fourier components or Stokes' order or cnoidal theory.
         crest_elevation: [unit length]
        """
        # get user data
        wdata = get_wave_data(case_data)
        Lw = wdata.pop()
        WCe = wdata.pop()
        title = wdata.pop()
        #
        conn = create_connection(self.db_file)
        with conn:
            # [name, wave_type, title, criteria_id]
            wave_data = (name, 'regular', title)
            wave_id = self._push_wave(conn, wave_data)
            #  [wave_id, 
            #   height, period, water_depth, wave_theory, 
            #   wave_order, wave_length, crest_elevation]
            wave_reg = (wave_id, *wdata, None, Lw, WCe)
            self._push_data(conn, wave_reg)
        #
        # TODO : solve wave
        self.run_wave(wave_name=name)
        #print('--')

    #
    def __getitem__(self, name: int | str) -> tuple:
        """
        case_name : Wave name
        """
        conn = create_connection(self.db_file)
        with conn:
            data = pull_wave(conn, wave_name=name)
        #
        # TODO : WCe
        # [number, name, w_type, Hw, Tw, d, w_theory, w_order, Lw, WCe]
        return RegWaveItem(number=data[0], name=data[1],
                           db_file=self.db_file,
                           infinite_depth=False) # TODO: fix here
                           #surface_points=self._hpoints,
                           #depth_points=self._vpoints)
        #Hw=data[3], Tw=data[4], d=data[5],
        #theory=data[6], db_file=self.db_file,
        #Lw=data[8], order=data[7])

    #
    # ------------------------------
    #
    @property
    def _label(self):
        """ """
        conn = create_connection(self.db_file)
        with conn:
            table = "SELECT * FROM Wave \
                     WHERE type = 'regular'"
            cur = conn.cursor()
            cur.execute(table)
            rows = cur.fetchall()
        #
        labels = set([item[1] for item in rows])
        return list(labels)

    #
    def run_wave(self, wave_name: str | int,
                 current: float = 0.0, c_type: int = 1,
                 infinite_depth: bool = False) -> None:
        """ """
        conn = create_connection(self.db_file)
        with conn:
            wdata = self._pull_data(conn, wave_name=wave_name)
        #
        #(number, name, wave_type,
        # Hw, Tw, d, theory, order, Lw, WCe, WCtrough, Wskewness) = data
        #
        # default 5
        #if not wdata.order: order = 5
        #
        wave = wdata.wave(current=current, c_type=c_type,
                          infinite_depth=infinite_depth)
        #return wave
        wave.solve()
        #print('--')
        #
        #self._wave()
        #wave_length = (tau / self._z[ 1 ]) * self.d
        #
        # Calculate Surface
        surface = wave.get_surface(surface_points=self._hpoints)
        surface['wave_id'] = wdata.number
        #
        # Update wave length
        conn = create_connection(self.db_file)
        with conn:
            items = (wdata.number,)
            table = f"UPDATE WaveRegular \
                     SET \
                     wave_length = {wave.Lw}, \
                     wave_order = {wave.order}, \
                     wave_crest = {surface.eta.max()}, \
                     wave_trough = {surface.eta.min()}, \
                     wave_skewness = {surface.eta.max() / wave.Hw} \
                     WHERE wave_id = ?"
            cur = conn.cursor()
            cur.execute(table, items)
        #
        # update wave surface
        conn = create_connection(self.db_file)
        with conn:
            # push surface data
            surface.to_sql('WaveSurface', conn,
                           index_label=['wave_id', 'type', 'length',
                                        'eta', 'phase', 'time'],
                           if_exists='append', index=False)
            #
            # get surface data from sql
            surface = pull_surface(conn,number=wdata.number)
        #
        #
        kindf = wave.get_kinematics(depth_points=self._vpoints)
        kindf['wave_id'] = wdata.number
        #kindf['type'] = 'order_1'
        #
        depth_steps = np.arange(self._vpoints + 1) / self._vpoints
        kindf['surface_id'] = repmat(surface['number'].to_numpy(),
                                     depth_steps.size, 1).flatten('F')
        #
        kindf.drop(columns=['phase', 'x', 'time'], inplace=True, axis=1)
        kindf.rename(columns={'z': 'elevation'}, inplace=True)
        #
        conn = create_connection(self.db_file)
        with conn:
            kindf.to_sql('WaveKinematic', conn,
                         index_label=['wave_id', 'surface_id', 'type',
                                      'elevation', 'u', 'v', 'dphidt',
                                      'ut', 'vt', 'ux', 'uz',
                                      'pressure', 'Benoulli_check'],
                         if_exists='append', index=False)
        #
        print('---> here')
    #
    #
    # ------------------------------
    # SQL ops
    #
    def _create_table(self, conn) -> None:
        """ """
        # Wave main
        table = "CREATE TABLE IF NOT EXISTS WaveRegular (\
                    number INTEGER PRIMARY KEY NOT NULL,\
                    wave_id INTEGER NOT NULL REFERENCES Wave(number),\
                    height  DECIMAL NOT NULL,\
                    period DECIMAL NOT NULL,\
                    water_depth DECIMAL NOT NULL,\
                    wave_theory TEXT,\
                    wave_order DECIMAL,\
                    wave_length DECIMAL,\
                    wave_crest DECIMAL,\
                    wave_trough DECIMAL,\
                    wave_skewness DECIMAL);"
        create_table(conn, table)
        #
        # Wave surface
        table = "CREATE TABLE IF NOT EXISTS WaveSurface (\
                    number INTEGER PRIMARY KEY NOT NULL,\
                    wave_id INTEGER NOT NULL REFERENCES Wave(number), \
                    type TEXT NOT NULL, \
                    length DECIMAL NOT NULL, \
                    eta INTEGER NOT NULL, \
                    phase INTEGER, \
                    time INTEGER NOT NULL);"
        create_table(conn, table)
        #
        # Wave Kinematics
        table = "CREATE TABLE IF NOT EXISTS WaveKinematic (\
                    number INTEGER PRIMARY KEY NOT NULL,\
                    wave_id INTEGER NOT NULL REFERENCES Wave(number), \
                    surface_id INTEGER NOT NULL REFERENCES WaveSurface(number), \
                    type TEXT NOT NULL, \
                    elevation INTEGER NOT NULL, \
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
    def _push_data(self, conn, wave_data: tuple) -> None:
        """
        Create a new project into the projects table
        """
        table = 'INSERT INTO WaveRegular(wave_id, \
                                        height, period, water_depth,  \
                                        wave_theory, wave_order, \
                                        wave_length, wave_crest) \
                                        VALUES(?,?,?,?,?,?,?,?)'
        # push
        cur = conn.cursor()
        cur.execute(table, wave_data)

    #
    def _pull_data(self, conn, wave_name: str | int,
                   item: str = "*") -> dataclass:
        """
        wave_name:
        item: str
        
        Return
        [number, name, w_type, Hw, Tw, d, w_theory, w_order, Lw, WCe]
        """
        wave = self._pull_wave(conn, wave_name, wave_type='regular')
        #
        project = (wave[0],)
        sql = f'SELECT {item} FROM WaveRegular \
               WHERE wave_id = ?'
        cur = conn.cursor()
        cur.execute(sql, project)
        data = cur.fetchone()
        return RegWaveBasic(*wave[:3], *data[2:])
        #return [*wave[:3], *data[2:]]

    #
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
    #def solve(self, surface_points: int = 36,
    #          depth_points: int = 100):
    #    """ """
    #    #print('-->')
    #    for name in self._label:
    #        wave = self.__getitem__(name)
    #        wave(surface_points=surface_points,
    #             depth_points=depth_points)
    #    #print('-->')
    #    #1 / 0
    #
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
        return WaveCnoidal

    #
    #
    def grid(self, depth_points: int = 100,
             surface_points: int = 36):
        """
        To define the horizontal extremity of the volume
        to generate a grid for pre-calculation of fluid
        kinematics used in a seastate simulation analysis
        """
        self._hpoints = surface_points
        self._vpoints = depth_points


#
#
@dataclass
class RegWaveBasic:
    """ """
    number: int
    name: str | int
    wave_type: str
    Hw: float
    Tw: float
    d: float
    theory: str
    order: float|None
    Lw: float|None
    wave_crest: float|None
    wave_trough: float|None
    wave_skewness: float|None
    #
    nstep: int = 2
    niter: int = 40
    accuracy: float = 1e-6
    #
    #
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
    #@property
    def wave(self, current: float = 0.0, c_type: int = 1,
             infinite_depth: bool = False):
        """return wave"""
        if not self.order: order = 5

        if self.theory.lower() == 'stokes':
            wave_theory = WaveStokes
        elif self.theory.lower() == 'cnoidal':
            wave_theory = WaveCnoidal
        else:
            wave_theory = WaveFourier
        #
        wave = wave_theory(Hw=self.Hw, Tw=self.Tw,
                           Lw=self.Lw, d=self.d,
                           order=order, title=self.name,
                           nstep=self.nstep, niter=self.niter,
                           accuracy=self.accuracy,
                           current=current, c_type=c_type,
                           infinite_depth=infinite_depth)
        #wave.solve()
        return wave
    #

#
#
@dataclass
class RegWaveItem:  # (WaveItem)

    __slots__ = ['number', 'name', 'db_file', '_wave',
                 'finite_depth', '_hpoints', '_vpoints',]

    def __init__(self, number: int, name: str | int,
                 db_file: str, infinite_depth: bool) -> None:
                 #depth_points: int,
                 #surface_points: int,                 
                 #infinite_depth: bool) -> None:
        #Hw: float, Tw: float, d: float,
        #Lw: float | None = None, theory: str = 'fourier',
        #infinite_depth: bool = False,
        #current: float = 0.0, c_type: int = 1,
        #order: int = 5, nstep: int = 2,
        #niter: int = 40, accuracy: float = 1e-6) -> None:
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
        #self._hpoints = surface_points
        #self._vpoints = depth_points        
        #self.finite_depth = infinite_depth
        #
        ## default 5
        #if not order:
        #    order = 5
        #    #
        #if theory.lower() == 'stokes':
        #    self._wave = WaveStokes(Hw=Hw, Tw=Tw, Lw=Lw, d=d, title=name,
        #                            order=order, nstep=nstep, niter=niter,
        #                            accuracy=accuracy,
        #                            current=current, c_type=c_type,
        #                            infinite_depth=infinite_depth)
        ##elif theory.lower() == 'cnoidal':
        ##    WaveCnoidal.__init__(self)
        #else:
        #    self._wave = WaveFourier(Hw=Hw, Tw=Tw, Lw=Lw, d=d, title=name,
        #                             order=order, nstep=nstep, niter=niter,
        #                             accuracy=accuracy,
        #                             current=current, c_type=c_type,
        #                             infinite_depth=infinite_depth)
        ##
        #self.finite_depth = True
        #if infinite_depth:
        #    self.finite_depth = False
        ##

    #
    #def __call__(self, surface_points: int,
    #             depth_points: int):
    #    """ Solver """
    #    self._wave()
    #    #wave_length = (tau / self._z[ 1 ]) * self.d
    #    #
    #    # Calculate Surface
    #    surface = self._wave.get_surface(surface_points)
    #    #surface['wave_name'] = self.name
    #    surface['wave_id'] = self.number
    #    #        
    #    # Update wave lenght
    #    conn = create_connection(self.db_file)
    #    with conn:
    #        items = (self.number,)
    #        table = f"UPDATE WaveRegular \
    #                 SET \
    #                 wave_length = {self._wave.Lw}, \
    #                 wave_order = {self._wave.order}, \
    #                 wave_crest = {surface.eta.max()}, \
    #                 wave_trough = {surface.eta.min()}, \
    #                 wave_skewness = {surface.eta.max() / self._wave.H} \
    #                 WHERE wave_id = ?"
    #        cur = conn.cursor()
    #        cur.execute(table, items)
    #    #
    #    # update wave surface
    #    conn = create_connection(self.db_file)
    #    with conn:
    #        # push surface data
    #        surface.to_sql('WaveSurface', conn,
    #                       index_label=['wave_id', 'type', 'length', 'eta', 'phase', 'time'],
    #                       if_exists='append', index=False)
    #        #
    #        # get surface data from sql 
    #        surface = pull_surface(conn, number=self.number)
    #    #
    #    #
    #    kindf = self._wave.get_kinematics(depth_points=depth_points)
    #    #
    #    #kpoints = depth_points
    #    #kindf = get_kinematic(self.order, self._z, self._B, self._Tanh,
    #    #                      self.d, surface, kpoints,
    #    #                      self.finite_depth)
    #    #
    #    kindf['wave_id'] = self.number
    #    #kindf['type'] = 'order_1'
    #    #
    #    depth_steps = depth_steps = np.arange(depth_points + 1) / depth_points
    #    kindf['surface_id'] = repmat(surface['number'].to_numpy(),
    #                                 depth_steps.size, 1).flatten('F')
    #    #
    #    kindf.drop(columns=['phase', 'x', 'time'], inplace=True, axis=1)
    #    kindf.rename(columns={'z': 'elevation'}, inplace=True)
    #    #
    #    conn = create_connection(self.db_file)
    #    with conn:
    #        kindf.to_sql('WaveKinematic', conn,
    #                     index_label=['wave_id', 'surface_id', 'type',
    #                                  'elevation', 'u', 'v', 'dphidt',
    #                                  'ut', 'vt', 'ux', 'uz',
    #                                  'pressure', 'Benoulli_check'],
    #                     if_exists='append', index=False)
    #    #
    #    #print('---> here')
    #
    #
    @property
    def Lw(self):
        """ Wave Length"""
        conn = create_connection(self.db_file)
        with conn:
            data = pull_wave(conn, wave_name=self.name)
            wave_length = data[9]
        #
        #if not wave_length:
        #    #self.solver()
        #    #wave_length = (tau / self._z[ 1 ]) * self.d
        #    #conn = create_connection(self.db_file)
        #    # Update wave lenght
        #    #with conn:
        #    #    items = (self.name, )
        #    #    table = f"UPDATE Wave \
        #    #             SET wave_length = {wave_length} \
        #    #             WHERE name = ?"
        #    #    cur = conn.cursor()
        #    #    cur.execute(table, items)            
        #    self.__call__()
        #    data = self._pull_wave(conn)
        #    wave_length = data[9]
        #    #
        return wave_length

    #
    #@property
    #def wave(self):
    #    """ """
    #    conn = create_connection(self.db_file)
    #    with conn:
    #        data = pull_wave(conn, wave_name=self.name)
    #    #
    #    return RegWaveBasic(number=data[0],
    #                        name=data[1],
    #                        Hw=data[3],
    #                        Tw=data[4],
    #                        d=data[5],
    #                        theory=data[6])
    #
    #
    #
    def surface(self):
        """
        Free surface (wave) profile
        """
        conn = create_connection(self.db_file)
        #
        #df = DBframework()
        #header = ['wave_name', 'number', 'wave_id', 'type', 'x', 'eta', 'phase']
        #surface = df.DataFrame(data=data, columns=header)
        #surface = surface[['wave_name', 'type', 'x', 'eta', 'phase']]
        #
        #try:
        with conn:
            #(number, name, wave_type,
            # Hw, Tw, d, theory, order, Lw, WCe, WCtrough, Wskewness)            
            wdata = pull_wave(conn, wave_name=self.name)
            surface = pull_surface(conn, number=self.number)
        #surface.columns
        #except AttributeError:
        #    #surface = self.get_surface(surface_points)
        #    #surface =  surface[['wave_id', 'type', 'x', 'eta', 'phase']]
        #    #
        #    #conn = create_connection(self.db_file)
        #    #with conn:
        #    #    surface.to_sql('WaveSurface', conn,
        #    #                 index_label=['wave_id', 'type', 'x', 'eta', 'phase'], 
        #    #                 if_exists='append', index=False)
        #    #    surface['wave_name'] = self.name
        #    #
        #    self.__call__(surface_points=surface_points)
        #    with conn:
        #        surface = self._pull_surface(conn)
        #        #
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
        #surface['wave_id'] = self.number
        #surface['type'] = 'order_1'
        #
        #conn = create_connection(self.db_file)
        #with conn:
        #    surface.to_sql('WaveSurface', conn,
        #                 index_label=['wave_id', 'type', 'x', 'eta', 'phase'], 
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
        #wave = self.wave
        #
        return SurfaceResults(surface=surface,
                              Hw=wdata[3], Tw=wdata[4], d=wdata[5], 
                              finite_depth=False) # TODO : fix this

    #
    def kinematics(self):
        """ Get wave kinematics """
        #
        conn = create_connection(self.db_file)
        #
        #self.__call__()
        #try:
        with conn:
            kindf = pull_kinematics(conn, number=self.number)
        kindf.columns
        #except AttributeError:
        #    #surface = self.surface(surface_points)
        #    #
        #    #n = self.order
        #    #kd = self._z[1]
        #    #
        #    #crest = surface.eta
        #    #
        #    #depth = np.arange(depth_points + 1) / depth_points
        #    #
        #    #zdepth =  crestmax * depth  # --> fix sign
        #    #zdepth = surfacePoints(d=self.d, points=depth_points, eta=crest)
        #    #
        #    #get_kinematicX(n, self._z, self._Y, self._B, self._Tanh,
        #    #               surface_points, depth_points, self.finite_depth)
        #    #
        #    #1 / 0
        #    #kpoints = depth_points
        #    #kindf = get_kinematic(self.order, self._z, self._B, self._Tanh,
        #    #                      self.d, surface, kpoints,
        #    #                      self.finite_depth)
        #    #
        #    self.__call__(surface_points=surface_points,
        #                  depth_points=depth_points)
        #    with conn:
        #        kindf = self._pull_kinematics(conn)
        #        #
        #with conn:
        surface = self.surface()
        #
        #return KinematicResults(surface, kindf, depth_points+1)
        return KinematicResults(surface=surface,
                                kindata=kindf)
                                #depth_points=self._vpoints)
        #
        #return dataframe


#
def pull_surface(conn, number) -> list | None:
    """get wave surface"""
    #
    items = (number,)
    table = 'SELECT Wave.name, WaveSurface.* \
             FROM Wave, WaveSurface \
             WHERE WaveSurface.wave_id = ? \
             AND Wave.number = WaveSurface.wave_id'
    #try:
    cur = conn.cursor()
    cur.execute(table, items)
    data = cur.fetchall()
    #
    df = DBframework()
    header = ['wave_name', 'number', 'wave_id', 'type',
              'length', 'eta', 'phase', 'time']
    surface = df.DataFrame(data=data, columns=header)
    return surface[['number', 'wave_name', 'wave_id',
                        'type', 'length', 'eta', 'phase', 'time']]
#
#
def pull_wave(conn, wave_name:str|int,
              item: str = "*",
              wave_type:str='regular'):
    """ """
    items = (wave_name, wave_type)
    table = 'SELECT * FROM Wave \
             WHERE name = ? AND type = ?'
    cur = conn.cursor()
    cur.execute(table, items)
    wave = data = cur.fetchone()
    #return data
    project = (wave[0],)
    sql = f'SELECT {item} FROM WaveRegular \
            WHERE wave_id = ?'
    cur = conn.cursor()
    cur.execute(sql, project)
    data = cur.fetchone()
    return [*wave[:3], *data[2:]]
#
#
def pull_kinematics(conn, number:int):
    """read kin from sql"""
    items = (number,)
    table = 'SELECT Wave.name, WaveSurface.length, WaveSurface.phase, \
             WaveKinematic.* \
             FROM Wave, WaveSurface, WaveKinematic \
             WHERE WaveSurface.wave_id = ? \
             AND Wave.number = WaveKinematic.wave_id \
             AND WaveSurface.number = WaveKinematic.surface_id'
    #try:
    cur = conn.cursor()
    cur.execute(table, items)
    data = cur.fetchall()
    #
    df = DBframework()
    header = ['wave_name', 'length', 'phase', 'index',
              'wave_id', 'surface_id', 'type',
              'elevation', 'u', 'v', 'dphidt',
              'ut', 'vt', 'ux', 'uz',
              'pressure', 'Benoulli_check']
    surface = df.DataFrame(data=data, columns=header)
    #
    header = ['wave_name', 'wave_id', 'surface_id', 'type',
              'phase', 'length',
              'elevation', 'u', 'v', 'dphidt',
              'ut', 'vt', 'ux', 'uz',
              'pressure', 'Benoulli_check']
    return surface[header]
 #
 #