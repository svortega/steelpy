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
#import numpy as np


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
            #surface.rename(columns={'theta': 'length'}, inplace=True)
            # push surface data
            surface.to_sql('WaveSurface', conn,
                           index_label=['wave_id', 'type', 'theta',
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
        #depth_steps = np.arange(self._vpoints + 1) / self._vpoints
        kindf['surface_id'] = repmat(surface['number'].to_numpy(),
                                     self._vpoints, 1).flatten('F')
        #
        #kindf.drop(columns=['phase', 'length', 'time', 'theta'], inplace=True, axis=1)
        kindf.rename(columns={'z': 'elevation'}, inplace=True)
        #
        conn = create_connection(self.db_file)
        with conn:
            kindf.to_sql('WaveKinematic', conn,
                         index_label=['wave_id', 'surface_id', 'type',
                                      'elevation', 'u', 'v', 'dphidt',
                                      'ut', 'vt', 'ux', 'uz',
                                      'dudt', 'dvdt', 
                                      'pressure', 'Bernoulli_check'],
                         if_exists='append', index=False)
        #
        print('---> here')
    #
    #
    # ------------------------------
    # SQL ops
    #
    def _new_table(self, conn) -> None:
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
                    theta DECIMAL NOT NULL, \
                    eta DECIMAL NOT NULL, \
                    phase DECIMAL, \
                    time DECIMAL NOT NULL);"
        create_table(conn, table)
        #
        # Wave Kinematics
        table = "CREATE TABLE IF NOT EXISTS WaveKinematic (\
                    number INTEGER PRIMARY KEY NOT NULL,\
                    wave_id INTEGER NOT NULL REFERENCES Wave(number), \
                    surface_id INTEGER NOT NULL REFERENCES WaveSurface(number), \
                    type TEXT NOT NULL, \
                    elevation DECIMAL NOT NULL, \
                    u DECIMAL NOT NULL, \
                    v DECIMAL NOT NULL, \
                    dphidt DECIMAL NOT NULL, \
                    ut DECIMAL NOT NULL, \
                    vt DECIMAL NOT NULL, \
                    ux DECIMAL NOT NULL, \
                    uz DECIMAL NOT NULL, \
                    dudt DECIMAL NOT NULL, \
                    dvdt DECIMAL NOT NULL, \
                    pressure DECIMAL NOT NULL, \
                    Bernoulli_check DECIMAL NOT NULL);"
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
class RegWaveItem:
    __slots__ = ['number', 'name', 'db_file', '_wave',
                 'finite_depth', '_hpoints', '_vpoints',]

    def __init__(self, number: int, name: str | int,
                 db_file: str, infinite_depth: bool) -> None:
        """
        """
        self.number = number
        self.name = name
        self.db_file = db_file
    #
    #
    @property
    def Lw(self):
        """ Wave Length"""
        conn = create_connection(self.db_file)
        with conn:
            data = pull_wave(conn, wave_name=self.name)
            wave_length = data[9]
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
        with conn:
            #(number, name, wave_type,
            # Hw, Tw, d, theory, order, Lw, WCe, WCtrough, Wskewness)            
            wdata = pull_wave(conn, wave_name=self.name)
            surface = pull_surface(conn, number=self.number)
        #
        return SurfaceResults(surface=surface,
                              Hw=wdata[3], Tw=wdata[4], d=wdata[5], 
                              finite_depth=False) # TODO : fix this

    #
    def kinematics(self):
        """ Get wave kinematics """
        #
        conn = create_connection(self.db_file)
        with conn:
            kindf = pull_kinematics(conn, number=self.number)
        kindf.columns
        surface = self.surface()
        #return KinematicResults(surface, kindf, depth_points+1)
        return KinematicResults(surface=surface,
                                kindata=kindf)
                                #depth_points=self._vpoints)
#
#
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
    table = 'SELECT Wave.name, WaveSurface.theta, WaveSurface.phase, \
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
              'dudt', 'dvdt', 
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