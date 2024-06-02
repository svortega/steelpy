# 
# Copyright (c) 2019 steelpy
#
# 
# Python stdlib imports
from __future__ import annotations
from collections.abc import Mapping
#from collections import defaultdict
from dataclasses import dataclass
from typing import NamedTuple
#from math import prod

# package imports
from steelpy.utils.units.main import Units
#
from steelpy.metocean.process.bsotm import BSOTM
from steelpy.metocean.process.parameters import HydroDesign
from steelpy.metocean.wave.regular.main import RegWaveItem
from steelpy.metocean.hydrodynamic.wkf import WKFitem
from steelpy.metocean.current.main import CurrentItem
from steelpy.metocean.hydrodynamic.cbf import CBFitem
from steelpy.metocean.hydrodynamic.marine_growth import MGitem
from steelpy.metocean.hydrodynamic.morison import CdCmitem
from steelpy.metocean.hydrodynamic.element_segment import WIPitem
#
from steelpy.utils.sqlite.utils import create_connection, create_table
#
from steelpy.utils.dataframe.main import DBframework
#
#
#
class HydroCondition(Mapping):
    """
    FE Metocean Combination Class
    
    Combination
        |_ name
        |_ number
        |_ surface
        |_ gravity
        |_ water_density
        |_ air_density
        |_ fluid_viscosity
        |_ air_viscosity
        |_ zone
    
    **Parameter**:  
      :number:  integer internal number 
      :name:  string node external name
    """
    __slots__ = ['db_file', '_properties', '_design', '_criteria']

    def __init__(self, criteria: str,
                 properties, db_file: str):
        """
        """
        self._criteria = criteria
        self.db_file = db_file
        self._properties = properties
        self._design = HydroDesign(self.db_file)
        #
        # create table
        conn = create_connection(self.db_file)
        with conn:
            self._create_table(conn)
    #
    @property
    def _labels(self):
        """ """
        query = (self._criteria, )
        table = 'SELECT Condition.name FROM Condition \
                 WHERE criteria_id = ?'
        conn = create_connection(self.db_file)
        with conn:        
            cur = conn.cursor()
            cur.execute(table, query)
            items = cur.fetchall()
        return [item[0] for item in items]    
    #
    #
    def __setitem__(self, comb_name: str|int, title: str) -> None:
        """
        title, marine_growth, Buoyancy(False/True), conductor_shielding
        """
        #
        #values = self.get_setup(setupval)
        #comb_title = values.pop(0)
        #self._combination[comb_name] = CombTypes(comb_name, comb_title, self)
        #self._combination[comb_name].setup = values
        #
        #
        comb_data = [comb_name, title, self._criteria,   # name, title, criteria, 
                     0, 0, None, 'off',  # wave_id, wave_direction, wave_kfnumber, buoyancy
                     0, 0, None, 'on',   # current_id, current_direction, current_bfnumber, stretching,
                     0, 0, None]         # wind_number, wind_direction, wind_area_number
                     #None, None]         # mg_id, cdcm_id
        conn = create_connection(self.db_file)
        with conn:
            self._push_data(conn, comb_data)        
    #
    def __getitem__(self, comb_name:str|int):
        """
        """
        return CombTypes(comb_name,
                         criteria=self._criteria, 
                         properties=self._properties,
                         design=self._design, 
                         db_file=self.db_file)
        #1 / 0
        #return self._combination[comb_name]
    #
    #def __delitem__(self, load_name:str|int):
    #    """
    #    """
    #    del self._combination[load_name]
    #
    def __len__(self) -> int:
        return len(self._labels)
    #
    def __iter__(self):
        """
        """
        return iter(self._labels)
    #
    #
    #def get_setup(self, values:str|list|dict):
    #    """
    #    [title, marine_growth, CdCm, Buoyancy]
    #    """
    #    output = [None, False, False, False]
    #
    #    if isinstance(values, str):
    #        output[0] = values
    #
    #    elif isinstance(values, dict):
    #        raise NotImplementedError
    #
    #    elif isinstance(values, list):
    #        for idx, item in enumerate(values):
    #            output[idx] = item
    #
    #    else:
    #        raise IOError("input not valid")
    #
    #    return output
    #
    # ----------------------------------------
    # SQL ops
    #
    def _create_table(self, conn) -> None:
        """ """
        table = "CREATE TABLE IF NOT EXISTS Condition (\
                    number INTEGER PRIMARY KEY NOT NULL,\
                    name NOT NULL,\
                    title TEXT NOT NULL,\
                    criteria_id INTEGER NOT NULL REFERENCES Criteria(number),\
                    wave_id INTEGER NOT NULL REFERENCES Wave(number),\
                    wave_direction  DECIMAL NOT NULL,\
                    wave_kf_id INTEGER REFERENCES WaveKinFactor(number),\
                    buoyancy TEXT NOT NULL,\
                    current_id INTEGER REFERENCES Current(number),\
                    current_direction DECIMAL,\
                    current_bf_id INTEGER REFERENCES CurrentBlockageFactor(number),\
                    stretching TEXT,\
                    wind_id INTEGER REFERENCES Wind(number),\
                    wind_direction DECIMAL,\
                    wind_area_id INTEGER REFERENCES WindArea(number));"
        create_table(conn, table)
    #
    #
    def _push_data(self, conn, comb_data):
        """
        Create a new project into the projects table
        """
        #
        #cur = conn.cursor()
        #cur.execute("SELECT Current.name, Current.number FROM Current;")
        #items = cur.fetchall()
        #current = {item[0]:item[1] for item in items}        
        #
        #
        cur = conn.cursor()
        table = 'INSERT INTO Condition(name, title, criteria_id, \
                                            wave_id, wave_direction, wave_kf_id, buoyancy, \
                                            current_id, current_direction, current_bf_id, stretching, \
                                            wind_id, wind_direction, wind_area_id) \
                                            VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?)'
        # push
        cur = conn.cursor()
        cur.execute(table, tuple(comb_data))
        #print('-->')
    #
    #
    #
    # --------------------------------------
    #
    def df(self, df):
        """ """
        #1 / 0
        #grpname = df.groupby(['name'])
        dfcond = df.set_index('name')
        #
        for item in dfcond.itertuples():
        #for key, item in grpname:
            name = item.Index
            #
            comb = item.title #[['title']].values.tolist()
            #                  'conductor_shielding']].values.tolist()
            #
            wave = [item.wave_name, item.wave_direction,
                    item.wave_kinematics, item.crest_elevation]            
            #
            #wave = item[['wave_name', 'wave_direction',
            #             'wave_kinematics', 'crest_elevation']].values.tolist()
            #
            #
            current = [item.current_name, item.current_direction,
                       item.current_blockage, item.current_stretching]
            #
            #current = item[['current_name', 'current_direction',
            #                'current_blockage', 'current_stretching']].values.tolist()
            #
            try:
                wind = [item.wind_name, item.wind_direction]
            except AttributeError:
                wind = [None, None]
            #
            #wind = item[['wind_name', 'wind_direction']].values.tolist()
            #
            #
            properties = [item.MG, item.CdCm, item.conductor_shielding]         
            #
            #properties = item[['MG', 'CdCm', #'flooding',
            #                   'conductor_shielding']].values.tolist()
            #
            #
            parameters = [item.design_load, item.buoyancy, item.criterion]
            #
            #parameters = item[['design_load', 'buoyancy', 'criterion']].values.tolist()            
            #
            self.__setitem__(comb_name=name, title=comb)
            _combination = self.__getitem__(comb_name=name)
            _combination.wave = wave
            _combination.current = current
            #_combination.wind = wind
            _combination.properties = properties
            _combination.parameters = parameters
        #print('--')
        #1 / 0
    #
#
#
class ConditionBasic(NamedTuple):
    """
    """
    number:int
    name:str|int
    title: str
    criteria: str | int
    #
    wave_id:str|int = 0
    wave_direction:float = 0.0
    wave_kfnumber:float = 1.0
    #
    buoyancy : bool = True
    #
    current_id: str|int = 0
    current_direction:float = 0.0
    current_bfnumber:float  = 1.0
    current_stretching : bool = True
    #
    wind_number:str|int = 0
    wind_direction:float = 0.0
    wind_area_number: int = 0
#
class ParametersBasic(NamedTuple):
    """
    """
    number: int
    name: str|int
    title: str
    #
    design_load: str
    buoyancy : str|None
    criterion: str
    scaling_factor: float
#
class WaveBasic(NamedTuple):
    #number:int
    wave:str|int
    direction:float
    kinematic_factor:list
    crest_elevation:float|None
    #title:str
    #
    def kinematics(self):
        """vave kinemactis"""
        return self.wave.kinematics()
    #    print('here')
#
class WindBasic(NamedTuple):
    #number:int
    wind:str|int
    direction:float
    #wind_areas:list = []
    title:str
#
class CurrentBasic(NamedTuple):
    #number:int
    current: str|int
    direction:float
    blockage_factor:float
    stretching : bool
    #title:str
#
class PropertyBasic(NamedTuple):
    """ """
    number:int
    name:str|int
    title: str
    rho_w: float
    #
    marine_growth: tuple
    CdCm: tuple
    #WKF: tuple
    #conductor_shielding: tuple | None = None
    WIP: tuple
#
#
#
@dataclass
class CombTypes:
    """
    """
    __slots__ = ['name', 'number', '_db_file',
                 '_properties', '_design', '_criteria']
                # '_wave', '_current', '_wind', '_setup', 'number', 'title',
    
    def __init__(self, name:int|str,
                 criteria: str|int, 
                 properties, design, 
                 db_file: str):
        """
        """
        self.name = name
        self._criteria = criteria
        self._db_file = db_file
        self._properties= properties
        self._design = design
        #
        condition = self.condition
        self.number = condition.number
    #
    # -----------------------------
    #
    @property
    def condition(self):
        """ """
        conn = create_connection(self._db_file)
        with conn:        
            items = self._pull_condition(conn)
        #
        return ConditionBasic._make(items)
    #
    def _pull_condition(self, conn, item: str = '*'):
        """ """
        project = (self.name, self._criteria)
        table = f'SELECT {item} FROM Condition \
                  WHERE name = ? AND criteria_id = ?'
        cur = conn.cursor()
        cur.execute(table, project)
        condition = cur.fetchone()
        #
        #wave = self._pull_wave(conn, name=data[0])
        #
        #current = self._pull_current(conn, name=data[1])
        if not condition:
            raise IOError(f'Condition {self.name} not valid')
        #
        return condition
    #
    #
    @property
    def title(self):
        """ """
        conn = create_connection(self._db_file)
        with conn:        
            items = self._pull_condition(conn, 'title')
        return items[0]
    #
    #
    # -----------------------------
    #
    @property
    def setup(self):
        """ comb setup"""
        
        return self._setup

    @setup.setter
    def setup(self, values=list):
        """ comb setup"""
        # marine_growth
        try:
            mg = self._cls._marine_growth[values[0]]
        except KeyError:
            units = Units()
            self._cls._marine_growth['default'] = [self._cls._hydro.rho_w * units.kg / units.m**3,
                                                   0 * units.mm, 'constant']
            mg = self._cls._marine_growth['default']
        # CdCm
        #try:
        CdCm = self._cls._cdcm
        #except KeyError:
        #    1/0
        #
        buoyancy = values[2]
        self._setup = CombSetup(marine_growth=mg,
                                CdCm=CdCm,
                                hydro=self._cls._hydro,
                                buoyancy=buoyancy)
    #
    # -----------------------------
    #
    @property
    def wave(self):
        """
        """
        cond = self.condition
        conn = create_connection(self._db_file)
        with conn:
            # Wave data
            project = (cond.wave_id, self._criteria)
            table = 'SELECT * FROM Wave \
                     WHERE number = ? AND criteria_id = ?'
            cur = conn.cursor()
            cur.execute(table, project)      
            wdata = cur.fetchone()            
            #
            # Wave kinemactis factor
            #with conn:
            project = (cond.wave_kfnumber,)
            table = 'SELECT * FROM WaveKinFactor WHERE number = ?'
            cur = conn.cursor()
            cur.execute(table, project)
            wkfdata = cur.fetchone()            
            wkf = WKFitem(name=wkfdata[1], db_file=self._db_file)
        #
        # Regular wave
        if wdata[2].lower() in ['regular']:
            # get wave data
            with conn:
                project = (wdata[0],)
                table = 'SELECT * FROM WaveRegular WHERE wave_id = ?'
                cur = conn.cursor()
                cur.execute(table, project)
                wregdata = cur.fetchone()
            # TODO : WCe tb included 
            wave = RegWaveItem(number=wdata[0], name=wdata[1],
                               Hw=wregdata[2], Tw=wregdata[3], d=wregdata[4], 
                               theory=wregdata[5], order=wregdata[6], 
                               Lw=wregdata[7], db_file=self._db_file)
            #
            return WaveBasic(wave=wave, 
                             direction = cond.wave_direction, 
                             kinematic_factor=wkf, 
                             crest_elevation=wregdata[8])
        else:
            raise NotImplementedError(f'wave type {wave[2]} not yet implemented')
    
    @wave.setter
    def wave(self, values: list|dict):
        """
        [wave_name, Direction(deg), Kinematics, title]
        """
        values = self._get_wvalues(values)
        conn = create_connection(self._db_file)
        with conn:
            wave = self._pull_wave(conn, wave_name=values[0])
            if not wave:
                raise IOError(f'Wave {values[0]} not valid')
            #
            wkf = self._pull_wkf(conn, name=values[2])
            if not wkf:
                print('fix wkf')
                1 / 0
            #
            query = (self.name, self._criteria, )
            #
            table = f"UPDATE Condition \
                     SET wave_id = {wave[0]}, \
                     wave_direction= {values[1]}, \
                     wave_kf_id = {wkf[0]}, \
                     buoyancy = '{values[3]}' \
                     WHERE name = ? \
                     AND criteria_id = ?"
            #
            cur = conn.cursor()
            cur.execute(table, query)
        #print('-->')
    #
    def _pull_wave(self, conn, wave_name: int|str,
                  item: str = '*'):
        """ """
        #cur = conn.cursor()
        #cur.execute("SELECT * FROM Wave;")
        project = (wave_name, self._criteria)
        table = f'SELECT {item} FROM Wave \
                 WHERE name = ? AND criteria_id = ?'
        cur = conn.cursor()
        cur.execute(table, project)      
        items = cur.fetchone()
        #wave = {item[0]:item[1] for item in items}
        return items
    #
    #
    def _pull_wkf(self, conn, name: int|str,
                  item: str = '*'):
        """ """
        project = (name,)
        table = f'SELECT {item} FROM WaveKinFactor \
                  WHERE name = ?'
        cur = conn.cursor()
        cur.execute(table, project)      
        items = cur.fetchone()
        return items
    #
    #
    @property
    def current(self):
        """
        """
        cond = self.condition
        #
        conn = create_connection(self._db_file)
        with conn:
            project = (cond.current_id, self._criteria)
            table = 'SELECT * FROM Current \
                     WHERE number = ? AND criteria_id = ?'
            cur = conn.cursor()
            cur.execute(table, project)      
            currdata = cur.fetchone()            
        #
        #
        # Current blockage factor
        with conn:
            project = (cond.current_bfnumber,)
            table = 'SELECT * FROM CurrentBlockageFactor WHERE number = ?'
            cur = conn.cursor()
            cur.execute(table, project)
            cbfdata = cur.fetchone()            
            cbf = CBFitem(name=cbfdata[1], db_file=self._db_file)        
        #
        curritem = CurrentItem(name=currdata[1],
                               db_file=self._db_file,
                               criteria=self._criteria)
        #
        #1 / 0
        return CurrentBasic(current=curritem,
                            direction=cond.current_direction,
                            blockage_factor=cbf,
                            stretching=cond.current_stretching)
                            #title=title)
    
    @current.setter
    def current(self, values: list|dict):
        """
        """
        values = self._get_cvalues(values)
        #try:
        conn = create_connection(self._db_file)
        with conn:            
            current = self._pull_current(conn, current_name=values[0])
            if not current:
                raise IOError(f'Current {values[0]} not valid')
            #
            cbf = self._pull_cbf(conn, name=values[2])
            if not cbf:
                print('fix cbf')
                1 / 0                
            #
            query = (self.name, self._criteria)
            #
            table = f"UPDATE Condition \
                     SET current_id = {current[0]}, \
                     current_direction= {values[1]}, \
                     current_bf_id = {cbf[0]}, \
                     stretching = '{values[3]}' \
                     WHERE name = ? \
                     AND criteria_id = ?"
            #
            cur = conn.cursor()
            cur.execute(table, query)
        #except KeyError:
        #    raise IOError(f'current {values[0]} not found')
        #
        
    #
    #
    def _pull_current(self, conn, current_name: int|str,
                  item: str = '*'):
        """ """
        query = (current_name, self._criteria)
        table = f'SELECT {item} FROM Current \
                 WHERE name = ? AND criteria_id = ?'
        cur = conn.cursor()
        cur.execute(table, query)
        items = cur.fetchone()
        #wave = {item[0]:item[1] for item in items}
        return items
    #
    def _pull_cbf(self, conn, name: int|str,
                  item: str = '*'):
        """ Current Blockage Factor"""
        project = (name,)
        table = f'SELECT {item} FROM CurrentBlockageFactor \
                 WHERE name = ?'
        cur = conn.cursor()
        cur.execute(table, project)      
        items = cur.fetchone()
        return items
    #
    #
    @property
    def wind(self):
        """
        """
        return self._wind
    
    @wind.setter
    def wind(self, values: list|dict):
        """
        """
        1 / 0
        values, title = self.get_values(values)
        self._wind = WindBasic(wind=values[0],
                               direction=values[1],
                               title=title)
    #
    #
    # -----------------------------
    # Hydro parameters
    #
    @property
    def properties(self):
        """ cdcm class"""
        #
        conn = create_connection(self._db_file)
        with conn:
            propitem = self._pull_property(conn)
            # rho
            project = (propitem[0],)
            table = 'SELECT Criteria.rho_water \
                     FROM Criteria, Condition  \
                     WHERE Condition.number = ? \
                     AND Condition.criteria_id = Criteria.number'
            cur = conn.cursor()
            cur.execute(table, project)      
            item = cur.fetchone()
            rho = float(item[0])
            #
            # MG
            #with conn:
            project = (propitem[3],)
            table = 'SELECT * FROM MarineGrowth WHERE number = ?'
            cur = conn.cursor()
            cur.execute(table, project)      
            mgdata = cur.fetchone()       
            mg = MGitem(name=mgdata[1], db_file=self._db_file)
            #
            # CdCm
            #with conn:
            project = (propitem[4],)        
            table = 'SELECT * FROM CdCm WHERE number = ?'
            cur = conn.cursor()
            cur.execute(table, project)      
            cdcmdata = cur.fetchone()        
            cdcm = CdCmitem(name=cdcmdata[1], db_file=self._db_file)
            #
            # element segmentation
            project = (propitem[5],)
            table = 'SELECT * FROM ElementSegment WHERE number = ?'
            cur = conn.cursor()
            cur.execute(table, project)      
            wipdata = cur.fetchone()
            wip = WIPitem(name=wipdata[1], db_file=self._db_file)
        #
        return PropertyBasic(number=propitem[0], 
                             name=propitem[1],
                             title=propitem[2],
                             rho_w=rho, 
                             marine_growth=mg, 
                             CdCm=cdcm,
                             WIP=wip)
                             #flooding=self._properties.flooding)
    
    @properties.setter
    def properties(self, values: list|tuple|dict):
        """ cdcm input"""
        
        if isinstance(values, list|tuple|dict):
            data = self._get_pdata(values)
        else:
            raise IOError('Input data not valid')
        #
        conn = create_connection(self._db_file)
        with conn:
            self._push_property(conn, data)
    #
    def _push_property(self, conn, data: tuple,
                       item: str = "*"):
        """ """
        #
        condition = self._pull_condition(conn)
        #
        MG = self._pull_mg(conn, name=data[0])
        #
        CdCm = self._pull_cdcm(conn, name=data[1])
        #
        WIP = self._pull_wip(conn, name=data[2])
        #
        profile = (condition[0], MG[0], CdCm[0], WIP[0], *data[3:])
        cur = conn.cursor()
        table = 'INSERT INTO Property(condition_id, \
                                      mg_id, cdcm_id, \
                                      element_segment_id, \
                                      flooding_id, cshielding_id, \
                                      title) \
                                      VALUES(?,?,?,?,?,?,?)'
        #
        # push
        cur = conn.cursor()
        cur.execute(table, profile)        
        #       
    #
    def _pull_property(self, conn, item: str = '*'):
        """ """
        # Condition
        cond = self.condition
        #
        project = (cond.number,)
        table = f'SELECT {item} FROM Property \
                  WHERE condition_id = ?'
        cur = conn.cursor()
        cur.execute(table, project)      
        prop = cur.fetchone()
        return [*cond[:3], *prop[2:-1]]
    #
    #
    def _pull_cdcm(self, conn, name: int|str,
                   item: str = '*'):
        """ """
        project = (name,)
        table = f'SELECT {item} FROM CdCm WHERE name = ?'
        cur = conn.cursor()
        cur.execute(table, project)      
        items = cur.fetchone()
        return items    
    #
    def _pull_mg(self, conn, name: int|str,
                 item: str = '*'):
        """ """
        project = (name,)
        table = f'SELECT {item} FROM MarineGrowth WHERE name = ?'
        cur = conn.cursor()
        cur.execute(table, project)      
        items = cur.fetchone()
        return items
    #
    def _pull_wip(self, conn, name: int|str,
                 item: str = '*'):
        """ """
        project = (name,)
        table = f'SELECT {item} FROM ElementSegment WHERE name = ?'
        cur = conn.cursor()
        cur.execute(table, project)      
        items = cur.fetchone()
        return items
    #
    # -----------------------------
    # Seastate parameters
    #
    @property
    def parameters(self):
        """ """
        conn = create_connection(self._db_file)
        with conn:
            param = self._pull_parameters(conn)
        #
        #return self._design
        return ParametersBasic._make(param)
    
    @parameters.setter
    def parameters(self, values: list|tuple|dict):
        """ """
        if isinstance(values, list|tuple|dict):
            data = self._get_ddata(values)
        else:
            raise IOError('Input data not valid')
        #
        conn = create_connection(self._db_file)
        with conn:
            self._push_parameters(conn, data)
    #
    def _push_parameters(self, conn, data: tuple,
                       item: str = "*"):
        """ """
        condition = self._pull_condition(conn)
        #
        profile = (condition[0], *data)
        cur = conn.cursor()
        table = 'INSERT INTO Parameter(condition_id, \
                                        design_load, buoyancy, \
                                        criterion, \
                                        scaling_factor, title) \
                 VALUES(?,?,?,?,?,?)'
        # push
        cur = conn.cursor()
        cur.execute(table, profile)
    #
    def _pull_parameters(self, conn, item: str = '*'):
        """ """
        # Condition
        cond = self.condition
        #
        project = (cond.number,)
        table = 'SELECT {:} FROM Parameter WHERE condition_id = ?'.format(item)
        cur = conn.cursor()
        cur.execute(table, project)      
        prop = cur.fetchone()
        #
        return [*cond[:3], *prop[2:-1]]
    #
    # -----------------------------
    # operations
    #
    def _get_wvalues(self, values):
        """ values : [wave_name, Direction(deg), Kinematics, Buoyancy(False/True)] """
        # TODO : values need update
        if isinstance(values, (list|tuple)):
            #
            # Direction
            try:
                values[1] = values[1].value
            except (AttributeError, IndexError):
                #values.append(False)
                pass
            #
            # Kinematics
            try:
                values[2]
            except IndexError:
                values.append(1.0)
            #
            # Buoyancy(False/True)
            try:
                if values[3]:
                    values[3] = 'on'
                else:
                    values[3] = 'off'
            except IndexError:
                values.append('off')
            #
            # crest_elevation
            #try:
            #    values[4] = values[4].value
            #except (AttributeError, IndexError):
            #    values.append(None)
            
        elif isinstance(values, dict):
            1 / 0
            
        else:
            raise IOError('Input data not valid')
        #
        return values
    #
    def _get_cvalues(self, values):
        """ [current_name,  Direction(deg), Blockage, Stretching] """
        if isinstance(values, (list|tuple)):
            #
            # Direction
            try:
                values[1] = values[1].value
            except (AttributeError, IndexError):
                #values.append(False)
                pass
            #
            # Blockage
            try:
                values[2]
            except IndexError:
                values.append(1.0)
            #
            # Stretching
            try:
                if values[3]:
                    values[3] = 'on'
                else:
                    values[3] = 'off'
            except IndexError:
                values.append('off')
            #
            
        elif isinstance(values, dict):
            1 / 0
            
        else:
            raise IOError('Input data not valid')
        #
        return values        
    #
    def _get_pdata(self, values):
        """ [marine_growth, CdCm, element_refinament, Flooding, conductor_shielding, title]"""
        #
        output = [None] * 6
        #
        if isinstance(values, (list|tuple)):
            for x, item in enumerate(values):
                output[x] = item
        
        elif isinstance(values, dict):
            1 / 0
        
        else:
            raise IOError('Input data not valid')
        #
        return output
    #
    def _get_ddata(self, values):
        """ [design_load, buoyancy, criterion, scaling_factor, title]"""
        #
        output = ['max_bs', None, 'local', 1.0, None]
        #
        if isinstance(values, (list|tuple)):
            for x, item in enumerate(values):
                output[x] = item
        
        elif isinstance(values, dict):
            1 / 0
        
        else:
            raise IOError('Input data not valid')
        #
        return output
    #
    # -----------------------------
    # load process
    #
    def load(self):
        """ convert to beam load"""
        wave = self.wave
        current = self.current
        hydro = self.properties
        rho_w = self.properties.rho_w
        #
        res = BSOTM(wave=wave, 
                    current=current,
                    properties=hydro,
                    rho_w=rho_w, 
                    condition=2)
        return res
    #
    #
    def beam_load(self, beams):
        """
        Generate beam loading
        beam : beam class
        """
        # get load
        wload = self.load()
        # get parameters
        parameters = self.parameters
        #
        dftemp = []
        #
        # ------------------------------------------
        # TODO: Multiprocess --> sqlite issue
        #
        #cpuno = os.cpu_count() - 1
        #
        #def myfunc(beam):
        #    beam = BeamItemSQL(name,
        #                       plane=self._plane,
        #                       db_file=self.db_file)
        #    Fwave = wload.Fwave(beam=beam)
        #    return Fwave.solve()
        #
        #beams = [BeamItemWave(name,
        #                      db_file=self.db_file)
        #         for name in labels]
        #
        #with ProcessPoolExecutor(max_workers=cpuno) as executor:
        #    for r in executor.map(myfunc, beams):
        #        dftemp.append(r)
        #
        #engine = Engine(wload=wload)
        #                plane=self._plane,
        #                wload=wload)
        #
        #dftemp = engine.run(beams)
        #
        #with ProcessPoolExecutor(max_workers=cpuno) as executor:
        #    for r in executor.map(engine, beams):
        #        dftemp.append(r)
        #
        #try:
        #    pool = Pool(cpuno)
        #    engine = Engine(db_file=self.db_file,
        #                    plane=self._plane,
        #                    wload=wload)
        #    dftemp = pool.map(engine, beams)
        #finally:
        #    pool.close()
        #    pool.join()
        #
        # ------------------------------------------
        #
        for beam in beams.values():
            # solve beam forces
            Fwave = wload.Fwave(beam=beam)
            solution = Fwave.solve()
            #dftemp.extend(solution)
            dftemp.append(solution)
            #print('-->')
        #
        # ------------------------------------------
        # create database
        #
        #header = ['element_type', 'element_name', 'element_id', 
        #          'type',
        #          'qx0', 'qy0', 'qz0', 'qt0', 
        #          'qx1', 'qy1', 'qz1', 'qt1',
        #          'L0', 'L1', 'BS', 'OTM', 
        #          'x', 'y', 'z']        
        #
        df = DBframework()
        #df_bload = df.DataFrame(data=dftemp, columns=header, index=None)
        df_bload = df.concat(dftemp)
        #
        #
        # ------------------------------------------
        # select data in database        
        #
        design_load = parameters.design_load
        design_load = design_load.split("_")
        value_type = design_load[0]
        load_type = design_load[1].upper()
        
        if parameters.criterion == 'local':
            #blgrp = df_bload.groupby(['element_name', 'element_id', 'element_type'])
            #for key, wload in blgrp:
            #    key, wload
            #
            #
            grpwave = df_bload.groupby(['element_name','time'])[['BS', 'OTM']].sum()
            #
            grpm = grpwave[load_type].abs().groupby('element_name')
            #
            if value_type.lower() == 'max':
                idx = grpm.idxmax()
                design = 'local_max'
            else:
                1 / 0
            #
            df_bload.set_index(['element_name','time'], inplace=True)
            df_bload = df_bload.loc[idx]
            df_bload.reset_index(inplace=True)
        else:
            raise NotImplementedError
        #
        #
        # ------------------------------------------
        # update database
        #
        header = ['element_id', 'design_load', 
                  'title', 'system', 'type',
                  'L0', 'qx0', 'qy0', 'qz0', 'qt0',
                  'L1', 'qx1', 'qy1', 'qz1', 'qt1',
                  'BS', 'OTM', 'time'] # , 'x', 'y', 'z'
        #
        # FIXME : load_id
        df_bload['design_load'] = design
        df_bload['title'] = f'MET_{self.name}'
        # FIXME : global or local? 
        df_bload['system'] = 'local'
        #df_bload.rename(columns={'load_type': 'type'}, inplace=True)
        df_bload =  df_bload[header]            
        #
        # ------------------------------------------
        # push data in database        
        #
        #with conn:
        #    df_bload.to_sql('LoadBeamLine', conn,
        #                    index_label=header, 
        #                    if_exists='append', index=False)
        #1 / 0
        return df_bload
    #    
#