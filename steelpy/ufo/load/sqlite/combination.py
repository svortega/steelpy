# 
# Copyright (c) 2009 steelpy
#
# 
# Python stdlib imports
from __future__ import annotations
from collections.abc import Mapping
#from collections import defaultdict
from collections import defaultdict
#from dataclasses import dataclass
#from typing import NamedTuple
#from math import prod

# package imports
from steelpy.ufo.load.sqlite.load_case import BasicLoadSQL
#from steelpy.ufo.load.process.utils import duplicates, indices
from steelpy.ufo.load.process.combination import LoadCombinationBasic, LoadCombinationRoot
from steelpy.ufo.load.sqlite.beam import BeamLoadItemSQL, BeamLoadGloabalSQL
# steelpy
from steelpy.utils.sqlite.utils import create_connection, create_table
from steelpy.utils.dataframe.main import DBframework
#
#
class LoadCombSQL(LoadCombinationBasic):
    """Load Combination Class SQL"""
    __slots__ = ['_title', '_number', 'db_file', '_name',
                 '_index', '_basic', '_combination', '_metocean',
                 '_mesh_id']
    #
    def __init__(self, db_file:str, #plane: NamedTuple,
                 mesh_id: int, name:int|str):
        """
        """
        super().__init__(name)
        #
        self.db_file = db_file
        #self._name = name
        self._mesh_id = mesh_id
        #
        self._basic = BasicLoadSQL(db_file=self.db_file,
                                   name=name,
                                   mesh_id=mesh_id)
        self._combination = {}
        #
        conn = create_connection(self.db_file)
        with conn: 
            self._new_table(conn)
    #   
    #
    def __setitem__(self, load_name:int|str, load_title:str) -> None:
        """
        """
        try:
            self._labels.index(load_name)
            raise Exception('    *** warning load combination name {:} already exist'
                            .format(load_name))
        except ValueError:
            self._labels.append(load_name)
            self._title.append(load_title)
            load_id = next(self.get_number())
            self._number.append(load_id)
            #
            self._combination[load_name] = CombTypeSQL(name=load_name,
                                                       title=load_title,
                                                       mesh_id=self._mesh_id, 
                                                       db_file=self.db_file)
            #
            conn = create_connection(self.db_file)
            with conn:
                load_no = self._push_combination(conn, load_name,
                                                 self._mesh_id, 
                                                 'combination', 
                                                 load_title)
                
                self._combination[load_name].number = load_no
    #
    def __getitem__(self, load_name:str|int):
        """
        """
        try:
            self._index = self._labels.index(load_name)
            return self._combination[load_name]
        except ValueError:
            raise KeyError(f"load combination {load_name} not defined")
    #
    def _new_table(self, conn):
        """ """
        table = "CREATE TABLE IF NOT EXISTS LoadCombination(\
                    number INTEGER PRIMARY KEY NOT NULL,\
                    load_id INTEGER NOT NULL REFERENCES Load(number),\
                    basic_id INTEGER REFERENCES Load(number),\
                    combination_id INTEGER REFERENCES Load(number),\
                    factor DECIMAL NOT NULL);"
        create_table(conn, table)
    #
    def _push_combination(self, conn, name:int|str,
                          mesh_id:int, level:str, title:str):
        """
        """
        project = (name, mesh_id, level, title, 'ufo')
        sql = ('INSERT INTO Load(name, mesh_id, level,\
                                 title, input_type)\
               VALUES(?,?,?,?,?)')
        #
        cur = conn.cursor()
        cur.execute(sql, project)
        #conn.commit()
        return cur.lastrowid
        

    #def solve_combinations(self, basic_res, memb_force):
    #    """
    #    """
    #    comb_res = {}
    #    memb_comb = {}
    #    bloads = self.to_basic()
    #    for lcomb, comb in bloads.items():
    #        #lcomb = load_combination[cname].title
    #        beam_load = {}
    #        for bname, factor in comb:
    #            try:
    #                comb_res[lcomb] += basic_res[bname] * factor
    #            except KeyError:
    #                comb_res[lcomb] = basic_res[bname] * factor
    #            #
    #            for mname, member in memb_force[bname].items():
    #                try:
    #                    beam_load[mname] += member * factor
    #                except KeyError:
    #                    beam_load[mname] =  member * factor
    #        memb_comb[lcomb] = beam_load
    #    return comb_res, memb_comb
    #
    # -----------------------------------------------
    #
    def Pn(self):
        """
        Returns:
            Global nodal force vector
        """                             
        nload = self._basic._nodes.force
        values = ['Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']
        P = self._get_combination(values, nload)
        return P
    #    
    #
    #
    def to_basic(self):
        """ """
        comb2basic = combination_basic(mesh_name=self._mesh_id,
                                       db_file=self.db_file)
        return comb2basic
#     
#
class CombTypeSQL:
    """
    Load Combination Type Class
    """
    __slots__ = ['name', 'number', 'title', 
                 '_basic', '_metocean', '_combination',
                 'db_file', '_mesh_id', '_beams']

    def __init__(self, name:str, title:str,
                 mesh_id: int, db_file:str):
        """
        """
        self._mesh_id = mesh_id
        self.name = name
        self.title = title
        self.db_file = db_file
        #
        self._basic = BasicCombSQL(bl_type="basic",
                                   comb_name=self.name,
                                   mesh_id=self._mesh_id, 
                                   db_file=self.db_file)
        
        self._metocean = BasicCombSQL(bl_type="metocean",
                                      comb_name=self.name,
                                      mesh_id=self._mesh_id, 
                                      db_file=self.db_file)
        
        self._combination = BasicCombSQL(bl_type="combination",
                                         comb_name=self.name,
                                         mesh_id=self._mesh_id, 
                                         db_file=self.db_file)
        #
        self._beams = BeamLoadGloabalSQL(mesh_id=self._mesh_id,
                                         db_file=self.db_file)
    #
    #@property
    #def _labels(self):
    #    """ """
    #    query = ('combination', self._mesh_id, )
    #    table = "SELECT Load.name \
    #              FROM Load, LoadCombination \
    #              WHERE Load.level = ? \
    #              AND LoadCombination.load_id = Load.number\
    #              AND Load.mesh_id = ? ;"
    #    conn = create_connection(self.db_file)
    #    with conn:        
    #        cur = conn.cursor()
    #        cur.execute(table, query)
    #        items = cur.fetchall()
    #    items = set([item[0] for item in items])
    #    return list(items)
    #
    # ----------------------------    
    #
    @property
    def combination(self):
        """
        """
        return self._combination

    #
    @property
    def basic(self):
        """
        """
        return self._basic
    
    #
    @property
    def metocean(self):
        """
        """
        return self._metocean
    #
    @property
    def mesh_name(self):
        """ mesh name """
        query = (self._mesh_id, )
        table = 'SELECT name \
                 FROM Mesh WHERE number = ?;'
        conn = create_connection(self.db_file)
        with conn:
            cur = conn.cursor()
            cur.execute(table, query)
            items = cur.fetchone()
        return items[0]
    #
    #
    #
    # -----------------------------------------------
    #
    def function(self, beam_name: str | int,
                 load_name: str | int,
                 steps: int, Pa: float):
        """ """
        #
        comb2basic = combination_basic(mesh_name=self._mesh_id,
                                       db_file=self.db_file)
        #comb2basic = self.to_basic()
        header = ['load_name', 'mesh_name']
        comb_grp = comb2basic.groupby(header)
        #1 / 0
        # get basic load
        beam = self._beams[beam_name]
        #
        comb_item = comb_grp.get_group((load_name, self.mesh_name,))
        #
        #
        bfunction = beam.function(steps=steps,
                                  Pa=Pa, factor=1)
        # Axial   [FP, blank, blank, Fu]
        # torsion [T, B, Psi, Phi, Tw]
        # Bending [V, M, theta, w]
        #
        # [Fx, Fy, Fz, Mx, My, Mz]
        # [V, M, w, theta]
        header = ['load_name', 'mesh_name',
                  'load_comment', 'load_type',
                  'load_level', 'load_system',
                  'element_name', 'length',
                  'axial', 'torsion', 'VM_inplane', 'VM_outplane']
        #
        #          'FP', 'blank1', 'blank2', 'Fu',
        #          'T', 'B', 'Psi', 'Phi', 'Tw',
        #          'Vy', 'Mz', 'theta_y', 'w_y',
        #          'Vz', 'My', 'theta_z', 'w_z']
        df = DBframework()
        bfunction = df.DataFrame(data=bfunction, columns=header, index=None)
        grp_funct = bfunction.groupby(['load_name', 'element_name'])
        #
        loadfun = []
        for item in comb_item.itertuples():
            funct_item = grp_funct.get_group((item.basic_load, beam_name, )).reset_index()
            funct_item = funct_item.groupby(['load_name', 'element_name', 'length'])
            funct_item = funct_item[['axial', 'torsion', 'VM_inplane', 'VM_outplane']].sum()
            funct_item *= item.factor
            loadfun.append(funct_item.reset_index())
        #
        loadfun = df.concat(loadfun, ignore_index=True)
        loadfun = loadfun.groupby(['element_name', 'length'])
        loadfun = loadfun[['axial', 'torsion', 'VM_inplane', 'VM_outplane']].sum()
        loadfun.reset_index(inplace=True)        
        #print('-->')
        return loadfun      
    #
    # ----------------------------
    #
    def to_basicXX(self):
        """ """
        # get basic load
        conn = create_connection(self.db_file)
        with conn:
            # Load data
            query = (self._mesh_id, )
            table = "SELECT Load.*, Mesh.name \
                     FROM Load, Mesh \
                     WHERE Load.mesh_id = ? \
                     AND Mesh.number = Load.mesh_id; \
                    "
            cur = conn.cursor()
            cur.execute(table, query)
            load_data = cur.fetchall()
            #
            # combination data
            query = (self._mesh_id, )
            table = "SELECT LoadCombination.* \
                     FROM Load, LoadCombination \
                     WHERE Load.number = LoadCombination.load_id \
                     AND Load.mesh_id = ?; \
                     "
            cur.execute(table, query)
            combination = cur.fetchall()            
        #
        load_comb = {item[0]: item[1]
                     for item in load_data if item[3] == 'combination'}
        load_basic = {item[0]: item[1]
                      for item in load_data if item[3] == 'basic'}        
        # [number, ]
        load_df = {item[1]: [item[0], item[4], item[7]]
                   for item in load_data }
        #
        cbasic = defaultdict(list)
        ccomb = defaultdict(list)
        for item in combination:
            if item[2]: # basic load
                cbasic[item[1]].append([item[2], item[4]])
            elif item[3]: # combination
                ccomb[item[1]].append([item[3], item[4]])
        #
        comb1 = get_factor(load=cbasic)
        comb2 = get_factor(load=ccomb)   
        #
        dftemp = []
        test = defaultdict(list)
        for key, items in comb1.items():
            namec = load_comb[key]
            for precomb, prefactor in items.items():
                nameb =  load_basic[precomb]
                test[namec].append([nameb, prefactor])
        #
        for key, items in comb2.items():
            name = load_comb[key]
            for comb_item, factor in items.items():
                for precomb, prefactor in comb1[comb_item].items():
                    nameb =  load_basic[precomb]
                    test[name].append([nameb, prefactor * factor])
        #
        test = get_factor(load=test)
        #
        #number = 0
        for key, items in test.items():
            number = load_df[key][0]
            title = load_df[key][1]
            mesh_name = load_df[key][2]
            for item in items.items():
                dftemp.append([key, number, 'combination',
                               title, mesh_name, *item])

        #
        #
        header = ['load_name', 'load_id','load_level',
                  'load_title', 'mesh_name', 
                  'basic_load', 'factor']
        db = DBframework()
        df_comb = db.DataFrame(data=dftemp, columns=header, index=None)
        #1 / 0
        return df_comb
#
def combination_basic(mesh_name:str|int, db_file:str):
    """ load combination to basic"""
    conn = create_connection(db_file)
    with conn:
        # Load data
        query = (mesh_name,)
        table = "SELECT Load.*, Mesh.name \
                 FROM Load, Mesh \
                 WHERE Load.mesh_id = ? \
                 AND Mesh.number = Load.mesh_id; \
                "
        cur = conn.cursor()
        cur.execute(table, query)
        load_data = cur.fetchall()
        #
        # combination data
        query = (mesh_name,)
        table = "SELECT LoadCombination.* \
                 FROM Load, LoadCombination \
                 WHERE Load.number = LoadCombination.load_id \
                 AND Load.mesh_id = ?; \
                 "
        cur.execute(table, query)
        combination = cur.fetchall()
    #
    load_comb = {item[0]: item[1]
                 for item in load_data if item[3] == 'combination'}
    load_basic = {item[0]: item[1]
                  for item in load_data if item[3] == 'basic'}
    # [number, ]
    load_df = {item[1]: [item[0], item[4], item[7]]
               for item in load_data}
    #
    cbasic = defaultdict(list)
    ccomb = defaultdict(list)
    for item in combination:
        if item[2]:  # basic load
            cbasic[item[1]].append([item[2], item[4]])
        elif item[3]:  # combination
            ccomb[item[1]].append([item[3], item[4]])
    #
    comb1 = get_factor(load=cbasic)
    comb2 = get_factor(load=ccomb)
    #
    comb2basic = defaultdict(list)
    for key, items in comb1.items():
        comb_name = load_comb[key]
        for precomb, prefactor in items.items():
            basic_name = load_basic[precomb]
            comb2basic[comb_name].append([basic_name, prefactor])
    #
    for key, items in comb2.items():
        name = load_comb[key]
        for comb_item, factor in items.items():
            for precomb, prefactor in comb1[comb_item].items():
                nameb = load_basic[precomb]
                comb2basic[name].append([nameb, prefactor * factor])
    comb2basic = get_factor(load=comb2basic)
    #
    df_temp = []
    for key, items in comb2basic.items():
        number = load_df[key][0]
        title = load_df[key][1]
        mesh_name = load_df[key][2]
        for item in items.items():
            df_temp.append([key, number, 'combination',
                            title, mesh_name, *item])

    # --------------------------
    # to dataframe
    header = ['load_name', 'load_id', 'load_level',
              'load_title', 'mesh_name',
              'basic_load', 'factor']
    db = DBframework()
    df = db.DataFrame(data=df_temp, columns=header, index=None)
    return df
#
def get_factor(load):
    """ """
    comb = {}
    for key, items in load.items():
        factors = {}
        for item in items:
            try:
                factors[item[0]] += item[1]
            except KeyError:
                factors[item[0]] = item[1]
        #
        #comb[load_name[key]] = factors
        comb[key] = factors
    return comb
#
class BasicCombSQL(LoadCombinationRoot):
    """
    Basic Load Combination Class
    """
    __slots__ = ['_type', '_number', '_labels', 
                 'db_file', '_comb_name', '_mesh_id']

    def __init__(self, bl_type:str, comb_name:str|int,
                 mesh_id: int, db_file:str):
        """
        """
        super().__init__(name=comb_name)
        #self._name = comb_name
        self._labels: list[str,int] = []
        self._number: list[str,int] = []
        self._type = bl_type
        self.db_file = db_file
        self._mesh_id = mesh_id
    #
    #
    def __setitem__(self, load_name:int|str,
                    factor: float) -> None:
        """
        """
        self._labels.append(load_name)
        conn = create_connection(self.db_file)
        with conn:
            if self._type == 'basic':
                try:
                    lc_number = None
                    bl_number = get_load_number(conn, load_name,
                                                self._mesh_id,
                                                'basic')
                except TypeError:
                    raise IOError(f' Basic Load {load_name} not found')
                
            elif self._type == 'combination':
                try:
                    bl_number = None
                    lc_number = get_load_number(conn, load_name,
                                                self._mesh_id,
                                                'combination')
                except TypeError:
                    raise IOError(f' Load Combination {load_name} not found')
            
            else:
                raise IOError(f'load type {self._type} not valid')
            #
            load_id = get_load_number(conn, self._name,
                                      self._mesh_id, 'combination')
            #
            self._push_combination(conn, load_id, bl_number,
                                   lc_number, factor)
        #print('--?')
    #
    def __getitem__(self, load_name:int|str):
        """
        """
        try:
            comb_name = self._name
            db_file = self.db_file
            conn = create_connection(db_file)
            # get beam load
            with conn:
                comb_number = get_load_number(conn, comb_name,
                                              self._mesh_id, 'combination')
                if self._type == 'basic':
                    basicn = get_load_number(conn, load_name,
                                             self._mesh_id, 'basic')
                    
                    factor = get_load_factor(conn, comb_number, 
                                             'basic_id', basicn)
                
                elif self._type == 'combination':
                    combn = get_load_number(conn, load_name,
                                            self._mesh_id, 'combination')
                    
                    factor = get_load_factor(conn, comb_number,
                                             'combination_id', combn)
                
                else:
                    raise IOError(f'load type {self._type} not valid')
            return factor
        
        except TypeError:
            raise KeyError(f"{load_name}")

    #
    def _push_combination(self, conn, load_id:int, bl_name:int,
                          lc_name:int, factor:float):
        """
        """
        project = (load_id, bl_name, lc_name, factor)
        sql = 'INSERT INTO LoadCombination(\
               load_id, basic_id, combination_id, factor)\
               VALUES(?,?,?,?)'
        cur = conn.cursor()
        cur.execute(sql, project)
    #
    # ----------------------------
    #
    @property
    def load_type(self):
        """
        """
        return self._type
    #
    #
#   
#
def get_load_list(conn, mesh_id: int):
    """ """
    query = (mesh_id, )
    table = "SELECT * FROM Load WHERE mesh_id = ?"
    #
    cur = conn.cursor()
    cur.execute(table, query)
    loads = cur.fetchall()
    return loads
#
def get_load_number(conn, load_name:int|str,
                    mesh_id:int, level:str):
    """ """
    
    query = (load_name, mesh_id, level)
    table = "SELECT * FROM Load \
             WHERE name = ? \
             AND mesh_id = ? \
             AND level = ?"
    #
    cur = conn.cursor()
    cur.execute(table, query)
    loads = cur.fetchone()
    return loads[0]
#
def get_load_factor(conn, lc_id, load_type:str, 
                    load_number: int):
    """ """
    query = (lc_id, load_number,)
    table = f"SELECT * FROM LoadCombination\
              WHERE load_id = ? \
              AND {load_type} = ?"
    #
    cur = conn.cursor ()
    cur.execute(table, query)
    loads = cur.fetchall()
    factor = sum([item[4] for item in loads])
    return factor
#
def get_comb_factor(conn, load_id, lc_number):
    """ """
    query = (load_id, lc_number, )
    table = f"SELECT * FROM LoadCombination\
            WHERE load_id = ? \
            AND combination_id = ? "
    cur = conn.cursor ()
    cur.execute(table, query)
    loads = cur.fetchone()
    return loads[4]
#
def get_comb_number(conn, load_name:int,
                    mesh_id: int):
    """ """
    query = (load_name, mesh_id, )
    table = f"SELECT * FROM Load\
            WHERE name = ? \
            AND type = 'combination' \
            AND mesh_id = ?"
    cur = conn.cursor()
    cur.execute(table, query)
    loads = cur.fetchone()
    return loads[0]