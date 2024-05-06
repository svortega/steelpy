# 
# Copyright (c) 2019 steelpy
#
# 
# Python stdlib imports
from __future__ import annotations
from collections.abc import Mapping
#from collections import defaultdict
from collections import defaultdict
#from dataclasses import dataclass
from typing import NamedTuple
#from math import prod

# package imports
from steelpy.ufo.load.sqlite.load_case import BasicLoadSQL
from steelpy.ufo.load.process.operations import duplicates, indices
from steelpy.ufo.load.process.combination import LoadCombinationBasic
# steelpy
from steelpy.utils.sqlite.utils import create_connection, create_table
from steelpy.utils.dataframe.main import DBframework
#
#
class LoadCombSQL(LoadCombinationBasic):
    
    __slots__ = ['_labels', '_title', '_number', 'db_file', '_plane', 
                 '_index', '_basic', '_combination', '_metocean',
                 '_component']
    #
    def __init__(self, db_file:str, plane: NamedTuple,
                 component: int):
        """
        """
        super().__init__()
        #
        self.db_file = db_file
        self._plane=plane
        self._component = component
        #
        self._basic = BasicLoadSQL(db_file=self.db_file,
                                   component=component,
                                   plane=self._plane)
        self._combination = {}
        #
        conn = create_connection(self.db_file)
        with conn: 
            self._create_table(conn)
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
                                                       component=self._component, 
                                                       db_file=self.db_file)
            #
            conn = create_connection(self.db_file)
            with conn:
                load_no = self._push_combination(conn, load_name,
                                                 self._component, 
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
            raise IOError("load combination {:} not defined".format(load_name))
    #
    def _create_table(self, conn):
        """ """
        table = "CREATE TABLE IF NOT EXISTS LoadCombination(\
                number INTEGER PRIMARY KEY NOT NULL,\
                load_id INTEGER NOT NULL REFERENCES Load(number),\
                bl_number INTEGER REFERENCES Load(number),\
                lc_number INTEGER REFERENCES Load(number),\
                factor DECIMAL NOT NULL);"
        create_table(conn, table)
    #
    def _push_combination(self, conn, name:int|str,
                          component:int, level:str, title:str):
        """
        """
        project = (name, component, level, title)
        sql = 'INSERT INTO Load(\
               name, component_id, level, title)\
               VALUES(?,?,?,?)'
        #
        cur = conn.cursor()
        cur.execute(sql, project)
        #conn.commit()
        return cur.lastrowid
        
    #
    def to_basic2(self):
        """ """
        #
        1/0
        # get basic load
        conn = create_connection(self.db_file)
        cur = conn.cursor()
        cur.execute("SELECT load_id, bl_number, factor\
                      FROM LoadCombination\
                      WHERE bl_number IS NOT NULL")
        basic_loads = cur.fetchall()
        # get basic load factors
        blc = defaultdict(list)
        for basic in basic_loads:
            blc[basic[0]].append([basic[1], basic[2]])
        #
        # get load combination and convert to basic loads
        cur.execute("SELECT load_id, lc_number, factor\
                      FROM LoadCombination\
                      WHERE lc_number IS NOT NULL")
        comb_loads = cur.fetchall()
        for comb in comb_loads:
            for bl in blc[comb[1]]:
                blc[comb[0]].append([bl[0], bl[1]*comb[2]])
        #
        # organize basic load in user load name
        load_list = get_load_list(conn)
        lnumber = {item[0]:item[2] for item in load_list}
        basic_loads = defaultdict(list)
        for key, item in blc.items():
            lname = lnumber[key]
            tlist = list(zip(*item))
            dup = duplicates(tlist[0])
            dup = indices(tlist[0], dup)
            if dup:
                # duplicates
                for name, load in dup.items():
                    blname = lnumber[name]
                    factor = [tlist[1][index] for index in load]
                    #basic_loads[key].append([name, sum(factor)])
                    basic_loads[ lname ].append( [ blname, sum ( factor ) ] )
                    #dftemp.append([key, item.number, 'combination', item.title, blname, sum(factor)])
                # singles
                sload = set(tlist[0]) - set(dup)
                for name in sload:
                    blname = lnumber[name]
                    index = tlist[0].index(name)
                    #basic_loads[key].append([name, tlist[1][index]])
                    basic_loads[lname].append([blname, tlist[1][index]])
            else:
                #basic_loads[key] = item
                basic_loads[lname] = [[lnumber[_bl[0]], _bl[1]] for _bl in item]
                #dftemp.extend([[key, item.number, 'combination', item.title, lnumber[_bl[0]], _bl[1]]
                #               for _bl in item])
        #
        #header = ['load_name', 'load_id','load_type', 'load_title', 'basic_load', 'factor']
        #dfcomb = db.DataFrame(data=dftemp, columns=header, index=None)
        #return dfcomb #, basic_loads
        return basic_loads
    #
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
    #
    def function(self, steps:int,
                 Fb_local):
        """ """
        #print('load function')
        basic = self._basic
        dfcomb = self.to_basic()
        header = ['load_name', 'component_name']
        combgrp = dfcomb.groupby(header)
        Fbgrp = Fb_local.groupby(header)
        #
        loadfun = []
        for key, items in combgrp:
            Fb = Fbgrp.get_group(key)
            Fb.set_index('element_name', inplace=True)
            for item in items.itertuples():
                #print(item)
                lbasic = basic[item.basic_load]
                lfout = lbasic._beam.load_function(Fb=Fb,
                                                   factor=item.factor,
                                                   steps=steps)
                loadfun.extend(lfout)
        #
        header = ['load_name', 'basic_load',
                  'component_name',
                  'load_comment', 'load_type', 
                  'load_level', 'load_system',
                  'element_name', 'node_end']
                  #'axial', 'torsion', 'VM_inplane', 'VM_outplane']
        values = ['axial', 'torsion', 'VM_inplane', 'VM_outplane']
        #
        df = DBframework()
        dfload = df.DataFrame(data=loadfun,
                              columns=header + values,
                              index=None)       
        #
        #
        header = ['load_name', 'component_name',
                  #'load_comment', 'load_type',
                  'load_level', 'load_system',
                  'element_name', 'node_end']
        dfload.drop(columns=['basic_load', 'load_comment', 'load_type'],
                    inplace=True)
        dfload = dfload.groupby(header, as_index=False)[values].sum()
        dfload['load_level'] = 'combination'
        return dfload
    #
    def Fn(self):
        """ """
        dfbasic = self._basic.Fn()
        dfcomb = self.to_basic()
        #
        values = ['Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz',
                  'x', 'y', 'z', 'rx', 'ry', 'rz']
        try:
            dfnew = self.update_combination(dfbasic, dfcomb,
                                            values=values)
        except ValueError:
            raise IOError('2nd order requires Load Combination as input')
        #
        colgrp = ['load_name', 'load_id', 
                  'load_level', 'load_title',
                  'load_system', 'component_name',
                  'node_name', 'node_index']
        #        
        dfnew = dfnew.groupby(colgrp, as_index=False)[values].sum()
        return dfnew
    #
    def ENL(self):
        """Equivalent Nodal Loads """
        #print('comb enl')
        ENLbasic = self._basic.ENL()
        cols = tuple(set(ENLbasic['load_name'].tolist()))
        if not cols:
            return ENLbasic
        #
        dfcomb = self.to_basic()
        grpcomb = dfcomb.groupby(['basic_load'])
        #
        # basic loading
        values = ['Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz',
                  'Psi', 'B', 'Tw']
        dfnew = self.update_combination(ENLbasic,
                                        grpcomb.get_group(cols),
                                        values=values)
        #
        colgrp = ['load_name', 'load_id', 
                  'load_level', 'load_title',
                  'load_system', 'component_name',
                  'element_name', 
                  'node_name', 'node_index']        
        dfnew = dfnew.groupby(colgrp, as_index=False)[values].sum()
        #
        return dfnew
    #
    def update_combination(self, dfbasic, dfcomb,
                            values:list[str]):
        """
        Update node displacements to include lcomb
        """
        db = DBframework()
        # group basic load by name
        header = ['component_name','load_name']
        blgrp = dfbasic.groupby(header)
        combgrp = dfcomb.groupby(header)
        dftemp = []
        for key, combfactors in combgrp:
            for row in combfactors.itertuples():
                try: # check for beam load only
                    comb = blgrp.get_group((key[0], row.basic_load,)).copy()
                except KeyError:
                    continue
                comb.loc[:, values] *= row.factor
                comb['load_level'] = 'combination'
                comb['load_name'] = row.load_name
                comb['load_id'] = row.load_id
                comb['load_title'] = row.load_title
                dftemp.append(comb)
        #
        dftemp = db.concat(dftemp, ignore_index=True)
        return dftemp
    #
#     
#
class CombTypeSQL:
    """
    """
    __slots__ = ['name', 'number', 'title', 
                 '_basic', '_metocean', '_combination',
                 'db_file', '_component']

    def __init__(self, name:str, title:str,
                 component: int, db_file:str):
        """
        """
        self._component = component
        self.name = name
        self.title = title
        self.db_file = db_file
        #
        self._basic = BasicCombSQL(bl_type="basic",
                                   comb_name=self.name,
                                   component=self._component, 
                                   db_file=self.db_file)
        
        self._metocean = BasicCombSQL(bl_type="metocean",
                                      comb_name=self.name,
                                      component=self._component, 
                                      db_file=self.db_file)
        
        self._combination = BasicCombSQL(bl_type="combination",
                                         comb_name=self.name,
                                         component=self._component, 
                                         db_file=self.db_file)        
        #
        #conn = create_connection(self.db_file)
        #with conn: 
        #    self._create_table(conn)
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
    def _component_name(self):
        """ component name """
        query = (self._component, )
        table = 'SELECT name \
                 FROM Component WHERE number = ?;'
        conn = create_connection(self.db_file)
        with conn:
            cur = conn.cursor()
            cur.execute(table, query)
            items = cur.fetchone()
        return items[0]
#
#
class BasicCombSQL(Mapping):
    """
    FE Metocean Combination Class
    """
    __slots__ = ['_labels', '_type', '_number',
                 'db_file', '_comb_name', '_component']

    def __init__(self, bl_type:str, comb_name:str|int,
                 component: int, db_file:str):
        """
        """
        self._comb_name = comb_name
        self._labels: list[str,int] = []
        self._number: list[str,int] = []
        self._type = bl_type
        self.db_file = db_file
        self._component = component
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
                                                self._component,
                                                'basic')
                except TypeError:
                    raise IOError(f' Basic Load {load_name} not found')
                
            elif self._type == 'combination':
                try:
                    bl_number = None
                    lc_number = get_load_number(conn, load_name,
                                                self._component,
                                                'combination')
                except TypeError:
                    raise IOError(f' Load Combination {load_name} not found')
            
            else:
                raise IOError(f'load type {self._type} not valid')
            #
            load_id = get_load_number(conn, self._comb_name,
                                      self._component, 'combination')
            #
            self._push_combination(conn, load_id, bl_number,
                                   lc_number, factor)
            
    #
    def __getitem__(self, load_name:int|str):
        """
        """
        try:
            comb_name = self._comb_name
            db_file = self.db_file
            conn = create_connection(db_file)
            # get beam load
            with conn:
                comb_number = get_load_number(conn, comb_name,
                                              self._component, 'combination')
                if self._type == 'basic':
                    basicn = get_load_number(conn, load_name,
                                             self._component, 'basic')
                    
                    factor = get_load_factor(conn, comb_number, 
                                             'bl_number', basicn)
                
                elif self._type == 'combination':
                    combn = get_load_number(conn, load_name,
                                            self._component, 'combination')
                    
                    factor = get_load_factor(conn, comb_number,
                                             'lc_number', combn)
                
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
               load_id, bl_number, lc_number, factor)\
               VALUES(?,?,?,?)'
        cur = conn.cursor()
        cur.execute(sql, project)
    #
    @property
    def load_type(self):
        """
        """
        return self._type
    #
    def __len__(self) -> int:
        return len(self._labels)

    def __contains__(self, value) -> bool:
        return value in self._labels
    #
    def __iter__(self):
        """
        """
        return iter(self._labels)
#   
#
def get_load_list(conn, component: int):
    """ """
    query = (component, )
    table = "SELECT * FROM Load WHERE component_id = ?"
    #
    cur = conn.cursor()
    cur.execute(table, query)
    loads = cur.fetchall()
    return loads
#
def get_load_number(conn, load_name:int|str,
                    component:int, level:str):
    """ """
    
    query = (load_name, component, level)
    table = "SELECT * FROM Load \
             WHERE name = ? \
             AND component_id = ? \
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
    loads = cur.fetchone ()
    return loads[4]
#
def get_comb_factor(conn, load_id, lc_number):
    """ """
    query = (load_id, lc_number, )
    table = f"SELECT * FROM LoadCombination\
            WHERE load_id = ? \
            AND lc_number = ? "
    cur = conn.cursor ()
    cur.execute(table, query)
    loads = cur.fetchone()
    return loads[4]
#
def get_comb_number(conn, load_name:int,
                    component: int):
    """ """
    query = (load_name, component, )
    table = f"SELECT * FROM Load\
            WHERE name = ? \
            AND type = 'combination' \
            AND component_id = ?"
    cur = conn.cursor()
    cur.execute(table, query)
    loads = cur.fetchone()
    return loads[0]