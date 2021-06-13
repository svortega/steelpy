# 
# Copyright (c) 2019-2021 steelpy
#
# 
# Python stdlib imports
from array import array
from collections.abc import Mapping
#from collections import defaultdict
from collections import Counter, defaultdict
#from dataclasses import dataclass
from typing import Union, Dict, List, Union
#from math import prod

# package imports
from steelpy.f2uModel.load.operations.operations import duplicates, indices
from steelpy.f2uModel.load.operations.combination import LoadCombinationBasic
from steelpy.f2uModel.results.sqlite.operation.process_sql import create_connection

#
#
class LoadCombSQL(LoadCombinationBasic):
    
    __slots__ = ['_labels', '_title', '_number', 'bd_file',
                 '_index', '_basic', '_combination', '_metocean']
    #
    def __init__(self, bd_file:str):
        """
        """
        super().__init__()
        self.bd_file = bd_file
        #
        self._basic = BasicLoadSQL(self)
        self._combination = LoadCombinationSQL(self)
    #
    def __setitem__(self, load_name:int, load_title:str) -> None:
        """
        """
        try:
            self._labels.index(load_name)
            raise Exception('    *** warning load combination name {:} already exist'
                            .format(load_name))
        except ValueError:
            self._labels.append(load_name)
            self._title.append(load_title)
            conn = create_connection(self.bd_file)
            with conn:
                load_number = self._push_load_combination(conn, load_name, load_title)
                self._number.append(load_number)
                conn.commit()
    #
    def __getitem__(self, load_name:Union[str,int]):
        """
        """
        try:
            self._index = self._labels.index(load_name)
            return CombTypeSQL(self)
        except ValueError:
            raise IOError("load combination {:} not defined".format(load_name))
    #
    def _push_load_combination(self, conn, load_name:int, load_title:str):
        """ """
        #
        project = (load_name, load_title, "combination")
        sql = 'INSERT INTO tb_Load(name, title, type) VALUES(?,?,?)'
        cur = conn.cursor()
        cur.execute(sql, project)
        return cur.lastrowid
    #
    def to_basic(self):
        """ """
        conn = create_connection(self.bd_file)
        cur = conn.cursor()
        cur.execute ( "SELECT load_number, bl_number, factor\
                      FROM tb_LoadCombIndex\
                      WHERE bl_number IS NOT NULL")
        basic_loads = cur.fetchall()
        # get basic load
        blc = defaultdict(list)
        for basic in basic_loads:
            blc[basic[0]].append([basic[1], basic[2]])
        # get load combination and convert to basic loads
        cur.execute ( "SELECT load_number, lc_number, factor\
                      FROM tb_LoadCombIndex\
                      WHERE lc_number IS NOT NULL")
        comb_loads = cur.fetchall()
        for comb in comb_loads:
            for bl in blc[comb[1]]:
                blc[comb[0]].append([bl[0], bl[1]*comb[2]])
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
                    basic_loads[ lname ].append ( [ blname, sum ( factor ) ] )
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
        return basic_loads
    #
    def solve_combinations(self, basic_res, memb_force):
        """
        """
        comb_res = {}
        memb_comb = {}
        bloads = self.to_basic()
        for lcomb, comb in bloads.items():
            #lcomb = load_combination[cname].title
            beam_load = {}
            for bname, factor in comb:
                try:
                    comb_res[lcomb] += basic_res[bname] * factor
                except KeyError:
                    comb_res[lcomb] = basic_res[bname] * factor
                #
                for mname, member in memb_force[bname].items():
                    try:
                        beam_load[mname] += member * factor
                    except KeyError:
                        beam_load[mname] =  member * factor
            memb_comb[lcomb] = beam_load
        return comb_res, memb_comb
#     
#
class CombTypeSQL:
    """
    """
    __slots__ = ['_cls']

    def __init__(self, cls):
        """
        """
        self._cls = cls

    #
    @property
    def load_combination(self):
        """
        """
        return self._cls._combination

    #
    @property
    def basic_load(self):
        """
        """
        return self._cls._basic
#
class BasicLoadSQL(Mapping):
    """
    FE Metocean Combination Class
    """
    __slots__ = ['_cls', '_labels', '_number']

    def __init__(self, cls):
        """
        """
        self._cls = cls
        self._labels: List[Union[str,int]] = []
        self._number: List[Union[str,int]] = []

    #
    def __getitem__(self, basic_load_name:int):
        """
        """
        #1/0
        try:
            index = self._cls._index
            comb_name = self._cls._labels[index]
            #load_title = self._title[index]
            # load_number = self._load_number.index(load_name)
            bd_file = self._cls.bd_file
            conn = create_connection(bd_file)
            comb_number = get_comb_number( conn, comb_name )
            basic_number = self._get_basic_load_number(conn, basic_load_name)
            # get beam load
            with conn:
                factor = self._get_basic_load_factor(conn, comb_number, basic_number)
            return factor
        except TypeError:
            raise KeyError(f"{basic_load_name}")

    def __setitem__(self, basic_load_name:int, factor: float) -> None:
        """
        """
        index = self._cls._index
        #comb_title = self._cls._title[index]
        comb_name = self._cls._labels[index]
        conn = create_connection(self._cls.bd_file)
        load_number = get_comb_number(conn, comb_name)
        self._labels.append(basic_load_name)
        self._number.append(comb_name)
        with conn:
            self._push_combination(conn, load_number, basic_load_name, factor)
            conn.commit()
    #
    def _push_combination(self, conn, load_number:int, bl_name:int,
                          factor:float):
        """
        """
        #
        #try:
        basic_number = self._get_basic_load_number(conn, bl_name)
        #
        project = (load_number, basic_number, None, factor)
        sql = 'INSERT INTO tb_LoadCombIndex(\
               load_number, bl_number, lc_number, factor)\
               VALUES(?,?,?,?)'
        cur = conn.cursor()
        cur.execute(sql, project)
    #
    def _get_basic_load_number(self, conn, load_name:int):
        """ """
        cur = conn.cursor()
        cur.execute("SELECT * FROM tb_Load\
                     WHERE name = {:} \
                     AND type = 'basic'".format(load_name))
        loads = cur.fetchone()
        return loads[0]
    #
    def _get_basic_load_factor(self, conn, load_number, basic_number):
        """ """
        cur = conn.cursor ()
        cur.execute ( "SELECT * FROM tb_LoadCombIndex\
                     WHERE load_number = {:} \
                     AND bl_number = {:}".format(load_number, basic_number))
        loads = cur.fetchone ()
        return loads[4]
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
        index = self._cls._index
        comb_name = self._cls._labels[index]
        basic_load = [self._labels[x] for x, item in enumerate(self._number)
                      if item == comb_name]
        #number = self._number
        #basic = set(self._labels)
        return iter(basic_load)
#
#
class LoadCombinationSQL(Mapping):
    
    __slots__ = ['_cls', '_labels', '_number']

    def __init__(self, cls):
        """
        """
        self._cls = cls
        self._labels: List[Union[str, int]] = []
        self._number: List[Union[str, int]] = []
    #
    def __setitem__(self, load_name:int, factor: float) -> None:
        """
        """
        index = self._cls._index
        #comb_title = self._cls._title[index]
        comb_name = self._cls._labels[index]
        conn = create_connection(self._cls.bd_file)
        load_number = get_comb_number(conn, comb_name)
        with conn:
            self._push_combination(conn, load_number,
                                   load_name, factor)
            #conn.commit()
        #
        self._labels.append(load_name)
        self._number.append(comb_name)
    #
    def __getitem__(self, load_name:int):
        """
        """
        try:
            index = self._cls._index
            comb_name = self._cls._labels[index]
            bd_file = self._cls.bd_file
            conn = create_connection(bd_file)
            comb_number = get_comb_number(conn, comb_name)
            basic_number = get_comb_number(conn, load_name)
            # get beam load
            with conn:
                factor = self._get_comb_load_factor(conn, comb_number, basic_number)
            return factor
        except TypeError:
            raise KeyError(f"{load_name}")
    #
    def _push_combination(self, conn, load_number:int,
                            lc_number:int, factor:float):
        """
        """
        try:
            comb_number = get_comb_number(conn, lc_number)
        except TypeError:
            raise IOError(f"combination {lc_number} not found")
        #
        project = (load_number,  None, comb_number, factor)
        sql = 'INSERT INTO tb_LoadCombIndex(\
               load_number, bl_number, lc_number, factor)\
               VALUES(?,?,?,?)'
        cur = conn.cursor()
        cur.execute(sql, project)    
    #
    def _get_comb_load_factor(self, conn, load_number, basic_number):
        """ """
        cur = conn.cursor()
        cur.execute("SELECT * FROM tb_LoadCombIndex\
                     WHERE load_number = {:} \
                     AND lc_number = {:}".format(load_number, basic_number))
        loads = cur.fetchone()
        return loads[4]
    #
    def __len__(self) -> int:
        return len(self._labels)
    #
    def __contains__(self, value) -> bool:
        return value in self._labels
    #
    def __iter__(self):
        """
        """
        index = self._cls._index
        comb_name = self._cls._labels[index]
        comb_load = [self._labels[x] for x, item in enumerate(self._number)
                      if item == comb_name]
        #comb = set(self._labels)
        return iter(comb_load)
#   
#
def get_comb_number(conn, load_name:int):
    """ """
    cur = conn.cursor()
    cur.execute("SELECT * FROM tb_Load\
                 WHERE name = {:} \
                 AND type = 'combination'".format(load_name))
    loads = cur.fetchone()
    return loads[0]
#
#
def get_load_list(conn):
    """ """
    cur = conn.cursor()
    cur.execute("SELECT * FROM tb_Load")
    loads = cur.fetchall()
    return loads
#
#