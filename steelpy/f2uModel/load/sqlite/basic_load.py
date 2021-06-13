#
# Copyright (c) 2009-2020 fem2ufo
# 

# Python stdlib imports
from array import array
from collections.abc import Mapping
#from dataclasses import dataclass
from typing import NamedTuple, Tuple, List, Iterator, Dict, Iterable, ClassVar, Union
#import re

# package imports
from steelpy.f2uModel.load.operations.actions import SelfWeight
from steelpy.f2uModel.load.sqlite.element import BeamDistributedSQL, BeamPointSQL
from steelpy.f2uModel.load.sqlite.node import  NodeLoadSQL
from steelpy.f2uModel.results.sqlite.operation.process_sql import create_connection
from steelpy.f2uModel.load.operations.basic_load import BasicLoadBasic
#
#
# ---------------------------------
#
#
class BasicLoadSQL(BasicLoadBasic):
    """
    FE Load Cases

    LoadType
        |_ name
        |_ number
        |_ basic
        |_ combination_level
        |_ time_series
        |_
        |_ temperature

    **Parameters**:
      :number:  integer internal number
      :name:  string node external name
    """
    __slots__ = ['bd_file', '_labels', '_number', 
                  '_title', 'gravity', # '_index',
                 '_nodal_load', '_beam_line',
                 '_beam_point', '_selfweight']

    #
    def __init__(self, bd_file:str):
        """
        """
        super().__init__()
        self.bd_file = bd_file
        # create node table
        #self._create_table()
        #
        #self._load = LoadTypeSQL(self)
        self._nodal_load = NodeLoadSQL(self)
        self._beam_line = BeamDistributedSQL(self)
        self._beam_point = BeamPointSQL(self)
        self._selfweight = SelfWeight()
    #
    #
    def __setitem__(self, load_name:int, load_title:str) -> None:
        """
        """
        try:
            self._labels.index(load_name)
            raise Warning('    *** warning load name {:} already exist'
                            .format(load_name))
        except ValueError:
            self._labels.append(load_name)
            self._title.append(load_title)
            conn = create_connection(self.bd_file)
            with conn:
                load_number = self._push_basic_load(conn, load_name, load_title)
                self._number.append(load_number)
                #conn.commit()
    #
    def __getitem__(self, load_name:Union[str,int]):
        """
        """
        try:
            self._index = self._labels.index(load_name)
            #print(load_name, self._index)
            return LoadTypeSQL(self)
        except ValueError:
            raise IOError("load case {:} no defined".format(load_name))
    #
    def _push_basic_load(self, conn, load_name:int, load_title:str):
        """ """
        project = (load_name, load_title, "basic")
        sql = 'INSERT INTO tb_Load(name, title, type) VALUES(?,?,?)'
        cur = conn.cursor()
        cur.execute(sql, project)
        return cur.lastrowid
#
#
class LoadTypeSQL:
    """
    """
    __slots__ = ['_cls']

    def __init__(self, cls):
        """
        """
        self._cls = cls
        #self._nodal_load = NodeLoadSQL(cls)
        #self._beam_line = BeamDistributedSQL(cls)
        #self._beam_point = BeamPointSQL(cls)
        #self._selfweight = SelfWeight()
    #
    @property
    def name(self):
        """ """
        return self._cls._labels[self._cls._index]

    #@name.setter
    #def name(self, name:str):
    #    """ """
    #    self._cls._title[self._cls._index] = name
    #
    @property
    def number(self):
        """ """
        return self._cls._number[self._cls._index]

    #@number.setter
    #def number(self, name:str):
    #    """ """
    #    self._cls._number[self._cls._index] = name
    #
    #
    @property
    def title(self):
        """ """
        return self._cls._title[self._cls._index]

    #@title.setter
    #def title(self, title:str):
    #    """ """
    #    self._cls._title[self._cls._index] = title
    #
    @property
    def selfweight(self):
        """
        The self weight form allows you to specify multipliers to
        acceleration due to gravity (g) in the X, Y, and Z axes.
        If switched on, the default self weight acts in the Y axis
        with a magnitude and sign of -1."""
        return self._cls._selfweight
    #
    @property
    def point_node(self):
        """
        """
        return self._cls._nodal_load

    @point_node.setter
    def point_node(self, values: List):
        """
        Point Load
        """
        for value in values:
            self._cls._nodal_load[value[0]] = value[1:]
    #
    #
    @property
    def line_beam(self):
        """
        Linear Varying Load (lvl) - Non Uniformly Distributed Load
        """
        return self._cls._beam_line
    
    @line_beam.setter
    def line_beam(self, values:List):
        """
        Linear Varying Load (lvl) - Non Uniformly Distributed Load
                value : [qx1, qy1, qz1, qx2, qy2, qz2, L1, L2]
        
                        |
             q1         | q2
        o------|        |----------o
        |                          |
        +  L1  +        +    L2    +
        """
        for value in values:
            self._cls._beam_line[value[0]] = value[1:]
    #
    #
    #@property
    #def udl_beam(self):
    #    """
    #    Uniformly Distributed Load (udl)
    #    """
    #    return self._cls._beam_line
    #
    #@udl_beam.setter
    #def udl_beam(self, values:List):
    #    """
    #    Uniformly Distributed Load (udl)
    #    value : [qx1, qy1, qz1, qx2, qy2, qz2, L1, L2]
    #
    #                    |
    #         q1         | q2
    #    o------|        |----------o
    #    |                          |
    #    +  L1  +        +    L2    +
    #    """
    #    for value in values:
    #        self._cls._beam_line[value[0]] = value[1:]
    #
    #
    @property
    def point_beam(self):
        """ Concentrated force """
        return self._cls._beam_point

    @point_beam.setter
    def point_beam(self, values: List):
        """
        Concentrated force
        """
        for value in values:
            self._cls._beam_point[value[0]] = value[1:]
    #
    @property
    def node(self):
        """
        """
        return self._cls._nodal_load
    #
    #@node.setter
    #def node(self, values:List):
    #    """
    #    Node Load = [node_number, 'point', x,y,z,mx,my,mz],
    #                [node_number, 'mass' , x,y,z]
    #    """
    #    if isinstance(values[1], str):
    #        self._nodal_load[values[0]] = values[2:]
    #    else:
    #        for value in values:
    #            self._nodal_load[value[0]] = value[2:]    
    #    
#
#
#