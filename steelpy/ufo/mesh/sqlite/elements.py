# 
# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
#from dataclasses import dataclass
#from array import array
from collections import Counter #, defaultdict
from collections.abc import Mapping
#from math import dist
#from typing import NamedTuple
from itertools import chain
import re

#
# package imports
from steelpy.ufo.mesh.sqlite.beam import BeamSQL
from steelpy.ufo.mesh.sqlite.utils import get_connectivity #, get_elements
#from steelpy.utils.sqlite.main import ClassBasicSQL
from steelpy.utils.sqlite.utils import create_connection, create_table
from steelpy.ufo.utils.element import get_element #, find_element_type
from steelpy.ufo.utils.element import ElementMain
#from steelpy.utils.dataframe.main import DBframework
#
#
#
class ElementsSQL(ElementMain):
    __slots__ = ['_beams', '_plane', 'db_file', '_component']
    
    def __init__(self, db_file:str,
                 #plane: NamedTuple,
                 component: int, 
                 db_system:str="sqlite") -> None:
        """
        """
        super().__init__(component)
        #self._component = component
        self._beams = BeamSQL(db_file=db_file,
                              component=component)
        #
        self.db_file = db_file
        # create table
        conn = create_connection(self.db_file)
        with conn:
            self._new_table(conn)
    #
    #
    @property
    def _labels(self):
        """ """
        query = (self._component, )
        table = 'SELECT Element.name FROM Element \
                 WHERE mesh_id = ? '
        conn = create_connection(self.db_file)
        with conn:        
            cur = conn.cursor()
            cur.execute(table, query)
            items = cur.fetchall()
        return [item[0] for item in items]
    #
    @property
    def _type(self):
        """ """
        query = (self._component, )
        table = 'SELECT Element.type FROM Element \
                 WHERE mesh_id = ? '
        conn = create_connection(self.db_file)
        with conn:        
            cur = conn.cursor()
            cur.execute(table, query)
            items = cur.fetchall()
        return [item[0] for item in items]
    #
    #
    # ---------------------------------
    #
    def __setitem__(self, element_name: int|str,
                    parameters: list|tuple|dict) -> None:
        """
        element_name
        parameters = [element_type, nodei..n,
                      material, section, roll_angle,
                      title, idx]
        """
        try:
            self._labels.index(element_name)
            raise Exception('element {:} already exist'.format(element_name))
        except ValueError:
            #
            parameters = get_element(parameters)
            element_type = parameters[0]
            #
            concept_id = parameters[6:]
            #
            if re.match(r"\b(shell(s)?|plate(s)?)\b", element_type, re.IGNORECASE):
                #return self._curve[mat_number]
                raise NotImplementedError('shell element tobe implemented')
            elif re.match(r"\b(beam)\b", element_type, re.IGNORECASE):
                self._beams[element_name] = parameters[1:6]
                number =  self._beams[element_name].number
            else:
                raise IOError(f' element type {element_type} not recognised')
            #
            # Update concept
            if concept_id:
                conn = create_connection(self.db_file)
                with conn:              
                    self._update_element(conn, number=number,
                                         title=concept_id[0],
                                         idx=concept_id[1])
    #
    def __getitem__(self, element_name: int|str):
        """
        """
        try:
            index = self._labels.index(element_name)
            element_type = self._type[index]
        except ValueError:
            raise KeyError(f' element name: {element_name} not valid')        
        #
        if re.match(r"\b(shell(s)?|plate(s)?)\b", element_type, re.IGNORECASE):
            raise NotImplementedError('shell element tobe implemented')
        elif re.match(r"\b(beam)\b", element_type, re.IGNORECASE):
            return self._beams[element_name]
        else:
            raise IOError(f' element type {element_type} not valid')    
    #
    #
    # ---------------------------------
    # SQL ops      
    #
    #
    def _new_table(self, conn) -> None:
        """ """
        elements = "CREATE TABLE IF NOT EXISTS Element(\
                            number INTEGER PRIMARY KEY NOT NULL,\
                            name NOT NULL,\
                            mesh_id INTEGER NOT NULL REFERENCES Mesh(number), \
                            type TEXT NOT NULL,\
                            material_id INTEGER NOT NULL REFERENCES Material(number),\
                            section_id INTEGER NOT NULL REFERENCES Section(number),\
                            roll_angle DECIMAL, \
                            title TEXT, \
                            concept_idx INTEGER);"
        #
        connectivity = "CREATE TABLE IF NOT EXISTS ElementConnectivity(\
                                number INTEGER PRIMARY KEY NOT NULL,\
                                element_id INTEGER NOT NULL REFERENCES Element(number),\
                                node_id INTEGER REFERENCES Node(number),\
                                node_end INTEGER NOT NULL);"
        #
        univectors = "CREATE TABLE IF NOT EXISTS ElementDirectionCosine(\
                                number INTEGER PRIMARY KEY NOT NULL,\
                                element_id INTEGER NOT NULL REFERENCES Element(number),\
                                x DECIMAL,\
                                y DECIMAL,\
                                z DECIMAL,\
                                axis INTEGER);"
        #
        offset = "CREATE TABLE IF NOT EXISTS ElementEccentricity(\
                            number INTEGER PRIMARY KEY NOT NULL,\
                            element_id INTEGER NOT NULL REFERENCES Element(number),\
                            node_id INTEGER REFERENCES Node(number),\
                            node_end INTEGER NOT NULL,\
                            system TEXT NOT NULL,\
                            x DECIMAL,\
                            y DECIMAL,\
                            z DECIMAL);"
        #
        #conn = create_connection(self.db_file)
        create_table(conn, elements)
        create_table(conn, connectivity)
        create_table(conn, offset)
        create_table(conn, univectors)
    #
    #
    def _update_element(self, conn, number: int,
                        title: str, idx: int):
        """ title, idx """
        #item = 'concept'
        #
        query = (title, idx, number,)
        table = f"UPDATE Element \
                 SET title = ?, \
                     concept_idx = ? \
                 WHERE number = ?;"
        #
        conn = create_connection(self.db_file)
        with conn:          
            cur = conn.cursor()
            comp = cur.execute(table, query)
        #
        if not comp:
            raise IOError(f' element {number} not valid')
    #
    #
    @property
    def get_connectivities(self):
        """ """
        query = (self._component)
        table = f"SELECT name FROM Element \
                WHERE mesh_id = ?;"
        #
        conn = create_connection(self.db_file)
        with conn:
            cur = conn.cursor()
            cur.execute(table, query)
            elements = cur.fetchall()
        #
        connodes = []
        for element in elements:
            connodes.append(get_connectivity(conn, element[0]))
        conn.close()
        return connodes
    #    
    #
    # ---------------------------------
    #
    #def iter_elements(self, arraysize=1000):
    #    """
    #    """
    #    conn = create_connection(self.db_file)
    #    cur = conn.cursor()
    #    # TODO: check if direction cosines given
    #    cur.execute("SELECT Element.name, Element.number, Element.type,\
    #                Element.roll_angle, Element.material, Element.section\
    #                FROM Element;" )
    #    #
    #    try:
    #        while True:
    #            elements = cur.fetchmany(arraysize)
    #            if not elements:
    #                break
    #            for element in elements:
    #                #cur.execute("SELECT ElementConnectivity.node_end, ElementConnectivity.node_name\
    #                #            FROM ElementConnectivity\
    #                #            WHERE ElementConnectivity.element_name = {:};".format(element[0]))
    #                #row = cur.fetchall()
    #                #connodes = [x for _, x in sorted(row)]
    #                connodes = get_connectivity(conn, element[0])
    #                data = [*element[0:6], connodes, self.db_file]
    #                yield BeamElement(data)
    #    except Exception as e:
    #        print(e)
    #    finally:
    #        conn.close()
    #
    # ---------------------------------
    #
    def update_item(self, element_id:int,
                    item:str, value:float|int):
        """ """
        1 / 0
        conn = create_connection(self.db_file)
        with conn:
            update_element_item(conn, element_id, item, value)
            #conn.commit()
    #
    @property
    def get_free_nodes(self):
        """
        find nodes not sharing elements
        """
        connectivities = self.get_connectivities
        #connectivities = [conn for conn in connectivities.values()]
        #column
        flat = list(chain.from_iterable(connectivities))
        return [k for k, v in Counter(flat).items() if v == 1]
    #
    # ---------------------------------
    #
    #@property
    def beam(self, values:None|list=None,
              df=None):
        """
        """
        items = self._beam(values, df)
        return self._beams
    #    if values:
    #        if isinstance(values, list):
    #            for item in values:
    #                element_name = item[0]
    #                try:
    #                    self._labels.index(element_name)
    #                    raise Exception('element {:} already exist'.format(element_name))
    #                except ValueError:
    #                    element_type = 'beam'
    #                    #mnumber = next(self.get_number())
    #                    # default
    #                    self._labels.append(element_name)
    #                    #self._number.append(mnumber)
    #                    self._type.append(element_type)
    #                    self._beams[element_name] = item[1:]
    #    #
    #    # dataframe input
    #    try:
    #        df.columns
    #        self._beams.df = df
    #    except AttributeError:
    #        pass
    #    #
    #    return self._beams
    #    #return ElementType(item_type='beam',
    #    #                   cls_type=self._beams, cls=self)
    #
    #    
    # ---------------------------------
    #
    #@property
    #def df(self):
    #    """nodes in dataframe format"""
    #    #print('elements df out')
    #    conn = create_connection(self.db_file)
    #    with conn:
    #        data = get_elements(conn, component=self._component)
    #    #
    #    header = ['name', 'type', 'material', 'section',
    #              'node_1', 'node_2', 'node_3', 'node_4',
    #              'roll_angle', 'title']
    #    return data[header]
    #
    #@df.setter
    #def df(self, df):
    #    """nodes in dataframe format"""
    #    #df = get_element_df(df)
    #    columns = list(df.columns)
    #    for key in columns:
    #        if re.match(r"\b((element(s)?(_|-|\s*)?)?type)\b", key, re.IGNORECASE):
    #            df['type'] = df[key].apply(lambda x: find_element_type(x))
    #            break
    #    #
    #    group = df.groupby("type")
    #    #
    #    try:
    #        group = group.get_group("beam")
    #        self._beams.df = group
    #    except KeyError:
    #        raise IOError('Element df not valid')
#  
#
#
class ElementType(Mapping):
    __slots__ =  ['_type', '_labels', #'_number',
                  '_cls_type', '_item_type']
    
    def __init__(self, item_type: str, cls_type, cls):
        """
        """
        self._cls_type = cls_type
        self._item_type = item_type
        self._labels = cls._labels
        self._type = cls._type
    #
    def __setitem__(self, item_name:str|int,
                    properties:list[str|float]) -> None:
        """
        item_name : element number
        properties : [material, section, node1, node2, roll_angle]
        """
        self._labels.append(item_name)
        self._type.append(self._item_type)
        self._cls_type[item_name] = properties
    
    def __getitem__(self, item_name:str|int):
        """
        """
        index = self._labels.index(item_name)       
        return self._cls_type[item_name]
    #
    #
    #
    def __len__(self) -> float:
        return len(self._cls_type._labels)

    def __iter__(self):
        """
        """
        return iter(self._cls_type._labels)

    def __contains__(self, value) -> bool:
        return value in self._cls_type._labels
    #
    #
    def __str__(self) -> str:
        """ """
        #print('--')
        return self._cls_type.__str__()
    #
    #
    @property
    def df(self):
        """ """
        return self._cls_type.df
    #
    #
    #@property
    def get_connectivities(self):
        """ """
        return self._cls_type.get_connectivities 
#
#
