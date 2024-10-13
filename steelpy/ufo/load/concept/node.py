#
# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
from array import array
from collections.abc import Mapping
import re
#
# package imports
# 
from steelpy.ufo.load.process.node import (NodeLoadBasic,
                                            get_nodal_load,
                                            PointNode)
from steelpy.utils.dataframe.main import DBframework


# ---------------------------------
#
#
class NodeLoadMaster(NodeLoadBasic):
    """
    FE Node Load class
    
    NodeLoad
        |_ name
        |_ number
        |_ type
        |_ complex
        |_ point [x, y, z, mx, my, mz]
        |_ acceleration [th[0], th[1], th[2],..., th[n]]
        |_ displacement [th[0], th[1], th[2],..., th[n]]
        |_ mass
    
    **Parameters**:  
      :number:  integer internal number 
      :name:  string node external name
    """
    __slots__ = ['_labels', '_title', '_complex', '_load_id', '_system',
                 '_type', '_fx', '_fy', '_fz', '_mx', '_my', '_mz']

    def __init__(self, load_type:str) -> None:
        """
        """
        super().__init__()
        #
        self._labels: list = []
        self._title: list = []
        self._complex: array = array("I", [])
        self._load_id: list = []
        self._system: array = array("I", [])        
        # real
        self._type = load_type
        self._fx: array = array('f', [])
        self._fy: array = array('f', [])
        self._fz: array = array('f', [])
        self._mx: array = array('f', [])
        self._my: array = array('f', [])
        self._mz: array = array('f', [])
    #
    #
    def __getitem__(self, node_name: int|str) -> list:
        """
        """
        idx_list: list = [x for x, _item in enumerate(self._labels)
                          if _item == node_name]
        #
        points: list = []
        if self._type in ['load']:
            for idx in idx_list:
                points.append(PointNode(self._fx[idx], self._fy[idx], self._fz[idx],
                                        self._mx[idx], self._my[idx], self._mz[idx],
                                        self._labels[idx], self._title[idx], self._load_id[idx],
                                        self._system[idx], self._complex[idx], self._type))
        
        elif self._type in ['displacement']:
            for idx in idx_list:
                points.append(DispNode(self._fx[idx], self._fy[idx], self._fz[idx],
                                       self._mx[idx], self._my[idx], self._mz[idx],
                                       self._labels[idx], self._title[idx], self._load_id[idx],
                                       self._system[idx], self._complex[idx], self._type))
        
        elif self._type in ['mass']:
            raise NotImplementedError()
        
        else:
            raise IOError(f'load type: {self._type} not valid')
            
        return points
    #
    #
    def __delitem__(self, node_name: int|str) -> None:
        """
        """
        indexes = [i for i, x in enumerate(self._labels)
                   if x == node_name]
        indexes.sort(reverse=True)
        for _index in indexes:
            self._fx.pop(_index)
            self._fy.pop(_index)
            self._fz.pop(_index)
            self._mx.pop(_index)
            self._my.pop(_index)
            self._mz.pop(_index)
            #
            self._labels.pop(_index)
            self._title.pop(_index)
            self._system.pop(_index)
            self._complex.pop(_index)
    #
    #
    @property
    def df(self):
        """nodes in dataframe format"""
        db = DBframework()
        # TODO : merge nodes
        #title = []
        data = {'load_name': self._load_id,
                'load_type': ['basic' for item in self._labels],
                'load_id': [idx + 1 for idx, item in enumerate(self._labels)],
                'load_system': self._system, #['global' if item == 0 else 'local'
                           #for item in self._system],
                'load_comment': self._title,
                'node_name':self._labels, 
                'Fx':self._fx, 'Fy':self._fy, 'Fz':self._fz, 
                'Mx':self._mx, 'My':self._my, 'Mz':self._mz}
        #
        # beam point case
        try:
            data.update({'L1': self._L1})
        except AttributeError:
            pass
        #
        # beam node load conversion
        try:
            data.update({'element_name': self._beam})
        except AttributeError:
            pass
        #      
        #
        df_nload = db.DataFrame(data=data, index=None)
        #df = df[['load_name', 'load_type', 'load_id', 'load_system', 'load_comment',
        #         'node_name', 'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']]          
        return df_nload

    @df.setter
    def df(self, values):
        """nodes in dataframe format"""
        #update
        self._labels.extend(values.node_name.tolist())
        self._load_id.extend(values.load_name.tolist())
        self._title.extend(values.load_title.tolist())
        self._system.extend([1 if item in ['local', 'member', 1] else 0
                             for item in values.system.tolist()])
        self._complex.extend([0 for item in values.system.tolist()])
        #
        self._fx.extend(values.Fx.tolist())
        self._fy.extend(values.Fy.tolist())
        self._fz.extend(values.Fz.tolist())
        self._mx.extend(values.Mx.tolist())
        self._my.extend(values.My.tolist())
        self._mz.extend(values.Mz.tolist())
        #
        # beam point case
        try:
            self._L1.extend(values.L1.tolist())
        except AttributeError:
            pass
        #
        # beam node load conversion
        try:
            self._beam.extend(values.element_name.tolist())
        except AttributeError:
            pass
        #
        #print('nodes df out')
    
#
#
# ---------------------------------
#
#
class NodeLoadItemIM(NodeLoadBasic):
    __slots__ = ['_title', '_labels', '_type', '_number',
                 '_f2u_points', #'_load_name',
                 '_load', '_displacement', '_mass', '_node']

    def __init__(self, load_name, load_title, points) -> None:
        """
        """
        super().__init__()
        #
        self._labels = []
        self._type = []
        #
        self._node = NodeItemIM(load_name=load_name,
                                load_title=load_title)
        self._f2u_points = points
    #
    def __setitem__(self, node_name: int|str,
                    node_load: list) -> None:
        """
        """
        load_type = node_load[0]
        #load_id = next(self.get_number())
        self._labels.append(node_name)
        #
        if re.match(r"\b(point|load|node)\b", load_type, re.IGNORECASE):
            #self.node_load._load = node_load[1:]
            self._node._load[node_name] = node_load[1:]
        
        elif re.match(r"\b(mass)\b", load_type, re.IGNORECASE):
            self._node._mass[node_name] = node_load[1:]
        
        elif re.match(r"\b(disp(lacement)?)\b", load_type, re.IGNORECASE):
            self._node._displacement[node_name] = node_load[1:]
        
        else:
            raise IOError(f'Load type {load_type} not recognized')
        #
    #
    def __getitem__(self, node_name: int | str):
        """
        """
        try:
            point = self._f2u_points[node_name]
            #node_name =self._f2u_points.get_name(node_name)
            if not node_name in self._labels:
                self._labels.append(node_name)
            #
            return self._node(node_name)
        except KeyError:
            raise IOError(f"Point {node_name} not found")

    #
    #def __contains__(self, value) -> bool:
    #    return value in self._labels

    #def __len__(self) -> int:
    #    return len(self._labels)

    #def __iter__(self):
    #    """
    #    """
    #    items = list(dict.fromkeys(self._labels))
    #    #items = list(set(self._labels))
    #    return iter(items)

    #
    def __str__(self, units: str = "si") -> str:
        """ """
        output = ""
        output += self._node.__str__()
        return output
    #
    #
    #
    @property
    def load(self):
        """
        """
        return self._node._load
    #
    #@load.setter
    #def load(self, value):
    #    """
    #    """
    #    self._load

#
# ---------------------------------
#
class NodeItemIM:
    __slots__ = ['_load', #'_load_name',
                 '_displacement', '_mass', '_node_id']

    def __init__(self, load_name:str|int, load_title:str):
        """
        """
        #self._load_name = load_name
        self._load = NodeLoadIM(load_name, load_title, "load")
        self._displacement = NodeLoadIM(load_name, load_title, "displacement")
        self._mass = NodeLoadIM(load_name, load_title, "mass")
    #
    def __call__(self, node_id):
        self._node_id = node_id
        return self
    #
    #
    @property
    def load(self):
        """
        """
        try:
            point_id = self._node_id
            return self._load[point_id]
        except :
            raise IndexError

    @load.setter
    def load(self, values: list):
        """
        Point Load
        """
        node_name = self._node_id
        if isinstance(values, dict):
            values.update({'type': 'load',})
            self._load[node_name] = values
        
        elif isinstance(values[0], list):
            for item in values:
                item.insert(0, 'load')
                self._load[node_name] = item
        
        else:
            values.insert(0, 'load')
            self._load[node_name] = values

    #
    @property
    def mass(self):
        """
        """
        point_id = self._node_id
        return self._mass[point_id]

    @mass.setter
    def mass(self, values: list):
        """
        """
        node_name = self._node_id
        if isinstance(values, dict):
            values.update({'type': 'mass',})
            self._mass[node_name] = values
            
        elif isinstance(values[0], list):
            for item in values:
                item.insert(0, 'mass')
                self._mass[node_name] = item
        else:
            values.insert(0, 'mass')
            self._mass[node_name] = values

    #
    @property
    def displacement(self):
        """
        """
        point_id = self._node_id
        return self._displacement[point_id]

    @displacement.setter
    def displacement(self, values: list):
        """
        """
        node_name = self._node_id
        if isinstance(values, dict):
            values.update({'type': 'displacement',})
            self._displacement[node_name] = values
            
        elif isinstance(values[0], list):
            for item in values:
                item.insert(0, 'displacement')
                self._displacement[node_name] = item
        else:
            values.insert(0, 'displacement')
            self._displacement[node_name] = values

    #   
    #
    def __str__(self, units: str = "si") -> str:
        """ """
        output = ""
        # output += "--- Nodal Load \n"
        output += self._load.__str__()
        output += self._mass.__str__()
        output += self._displacement.__str__()
        return output
#
# ---------------------------------
#
class NodeLoadIM(NodeLoadMaster):
    """
    FE Node Load class
    
    NodeLoad
        |_ name
        |_ number
        |_ type
        |_ complex
        |_ point [x, y, z, mx, my, mz]
        |_ acceleration [th[0], th[1], th[2],..., th[n]]
        |_ displacement [th[0], th[1], th[2],..., th[n]]
        |_ mass
    
    **Parameters**:  
      :number:  integer internal number 
      :name:  string node external name
    """
    __slots__ = ['_labels', '_title', '_complex', '_load_id', '_system',
                 '_type', '_fx', '_fy', '_fz', '_mx', '_my', '_mz',
                 '_load_name', '_load_title', '_system_flag']

    def __init__(self, load_name:str|int, load_title:str, 
                 load_type:str) -> None:
        """
        """
        super().__init__(load_type)
        self._load_name = load_name
        self._load_title = load_title
        #self._system_flag: int = 0 # Global system

    #
    def __setitem__(self, node_name: int|str,
                    point_load: list) -> None:
        """
        """
        #
        self._labels.append(node_name)
        self._load_id.append(self._load_name)
        #
        point_load = get_nodal_load(point_load)
        load_type = point_load.pop(0)
        load_title = point_load.pop(-1)        
        #
        self._system.append(self._system_flag)
        self._complex.append(0)
        #
        # point_load = get_nodal_load(point_load)
        self._fx.append(point_load[0])
        self._fy.append(point_load[1])
        self._fz.append(point_load[2])
        self._mx.append(point_load[3])
        self._my.append(point_load[4])
        self._mz.append(point_load[5])
        #
        self._title.append(load_title)
        #self._type.append(load_type)
        # print('--')
    #
    # @property
    # def get_items_sum(self):
    def get_group_nodes(self):
        """
        """
        items = []
        set_items = set(self._labels)
        for item in set_items:
            _index = self._labels.index(item)
            rows = self.__getitem__(item)
            sum_load = [0, 0, 0, 0, 0, 0]
            for _load in rows:
                sum_load[0] += _load.fx
                sum_load[1] += _load.fy
                sum_load[2] += _load.fz
                sum_load[3] += _load.mx
                sum_load[4] += _load.my
                sum_load[5] += _load.mz
            items.append(PointNode(*sum_load,
                                   self._labels[_index],
                                   self._title[_index],
                                   self._system[_index],
                                   self._complex[_index]))
        return items
    #
    # def update(self, other) -> None:
    #    """
    #    """
    #    #1/0
    #    pl = other._point
    #    try:
    #        #_name = pl.name
    #        index = pl._index
    #        # update
    #        self._labels.extend(pl._labels)
    #        self._title.extend(pl._title)
    #        self._system.extend(pl._system)
    #        self._distance.extend(pl._distance)
    #        self._complex.extend(pl._complex)
    #        #
    #        self._fx.extend(pl._fx)
    #        self._fy.extend(pl._fy)
    #        self._fz.extend(pl._fz)
    #        self._mx.extend(pl._mx)
    #        self._my.extend(pl._my)
    #        self._mz.extend(pl._mz)
    #        # print('?????')
    #    except AttributeError:
    #        pass
    #
#
