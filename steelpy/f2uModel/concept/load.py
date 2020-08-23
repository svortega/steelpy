# 
# Copyright (c) 2009-2020 fem2ufo
# 

# Python stdlib imports
from collections.abc import Mapping
from typing import Dict, Union, List


# package imports
from steelpy.f2uModel.load.main import Loading


class PointLoad:
    
    __slots__ = ["_cls", "_point", "_load_name"]
    
    def __init__(self, cls, load_name, point_load_type) -> None:
        """
        """
        self._cls = cls
        self._point = point_load_type
        self._load_name = load_name
    #
    #
    def __setitem__(self, load_name:Union[str,int],
                    load) -> None:
        """
        """
        node_name = self._cls._nodes[-1]
        self._point[node_name] = load
        #
        self._cls._labels.append(load_name)
        index = len(self._cls._labels) - 1
        try:
            self._cls._nodes[index]
        except IndexError:
            node = self._cls._nodes[-1]
            self._cls._nodes.append(node)
        print('---')
        self._cls._beams.append(-1)
        #self._cls._load_case.append(self._load_name)
    
    def __getitem__(self, load_name:Union[str,int]):
        """
        """
        #print('---')
        nodes = [ self._cls._nodes[_index] for _index, item in enumerate(self._cls._labels)
                  if self._cls._labels[_index] == load_name 
                  and self._cls._load_case[_index] == self._load_name
                  and self._cls._nodes[_index] != -1]
        return self._point[nodes[-1]]
        #return res
    

class PointLoadConcept:
    
    __slots__ = ["_cls", "_load_name"]
    
    def __init__(self, cls, load_name) -> None:
        """
        """
        self._cls = cls
        self._load_name = load_name
    
    #
    @property
    def load(self):
        """
        """
        point = f2u_load.basic[self._load_name]._nodal_load
        return PointLoad(self._cls, self._load_name, point)
    #
    @property
    def mass(self):
        """
        """
        return PointLoad(self._cls)
    

#
#
#
class BeamLoadSet:
    
    __slots__ = ["_cls", "_beam", "_load_name"]
    
    def __init__(self, cls, load_name, beam_load_type) -> None:
        """
        """
        self._cls = cls
        self._beam = beam_load_type
        self._load_name = load_name
    #
    def __setitem__(self, load_name:Union[str,int],
                    load) -> None:
        """
        """
        beam_name = self._cls._beams[-1]
        self._beam.coordinate_system = self._cls._system_flag
        self._beam[beam_name] = load
        self._beam.name = load_name
        self._cls._labels.append(load_name)
        index = len(self._cls._labels) - 1
        try:
            self._cls._beams[index]
        except IndexError:
            beam = self._cls._beams[-1]
            self._cls._beams.append(beam)
        #print('---')
        self._cls._nodes.append(-1)
    
    def __getitem__(self, load_name:Union[str,int]):
        """
        """
        beam = [self._cls._beams[_index] for _index, item in enumerate(self._cls._labels)
                if self._cls._labels[_index] == load_name 
                and self._cls._load_case[_index] == self._load_name
                and self._cls._beams[_index] != -1]
        return self._beam[beam[-1]]
        #print('---')
    #
    def __iter__(self):
        # return self._beam.__iter__()
        for key, item in self._beam.items():
            yield key, item


#
class BeamLoadConcept:
    
    __slots__ = ["_cls", "_load_name"]
    
    def __init__(self, cls, load_name) -> None:
        """
        """
        self._cls = cls
        self._load_name = load_name
    #
    @property
    def point_load(self):
        """
        """
        beam_load_type = f2u_load.basic[self._load_name]._beam_point
        return BeamLoadSet(self._cls, self._load_name, beam_load_type)
    
    #
    #
    @property
    def line_load(self):
        """
        """
        beam_load_type = f2u_load.basic[self._load_name]._beam_line
        return BeamLoadSet(self._cls, self._load_name, beam_load_type)
    #
    #
    @property
    def coordinate_system(self):
        if self._cls._system_flag != 0:
            return "local"
        return "global"
    
    @coordinate_system.setter
    def coordinate_system(self, system:Union[str,int]):
        """
        Coordinate system for load : global or local (member)
        """
        self._cls._system_flag = 0
        if system in ['local', 'member', 1]:
            self._cls._system_flag = 1
    #    
#
#
class LoadTypes:
    """
    """
    __slots__ = ["_cls", "name", "title", "number"]
    
    def __init__(self, cls, load_name):
        """
        """
        self._cls = cls
        self.name = load_name
        basic = f2u_load.basic[load_name]
        self.title = basic.title
        self.number = basic.number
    #
    @property
    def point(self):
        """
        """
        return PointLoadConcept(self._cls, self.name)
    #
    @point.setter
    def point(self, coord):
        """
        """
        node_name = f2u_points.get_point_name(coord)
        self._cls._nodes.append(node_name)
        self._cls._load_case.append(self.name)
    #
    @property
    def beam(self):
        """
        """
        return BeamLoadConcept(self._cls, self.name)
    
    @beam.setter
    def beam(self, element):
        """
        """
        self._cls._system_flag = 0 
        self._cls._beams.append(element.name)
        self._cls._load_case.append(self.name)   
    
#
class ConceptBasicLoad(Mapping):
    
    __slots__ = ["_basic", "_load_case", "_system_flag",
                 "_labels", "_nodes", "_beams", "_bload"]
    
    def __init__(self) -> None:
        """
        """
        self._load_case:List = []
        self._labels:List = []
        self._nodes:List = []
        self._beams:List = []
        # 0-global/ 1-local
        self._system_flag:int = 0           
    #
    def __setitem__(self, load_name:Union[str,int],
                    load_title:str) -> None:
        """
        """
        self._load_case.append(load_name)
        f2u_load.basic[load_name] = load_title
    
    def __getitem__(self, load_name:Union[str,int]):
        """ """
        return LoadTypes(self, load_name)
    #
    def __len__(self) -> float:
        return len(self._labels)

    def __iter__(self) : #-> Iterator
        """
        """
        return f2u_load.basic.__iter__()
        #for basic in f2u_load.basic:
        #    basic
        #return iter(self._labels)

    def __contains__(self, value) -> bool:
        return value in self._load_case    
#
class ConcepLoad:
    
    __slots__ = ["_basic", "f2u_load", "f2u_points", "f2u_concepts"]
    
    def __init__(self, points, concepts) -> None: # load, 
        """
        """
        global f2u_points, f2u_concepts, f2u_load
        f2u_load = Loading() # TODO: change name
        f2u_points = points
        f2u_concepts = concepts
        #
        #self._bload = Loading()
        self._basic = ConceptBasicLoad()
    #
    @property
    def basic(self):
        """
        """
        return self._basic
    #
    #@property
    #def combination(self):
    #    """
    #    """
    #    return self._combination
    #    