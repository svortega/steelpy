# 
# Copyright (c) 2009-2020 fem2ufo
# 

# Python stdlib imports
from typing import Dict, Union, List


# package imports


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
class BeamLine:
    
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
        beam = [ self._cls._beams[_index] for _index, item in enumerate(self._cls._labels)
                 if self._cls._labels[_index] == load_name 
                 and self._cls._load_case[_index] == self._load_name
                 and self._cls._beams[_index] != -1]
        return self._beam[beam[-1]]
        #print('---')


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
        return BeamLine(self._cls, self._load_name, beam_load_type)
    
    #
    #
    @property
    def line_load(self):
        """
        """
        beam_load_type = f2u_load.basic[self._load_name]._beam_line
        return BeamLine(self._cls, self._load_name, beam_load_type)
#
#
class LoadTypes:
    """
    """
    __slots__ = ["_cls", "_load_name"]
    
    def __init__(self, cls, load_name):
        """
        """
        self._cls = cls
        self._load_name = load_name
    #
    @property
    def point(self):
        """
        """
        return PointLoadConcept(self._cls, self._load_name)
    #
    @point.setter
    def point(self, coord):
        """
        """
        node_name = f2u_points.get_point_name(coord)
        self._cls._nodes.append(node_name)
        self._cls._load_case.append(self._load_name)
    #
    @property
    def beam(self):
        """
        """
        return BeamLoadConcept(self._cls, self._load_name)
    
    @beam.setter
    def beam(self, element):
        """
        """
        self._cls._beams.append(element.name)
        self._cls._load_case.append(self._load_name)   
    
#
class ConceptBasicLoad:
    
    __slots__ = ["_basic", "_load_case", 
                 "_labels", "_nodes", "_beams"]
    
    def __init__(self) -> None:
        """
        """
        self._load_case:List = []
        
        self._labels:List = []
        self._nodes:List = []
        self._beams:List = []
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
class ConcepLoad:
    
    __slots__ = ["_basic", "f2u_load", "f2u_points", "f2u_concepts"]
    
    def __init__(self, load, points, concepts) -> None:
        """
        """
        global f2u_load, f2u_points, f2u_concepts
        f2u_load = load
        f2u_points = points
        f2u_concepts = concepts
        #
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