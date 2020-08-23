# 
# Copyright (c) 2009-2020 fem2ufo
#

# Python stdlib imports
from dataclasses import dataclass
from typing import Tuple, Dict, List, ClassVar, Union

# package imports
from steelpy.f2uModel.mesh.main import Mesh
from steelpy.f2uModel.concept.main import Concepts
from steelpy.f2uModel.sections.main import Sections
from steelpy.f2uModel.material.main import Materials
from steelpy.f2uModel.properties.main import Properties
from steelpy.f2uModel.load.main import Loading
from steelpy.f2uModel.sql.main import f2uDB

#
def get_number(items: dict) -> int:
    """
    return maximum consecutive number of items
    """
    _number = [item.number for item in items.values()]
    try:
        return max(_number)
    except ValueError:
        return 0


# GEOMETRY Section
class f2uModel:
    """
    FE Geometry model class
    
    Parameters:
      :number: integer
      :name: string
      :nodes: list of node's class
      :elements: list of element's class
      :materials: list of material's class
      :sets: list of groups (elements, nodes, etc)
      :sections: list of section's class
      :vectors: list of guide points
      :eccentricities: list of eccentricities
      :joints: list of joint's class
      :hinges: list of hinges definitios
      :loads: list of load's class
      :data: FE model data
      :units: FE model units
      :soil: list of soil's class
      :hydrodynamics: hydrodynamic's data
      :boundaries: list of FE model boundary
    
    Parameters:  
      :number:  integer internal number 
      :name:  string node external name
    """
    __slots__ = ('component', 'type', 'data',
                 '_materials', '_sections', '_properties',
                 'sets', 'mesh', '_concept', '_sql',
                 #'_material_default', '_section_default', '_set_dafault', 
                 '_boundaries', '_nodes', '_load')

    def __init__(self, component: str) -> None:
        """
        """
        print("-- module : fem2ufo Version 5.00dev")
        print('{:}'.format(52*'-'))
        #
        #self.number: int = number
        self.component: str = component
        self.type: str = 'substructure'
        # set main components
        self._materials: ClassVar = Materials()
        self._sections: ClassVar = Sections()
        self._properties: ClassVar = Properties()
        # set mesh components
        self.mesh: ClassVar = Mesh()
        self.mesh.materials = self._materials
        self.mesh.sections = self._sections
        # loads component
        self._load = Loading()
        self._sql = f2uDB(model=self)
        #
        #self.sets: ClassVar = Groups(concepts=self._concept,
        #                             mesh= self.mesh)
        #self.data: Tuple = InputData()
        #
        # set concepts
        self._concept: ClassVar = Concepts(mesh= self.mesh,
                                           #load = self._load,
                                           properties= self._properties)
        #        
        #
        # start defaults
        #self._material_default: bool = None
        #self._section_default: bool = None
        #self._set_dafault: bool = None

    #
    @property
    def materials(self) -> ClassVar:
        """
        """
        return self._materials   

    #
    @property
    def sections(self) -> ClassVar:
        """
        """
        return self._sections

    #
    @property
    def properties(self) -> ClassVar:
        """
        """
        return self._properties
    #
    @property
    def concept(self) -> ClassVar:
        """
        """
        return self._concept

    #
    @property
    def groups(self):
        """ """
        return self.sets
    #
    @property
    def default_material(self, material_name: str) -> None:
        """
        """
        _material_default = self.material[material_name]

    #
    @property
    def default_section(self, section_name: str) -> None:
        """
        """
        _section_dafault = self.sets[section_name]

    #
    @property
    def default_group(self, group_name: str) -> None:
        """
        """
        _set_dafault = self.sets[group_name]
    #
    #
    @property
    def load(self) -> ClassVar:
        """
        """
        return self._load
    #
    @property
    def sql(self):
        """
        """
        return self._sql    
    #
    #
    #
    #
    def get_mesh(self) -> None:
        """
        """
        print('--- Meshing started')
        self._set_mesh()
        #self._set_boundary()
        self._set_load()
        print('--- Mesh Completed')
    #
    #
    def _set_mesh(self):
        """ """
        mesh = self.mesh
        number = mesh.elements.get_number()
        for key, memb in self._concept.beam.items():
            #print(key)
            total_length = memb.length
            _nodes = memb.connectivity
            node_res = _nodes[0].name # start
            node_end = _nodes[1].name
            for step in memb.step:
                _mnumber = next(number)
                step._mesh = _mnumber
                #print(_mnumber, key)
                try:
                    1/step.length.value
                    total_length -= step.length
                    coord = memb.find_coordinate(step.length)
                    new_node = self._concept.points.get_point_name(coord)
                    # elements [node1, node2, material, section]
                    mesh.elements[_mnumber] = ['beam', node_res, new_node,
                                               step.material.name, step.section.name]
                    node_res = new_node
                except ZeroDivisionError:
                    # elements [node1, node2, material, section]
                    mesh.elements[_mnumber] = ['beam', node_res, node_end,
                                               step.material.name, step.section.name]
            #print('-->')
        #
        print('end')
    #
    def _set_boundary(self):
        """ """
        mesh = self.mesh
        nodes = mesh._nodes
        cboundary = self._concept.boundary
        missing = []
        for key, value in cboundary._supports.items():
            key, value
        print ( '-->' )
    #
    def _set_load(self):
        """ """
        mesh = self.mesh
        nodes = mesh.nodes
        basic_load = self._load.basic
        concept_bload = self._concept.load._basic
        for load_name, lcase in concept_bload.items():
            basic_load[load_name] = lcase.title
            # line load process
            for bname, loads in lcase.beam.line_load:
                beam = self._concept.beam[bname]
                Lb = beam.length.value
                print('---> ',load_name, bname)
                #
                for load in loads:
                    label = load.name
                    waxial = linefit(load.qx1, load.qx2,
                                     Lb, load.L1, load.L2)
                    winplane = linefit(load.qy1, load.qy2,
                                       Lb, load.L1, load.L2)
                    woutplane = linefit(load.qz1, load.qz2,
                                        Lb, load.L1, load.L2)
                    # start loop beam steps
                    xi = 0
                    for step in beam.step:
                        elem_name = step._mesh
                        element = mesh.elements[elem_name]
                        Lbi = element.length_node2node(nodes)
                        xi += Lbi
                        qaxial = waxial.qi(xi)
                        qinp = winplane.qi(xi)
                        qoutp = woutplane.qi(xi)
                        # check load on segment
                        try:
                            Li = winplane.Li(xi, Lbi)
                        except RuntimeWarning:
                            continue # no load should be applied to this segment
                        # set load for mesh element
                        print(elem_name, label, Lb, xi, qinp, qoutp, Li)
                        basic_load[load_name].line_beam[elem_name]  = [qaxial[0], qinp[0], qoutp[0],
                                                                       qaxial[1], qinp[1], qoutp[1],
                                                                       Li[0], Li[1]]
            #
            # Point load process
            for bname, loads in lcase.beam.point_load:
                beam = self._concept.beam[bname]
                Lb = beam.length.value
                print('---> ',load_name, bname)
                for load in loads:
                    label = load.name
                    L1 = load.distance
                    pload = pointfit(Lb, L1)
                    # start loop beam steps
                    xi = 0                    
                    for step in beam.step:
                        elem_name = step._mesh
                        element = mesh.elements[elem_name]
                        Lbi = element.length_node2node(nodes)
                        xi += Lbi
                        # check load on segment
                        try:
                            Li = pload.Li(xi, Lbi)
                        except RuntimeWarning:
                            continue # no load should be applied to this segment                        
                        # set load for mesh element
                        print(elem_name, label, Lb, xi)
                        basic_load[load_name].point_beam[elem_name] = [Li, *load[:6]]
            #
        #
        #print('-->')
#
#
@dataclass
class linefit:
    __slots__ = ['q1', 'q2', 'L', 'L1', 'L2',
                 'L3', 'Lstart', 'Lstop', '_qi']

    def __init__(self, q1:float, q2:float,
                 L:float, L1:float, L2:float) -> None:
        """ """
        self.q1:float = q1
        self.q2:float = q2
        self._qi:float = q1
        #
        self.L:float = L
        self.L1:float = L1
        self.L2:float = L2
        self.L3 = self.L - self.L2
        self.Lstop = self.L - self.L2
    #
    @property
    def slope(self) -> float:
        """ """
        return (self.q2-self.q1)/(self.L3-self.L1)
    #
    def qi(self, x:float) -> List[float]:
        """ """
        q1 = self._qi
        if x > self.L3:
            self._qi = self.q2
            #q2 = round(self.q1 + self.slope * (self.L3-self.L1), 3)
        else:
            self._qi = round(self.q1 + self.slope * (x-self.L1), 3)
        #self._qi = q2
        return [q1, self._qi]
    #
    def Li(self, x:float, Lb:float) -> Union[Exception,List[float]]:
        """ """
        try:
            1/(self.L1 + self.L2)
            if x < self.L1: # no load for this step
                raise RuntimeWarning
            else:
                try:
                    1 / self.Lstop
                    try:
                        Lstart = self.Lstart
                    except AttributeError:
                        Lstart = Lb - (x - self.L1)
                        self.Lstart = 0
                    #
                    if x > self.L3:
                        self.Lstop = 0
                        return [Lstart, x - self.L3]
                    else:
                        return [Lstart, 0]
                except ZeroDivisionError: # no load after this step
                    raise RuntimeWarning
        except ZeroDivisionError:
            return [0, 0]
#
@dataclass
class pointfit:
    __slots__ = ['L', 'L1', 'Lstop']
    
    def __init__(self, L:float, L1:float) -> None:
        """ """
        self.L:float = L
        self.L1:float = L1
        self.Lstop:float = L1
    #
    def Li(self, x:float, Lb:float) -> Union[Exception,float]:
        """ """
        if x < self.L1: # no load for this step
            raise RuntimeWarning
        else:
            try:
                1 / self.Lstop
                self.Lstop = 0
                return Lb - (x - self.L1)
            except ZeroDivisionError:  # no load after this step
                raise RuntimeWarning

#