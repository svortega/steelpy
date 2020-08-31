# 
# Copyright (c) 2009-2020 fem2ufo
#

# Python stdlib imports
from dataclasses import dataclass
from typing import Tuple, Dict, List, ClassVar, Union

# package imports
from steelpy.process.units.main import Units
from steelpy.process.geometry.L3D import DistancePointLine3D

#
#
class Meshing:
    __slots__ = ["mesh", "concept", "load"]
    
    def __init__(self, mesh, concept, load) -> None:
        """
        """
        self.mesh = mesh
        self.concept = concept
        self.load = load
    #
    def get_mesh(self) -> None:
        """
        """
        #print('-- Meshing Component: {:}'.format(self.component))
        print('-- Meshing Component')
        self._set_mesh()
        self._set_boundary()
        self._set_load()
        print('-- Meshing Completed')
    #
    def _set_mesh(self):
        """ """
        print ( '--- Meshing Concepts' )
        mesh = self.mesh
        elements = mesh.elements
        elem_number = mesh.elements.get_number()
        cbeams = self.concept.beam
        cpoints = self.concept.points
        for key, beam in cbeams.items():
            total_length = beam.length
            _nodes = beam.connectivity
            node_res = _nodes[0].name # start
            node_end = _nodes[1].name
            for step in beam.step:
                _mnumber = next(elem_number)
                step._mesh = _mnumber
                print("concept:", key, "element :", _mnumber)
                try:
                    1/step.length.value
                    total_length -= step.length
                    coord = beam.find_coordinate(step.length)
                    try:
                        new_node = cpoints.get_point_name(coord)
                    except IOError:
                        new_node = cpoints.get_new_point(coord)
                    # elements [node1, node2, material, section]
                    elements[_mnumber] = ['beam', node_res, new_node,
                                          step.material.name, step.section.name]
                    node_res = new_node
                except ZeroDivisionError:
                    # elements [node1, node2, material, section]
                    elements[_mnumber] = ['beam', node_res, node_end,
                                          step.material.name, step.section.name]
            #print('-->')
        #
        #print('end')
    #
    def _set_boundary(self):
        """ """
        print ( '--- Meshing Boundaries' )
        units = Units()
        mesh = self.mesh
        nodes = mesh._nodes
        elements = mesh.elements
        elem_number = elements.get_number()
        boundaries = mesh._boundaries
        # concepts
        cboundary = self.concept.boundary
        cbeams = self.concept.beam
        cpoints = self.concept.points
        missing = []
        # find existing nodes
        for key, value in cboundary.items():
            support = value.support
            point = value.point
            try:
                node_id = nodes.get_point_name(point)
                boundaries.node[node_id] = support
                print("Boundary: ", key, " @ Node: ", node_id)
            except IOError:
                missing.append(key)
        #
        # if missing boundaries, find if coordinates along members
        if missing:
            for key, beam in cbeams.items():
                mnodes = beam.connectivity
                p1 = [mnodes[0].x.value, mnodes[0].y.value, mnodes[0].z.value]
                p2 = [mnodes[1].x.value, mnodes[1].y.value, mnodes[1].z.value]
                point_line = DistancePointLine3D(p1, p2)
                missing_found = []
                for boundary in missing:
                    Pb = cboundary[boundary].point
                    if point_line.is_on_segment(Pb):
                        left_dist = point_line.left_dist
                        #print(key, point_line.position, point_line.left_dist, point_line.right_dist)
                        missing_found.append(boundary)
                        total_length = 0
                        step_no = len(beam.step)
                        for step in beam.step:
                            memb = elements[step._mesh]
                            total_length += memb.length_node2node(nodes)
                            # TODO: capture warning when point misaligned due to tolerances
                            if total_length < left_dist:
                                continue
                            # get node coordinate
                            coord = beam.find_coordinate(left_dist*units.m)
                            try:
                                new_node = cpoints.get_point_name(coord)
                            except IOError:
                                new_node = cpoints.get_new_point(coord)
                            # set boundary
                            support = cboundary[boundary].support
                            boundaries.node[new_node] = support
                            print ( "Boundary: ", boundary," on Beam: ",key, " @ Node: ", new_node)
                            # existing element
                            mnodes = memb.connectivity
                            node_end = mnodes[-1]
                            mnodes[-1] = new_node
                            memb.connectivity = mnodes
                            # new element
                            _mnumber = next(elem_number)
                            elements[_mnumber] = ['beam', new_node, node_end,
                                                  step.material.name, step.section.name ]
                            # introduce new concept beam step to boundary coord
                            step_no += 1
                            cbeams[key].step[step_no].length = left_dist * units.m
                            cbeams[key].step[step_no].material = step.material.name
                            cbeams[key].step[step_no].section = step.section.name
                            cbeams[key].step[step_no]._mesh = _mnumber
                            print("concept:", key, "element :", _mnumber)
                            break
                        continue
                    #try:
                    #    1/test.on_segment
                    #
                    #except ZeroDivisionError:
                    #    print("--> off")
                    #print(test.dist)
                    #print ( "--> off" )
                #
                for item in missing_found:
                    missing.remove(item)
                # TODO: capture if missing not empty
                if not missing:
                    break
        #print('-->')
    #
    def _set_load(self):
        """ """
        print ( '--- Meshing Basic Load' )
        mesh = self.mesh
        nodes = mesh.nodes
        basic_load = self.load.basic
        concept_bload = self.concept.load._basic
        for load_name, lcase in concept_bload.items():
            basic_load[load_name] = lcase.title
            # Beam line load process
            for bname, loads in lcase.beam.line_load:
                beam = self.concept.beam[bname]
                Lb = beam.length.value
                print('---> ',load_name, bname)
                #
                for load in loads:
                    label = load.load_name
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
                        basic_load[load_name].line_beam.name = label
            # Beam point load process
            for bname, loads in lcase.beam.point_load:
                beam = self.concept.beam[bname]
                Lb = beam.length.value
                print('---> ',load_name, bname)
                for load in loads:
                    label = load.load_name
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
                        basic_load[ load_name ].point_beam.name = label
            # Nodal load process
            #pl = concept_bload.point
            basic_load[load_name].point_node.update(lcase.point.load)
        #
        #print('-->')
    #
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