#
# Copyright (c) 2019-2023 steelpy
#

# Python stdlib imports
from __future__ import annotations

import math
from dataclasses import dataclass
#from typing import Tuple, Dict, List, ClassVar, Union

# package imports
from steelpy.process.units.main import Units
from steelpy.process.geometry.L3D import DistancePointLine3D
from steelpy.f2uModel.mesh.main import Mesh
#
#
class Meshing:
    __slots__ = ["_mesh", "concept"]
    
    def __init__(self, concept, mesh) -> None:
        """
        """
        #self.mesh = mesh
        self.concept = concept
        self._mesh = mesh
        #materials = concept
        #self._mesh = Mesh(materials=concept.materials(),
        #                  sections=concept.sections(),
        #                  mesh_type=mesh_type,
        #                  db_file=db_file)
        #self._load = mesh.load()
        #print('--')
    #
    def get_mesh(self) -> None:
        """
        """
        print('-- Meshing Component')
        self._set_mesh()
        self._set_boundary()
        self._set_load()
        print('-- Meshing Completed')
        #return self._mesh
    #
    def _set_mesh(self):
        """ """
        print('--- Meshing Concepts')
        mesh = self._mesh
        elements = mesh.elements()
        #elem_number = elements.get_number()
        cbeams = self.concept.beams()
        #
        for key, beam in cbeams.items():
            total_length = beam.length
            p1, p2 = beam.connectivity
            node_res = self._get_node_name(p1[:3])
            node_end = self._get_node_name(p2[:3])
            for step in beam.step:
                mnumber = next(elements.get_number())
                step._mesh = mnumber
                print(f"concept: {key} --> element: {mnumber}")
                try:
                    1/step.length.value
                    total_length -= step.length
                    coord = beam.find_coordinate(step.length)
                    new_node = self._get_node_name(coord)
                    # elements [node1, node2, material, section, beta, title]
                    elements[mnumber] = ['beam', node_res, new_node,
                                          step.material.name, 
                                          step.section.name, 
                                          beam.beta, key]
                    node_res = new_node
                except ZeroDivisionError:
                    # elements [node1, node2, material, section, beta, title]
                    elements[mnumber] = ['beam', node_res, node_end,
                                          step.material.name, 
                                          step.section.name, 
                                          beam.beta, key]
            #print('-->')
        #print('end meshing')
    #
    def _get_node_name(self, coord):
        """ """
        nodes = self._mesh.nodes()
        try:
            return nodes.get_node_name(coord)
        except IOError:
            return nodes.get_new_node(coord)
    #
    def _set_boundary(self):
        """ """
        print('--- Meshing Boundaries')
        units = Units()
        # Mesh
        mesh = self._mesh
        nodes = mesh.nodes()
        elements = mesh.elements()
        boundaries = mesh.boundaries()
        supports = boundaries.supports()
        # concepts
        cboundary = self.concept.boundaries()
        cbeams = self.concept.beams()
        #
        missing = []
        # find existing nodes
        for key, value in cboundary.items():
            support = value.support
            point = value.point
            try:
                node_id = nodes.get_node_name(point)
                boundaries.node[node_id] = [*support[:6], key]
                print(f"Boundary: {key}  @ Node: {node_id}")
            except IOError:
                missing.append(key)
        #
        # if missing boundaries, find if coordinates along members
        if missing:
            for key, beam in cbeams.items():
                p1, p2 = beam.connectivity
                point_line = DistancePointLine3D(p1[:3], p2[:3])
                missing_found = []
                for boundary in missing:
                    Pb = cboundary[boundary].point
                    # Fixme : what format?
                    Pb = [Pb[0].value, Pb[1].value, Pb[2].value]
                    if point_line.is_on_segment(Pb):
                        left_dist = point_line.left_dist
                        #print(key, point_line.position, point_line.left_dist, point_line.right_dist)
                        missing_found.append(boundary)
                        total_length = 0
                        step_no = len(beam.step)
                        for step in beam.step:
                            memb = elements[step._mesh]
                            total_length += memb.L
                            # TODO: capture warning when point misaligned due to tolerances
                            print(f"beam length step: {total_length:1.4e} left: {left_dist:1.4e}")
                            if total_length < left_dist:
                                continue
                            # get node coordinate
                            coord = beam.find_coordinate(left_dist*units.m)
                            new_node = self._get_node_name(coord)
                            #try:
                            #    new_node = cpoints.get_point_name(coord)
                            #except IOError:
                            #    new_node = cpoints.get_new_point(coord)
                            # set boundary
                            support = cboundary[boundary].support
                            if support:
                                supports[new_node] = support
                                #boundaries.node[new_node] = support
                                print(f"Boundary: {boundary} on Beam: {key} @ Node: {new_node}")
                            # existing element
                            mnodes = memb.connectivity
                            node_end = mnodes[-1]
                            mnodes[-1] = new_node
                            memb.connectivity = mnodes
                            # new element
                            mnumber = next(elements.get_number())
                            elements[mnumber] = ['beam', new_node, node_end,
                                                  step.material.name,
                                                  step.section.name, beam.beta]
                            # introduce new concept beam step to boundary coord
                            step_no += 1
                            cbeams[key].step[step_no].length = left_dist * units.m
                            cbeams[key].step[step_no].material = step.material.name
                            cbeams[key].step[step_no].section = step.section.name
                            cbeams[key].step[step_no]._mesh = mnumber
                            print(f"concept: {key} --> element: {mnumber}")
                            break
                        continue
                #
                for item in missing_found:
                    missing.remove(item)
                # TODO: capture if missing not empty
                if not missing:
                    break
        #print(' end meshing boundary')
    #
    def _set_load(self):
        """ """
        print( '--- Meshing Basic Load' )
        # Mesh
        mesh = self._mesh
        nodes = mesh.nodes()
        elements = mesh.elements()
        # Load
        load = mesh.load()
        basic_load = load.basic()
        # Concept
        Concept = self.concept
        concept_load = Concept.load()
        Clbasic = concept_load.basic()
        Points = Concept.points()
        #
        for load_name, lcase in Clbasic.items():
            basic_load[load_name] = lcase.title
            node_load = basic_load[load_name].node()
            beam_load =  basic_load[load_name].beam()
            # Beam load process
            for bname, loads in lcase.beams.items():
                beam = self.concept.beam[bname]
                Lb = beam.length.value
                print(f'---> Load: {load_name} Beam: {bname} L: {Lb:4.2f}')
                # Beam line load process
                for load in loads.line:
                    label = load.load_name
                    print(f'Load Title: {label} - Line Load')
                    waxial = linefit(load.qx0, load.qx1,
                                     Lb, load.L0, load.L1)
                    winplane = linefit(load.qy0, load.qy1,
                                       Lb, load.L0, load.L1)
                    woutplane = linefit(load.qz0, load.qz1,
                                        Lb, load.L0, load.L1)
                    # start loop beam steps
                    xi = 0
                    for step in beam.step:
                        elem_name = step._mesh
                        element = elements[elem_name]
                        Lbi = element.length
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
                        print(f'Element: {elem_name} {Lbi:4.2f} {xi:4.2f} {qinp} {qoutp} {Li}')
                        beam_load[elem_name].line = [qaxial[0], qinp[0], qoutp[0],
                                                     qaxial[1], qinp[1], qoutp[1],
                                                     Li[0], Li[1], label]
                #
                # Beam point load process
                for load in loads.point:
                    label = load.load_name
                    print(f'Load Title: {label} - Point Load')
                    L1 = load.distance
                    pload = pointfit(Lb, L1)
                    # start loop beam steps
                    xi = 0                    
                    for step in beam.step:
                        elem_name = step._mesh
                        element = elements[elem_name]
                        Lbi = element.length
                        xi += Lbi
                        # check load on segment
                        try:
                            Li = pload.Li(xi, Lbi)
                        except RuntimeWarning:
                            continue # no load should be applied to this segment                        
                        # set load for mesh element
                        print(f'Element: {elem_name} L1: {Li:4.2f} {load[:6]}')
                        beam_load[elem_name].point = [Li, *load[:6], label]
            #
            # Nodal load process
            for node_id, loads in lcase.points.items():
                point = Points[node_id]
                try: # if node exist
                    node_name = nodes.get_node_name(point[:3])
                    node = nodes[node_name]
                    print(f'---> Load: {load_name} Point: {node_id} Node: {node.name}')
                    for load in loads.load:
                        label = load.load_name
                        print(f'Load Title: {label} - Nodal Load')
                        node_load[node.name].load = [*load[:6], load[7]]
                    for load in loads.mass:
                        1/0
                        node_load[node.name].mass = [*load[:3], load[7]]
                except IOError: # check if point load on a beam
                    #print(f'Load Title: {label} - Beam Point Load')
                    beams = elements.beams()
                    # TODO: avoid to loop all beams elements
                    for beam_name, beam in beams.items():
                        #Lb = beam.L
                        #p1, p2 = beam.nodes
                        #point_line = DistancePointLine3D(p1[:3], p2[:3])
                        # if point_line.is_on_segment(point[:3]):
                        #
                        if beam.intersect_point(point):
                            p1, p2 = beam.nodes
                            L1 = math.dist(p1[:3], point[:3])
                            L2 = math.dist(p2[:3], point[:3])
                            for load in loads.load:
                                label = load.load_name
                                if math.isclose(L1, 0):
                                    print(f'Load Title: {label} - Nodal Load')
                                    node = nodes[p1.name]
                                    node_load[node.name].load = [*load[:6], load[7]]
                                elif math.isclose(L2, 0):
                                    print(f'Load Title: {label} - Nodal Load')
                                    node = nodes[p2.name]
                                    node_load[node.name].load = [*load[:6], load[7]]
                                else:
                                    print(f'Element: {beam_name} L1: {L1:4.2f} {load[:6]}')
                                    beam_load[beam_name].point = [L1, *load[:6], load[7]]
                            #print('-->')
                            break
                    #
                    #raise Warning('')
        #
        #print('--> end')
    #
#
#
@dataclass
class linefit:
    __slots__ = ['q0', 'q2', 'L', 'L0', 'L1',
                 'L2', 'Lstart', 'Lstop', '_qi']

    def __init__(self, q0:float, q2:float,
                 L:float, L0:float, L1:float) -> None:
        """ """
        self.q0:float = q0
        self.q2:float = q2
        self._qi:float = q0
        #
        self.L:float = L
        self.L0:float = L0
        self.L1:float = L1
        self.L2 = self.L - self.L1
        self.Lstop = self.L - self.L1
    #
    @property
    def slope(self) -> float:
        """ """
        return (self.q2-self.q0)/(self.L2-self.L0)
    #
    def qi(self, x:float) -> list[float]:
        """ """
        q1 = self._qi
        if x > self.L2:
            self._qi = self.q2
            #q2 = round(self.q1 + self.slope * (self.L3-self.L1), 3)
        else:
            self._qi = round(self.q0 + self.slope * (x-self.L0), 3)
        #self._qi = q2
        return [q1, self._qi]
    #
    def Li(self, x:float, Lb:float) -> Union[Exception,List[float]]:
        """ """
        try:
            1/(self.L0 + self.L1)
            if x < self.L0: # no load for this step
                raise RuntimeWarning
            else:
                try:
                    1 / self.Lstop
                    try:
                        Lstart = self.Lstart
                    except AttributeError:
                        Lstart = Lb - (x - self.L0)
                        self.Lstart = 0
                    #
                    if x > self.L2:
                        self.Lstop = 0
                        return [Lstart, x - self.L2]
                    else:
                        return [Lstart, 0]
                except ZeroDivisionError: # no load after this step
                    raise RuntimeWarning
        except ZeroDivisionError:
            return [0, 0]
#
@dataclass
class pointfit:
    __slots__ = ['L', 'L0', 'Lstop']
    
    def __init__(self, L:float, L0:float) -> None:
        """ """
        self.L:float = L
        self.L0:float = L0
        self.Lstop:float = L0
    #
    def Li(self, x:float, Lb:float):
        """ """
        if x < self.L0: # no load for this step
            raise RuntimeWarning
        else:
            try:
                1 / self.Lstop
                self.Lstop = 0
                return Lb - (x - self.L0)
            except ZeroDivisionError:  # no load after this step
                raise RuntimeWarning

#