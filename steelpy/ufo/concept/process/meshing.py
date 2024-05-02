#
# Copyright (c) 2009 steelpy
#

# Python stdlib imports
from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import dataclass
#from typing import Tuple
#import os

# package imports
#
from steelpy.utils.units.main import Units
from steelpy.utils.geometry.L3D import DistancePointLine3D

#
#
class MeshingConcept:
    __slots__ = ["_mesh", "concept"]
    
    def __init__(self, concept) -> None:
        """
        """
        self.concept = concept
        #self._component = component
        #
        #filename = component + "_f2u.db"
        #path = os.path.abspath(filename)
        #db_file = path        
        #
        #self._mesh = Mesh(db_name=component)
        self._mesh = concept._mesh
        #print('--')
    #
    def get_mesh(self) -> None:
        """
        """
        print('-- Meshing Component')
        #self._set_properties()
        self._set_mesh()
        self._set_boundary()
        self._set_load()
        print('-- Meshing Completed')
        return self._mesh
    #
    def _set_properties(self):
        """ """
        # Material
        cmaterial = self.concept.materials()
        dfmat = cmaterial.df
        #
        group = dfmat.groupby("type")
        elastic = group.get_group("elastic")
        elastic = elastic.drop_duplicates(subset=['name'], keep='first')  
        self._mesh._materials._material._elastic.df = elastic
        #self._mesh._materials.df = dfmat
        #
        # Section
        csections = self.concept.sections()
        dfsect = csections.df
        dfsect = dfsect.drop_duplicates(subset=['name'], keep='first')        
        self._mesh._sections.df = dfsect
        #
        #print('-->')
        #1 / 0
    #
    def _set_mesh(self):
        """ """
        print('--- Meshing Concepts')
        mesh = self._mesh
        melements = mesh.element()
        #elem_number = elements.get_number()
        celements = self.concept.element()
        cbeams = celements.beam()
        #
        for key, beam in cbeams.items():
            total_length = beam.L
            p1, p2 = beam.nodes
            node_res = self._get_node_name(p1[:3])
            node_end = self._get_node_name(p2[:3])
            for x, step in enumerate(beam.step):
                idx = x + 1
                mnumber = next(melements.get_number())
                step._mesh = mnumber
                print(f"concept: {key} --> element: {mnumber}")
                try:
                    1/step.length #.value
                    total_length -= step.length
                    coord = beam.find_coordinate(step.length)
                    new_node = self._get_node_name(coord)
                    # elements [node1, node2, material, section, beta, 
                    #           title, idx]
                    melements[mnumber] = ['beam', node_res, new_node,
                                          step.material.name, 
                                          step.section.name, 
                                          beam.beta, key, idx ]
                    node_res = new_node
                except ZeroDivisionError:
                    # elements [node1, node2, material, section, beta, title]
                    melements[mnumber] = ['beam', node_res, node_end,
                                          step.material.name, 
                                          step.section.name, 
                                          beam.beta, key, idx]
            #print('-->')
        print('end meshing')
    #
    def _get_node_name(self, coord):
        """ """
        nodes = self._mesh.node()
        try:
            return nodes.get_point_name(coord)
        except IOError:
            return nodes.get_new_point(coord)
    #
    def _set_boundary(self):
        """ """
        print('--- Meshing Boundary')
        units = Units()
        # Mesh
        mesh = self._mesh
        mnodes = mesh.node()
        melements = mesh.element()
        mbeams = melements.beam()
        mboundaries = mesh.boundary()
        msupports = mboundaries.support()
        # concepts
        cboundary = self.concept.boundary()
        csupports = cboundary.support()
        celements = self.concept.element()
        cbeams = celements.beam()
        #
        # FIXME : check if new loops this works properly
        #
        missing = defaultdict(list)
        # find existing nodes
        for key, value in csupports.items():
            #support = value.support
            support = value.points
            for point in support.points:
                try:
                    node_id = mnodes.get_point_name(point)
                    msupports[node_id] = [*support[:6], key]
                    print(f"Boundary: {key}  @ Node: {node_id}")
                except IOError:
                    missing[key].append(point)
        #
        # if missing boundaries, find if coordinates along members
        #if missing:
        for boundary, points, in missing.items():
            missing_found = [item.name for item in points] #defaultdict(list)
            for key, cbeam in cbeams.items():
                p1, p2 = cbeam.nodes
                #p1, p2 = cbeam.connectivity
                point_line = DistancePointLine3D(p1[:3], p2[:3])
                #for boundary, points, in missing.items():
                #Pb = cboundary[boundary].point
                # Fixme : what format?
                #Pb = [Pb[0].value, Pb[1].value, Pb[2].value]
                for point in points:
                    Pb = point[:3]
                    if point_line.is_on_segment(Pb):
                        missing_found.remove(point.name)
                        left_dist = point_line.left_dist
                        #missing_found[boundary].append(idx)
                        total_length = 0
                        step_no = len(cbeam.step)
                        for step in cbeam.step:
                            beam = mbeams[step._mesh]
                            total_length += beam.L
                            # TODO: capture warning when point misaligned due to tolerances
                            print(f"beam length step: {total_length:1.4e} left: {left_dist:1.4e}")
                            if total_length < left_dist:
                                continue
                            # get node coordinate
                            coord = cbeam.find_coordinate(left_dist) #*units.m
                            new_node = self._get_node_name(coord)
                            # set boundary
                            #support = cboundary[boundary].support
                            #support = csupports[boundary]
                            #if support:
                            msupports[new_node] = boundary
                            print(f"Boundary: {boundary} on Beam: {key} @ Node: {new_node}")
                            #
                            # existing element
                            mnodes = beam.connectivity
                            node_end = mnodes[-1]
                            mnodes[-1] = new_node
                            beam.connectivity = mnodes
                            # new element
                            mnumber = next(melements.get_number())
                            melements[mnumber] = ['beam', new_node, node_end,
                                                  step.material.name,
                                                  step.section.name,
                                                  beam.beta, cbeam.name]
                            # introduce new concept beam step to boundary coord
                            step_no += 1
                            cbeams[key].step[step_no].length = left_dist * units.m
                            cbeams[key].step[step_no].material = step.material.name
                            cbeams[key].step[step_no].section = step.section.name
                            cbeams[key].step[step_no]._mesh = mnumber
                            print(f"concept: {key} --> element: {mnumber}")
                            break
                        #continue
                    if not missing_found:
                        break
                #
                if not missing_found:
                    break
                #
                #for item in reversed(missing_found):
                #    #for item in reversed(items):
                #    #missing.remove(item)
                #    missing[boundary].pop(item)
                # TODO: capture if missing not empty
                #if not missing[boundary]:
                #    #del missing[boundary]
                #    break
        print(' end meshing boundary')
    #
    def _set_load(self):
        """ """
        print( '--- Meshing Basic Load' )
        # Mesh
        mesh = self._mesh
        Mnodes = mesh.node()
        melements = mesh.element()
        Mbeams = melements.beam()
        # Mesh Load
        mload = mesh.load()
        Mlbasic = mload.basic()
        # Concept
        Concept = self.concept
        Celements = Concept.element()
        Cbeams = Celements.beam()
        Cloads = Concept.load()
        Clbasic = Cloads.basic()
        CPoints = Concept.point()
        #
        for Clb_name, Clb_item in Clbasic.items():
            # clone mesh load
            Mlbasic[Clb_name] = Clb_item.title
            mlb_node = Mlbasic[Clb_name].node()
            mlb_beam = Mlbasic[Clb_name].beam()
            #
            # Beam load process
            # TODO : update linefit
            for bcname, CBloads in Clb_item._beam.items():
                cbeam = Cbeams[bcname]
                Lc = cbeam.L #.value
                #
                #print(f'---> Load: {load_name} Beam: {bcname} L: {Lb:4.2f}')
                lcoord_system = CBloads.coordinate_system
                # Beam line load process
                for lbload in CBloads.line[bcname]:
                    #for lbload in rows:
                    label = lbload.load_name
                    print(f'Load Title: {label} - Line Load')
                    waxial = linefit(lbload.qx0, lbload.qx1,
                                     Lc, lbload.L0, lbload.L1)
                    winplane = linefit(lbload.qy0, lbload.qy1,
                                       Lc, lbload.L0, lbload.L1)
                    woutplane = linefit(lbload.qz0, lbload.qz1,
                                        Lc, lbload.L0, lbload.L1)
                    wtorsion = linefit(lbload.qt0, lbload.qt1,
                                        Lc, lbload.L0, lbload.L1)
                    # start loop beam steps
                    xi = 0
                    for step in cbeam.step:
                        bmid = step._mesh
                        beam = Mbeams[bmid]
                        Lbi = beam.L
                        xi += Lbi
                        qaxial = waxial.qi(xi)
                        qinp = winplane.qi(xi)
                        qoutp = woutplane.qi(xi)
                        qtorsion = wtorsion.qi(xi)
                        # check load on segment
                        try:
                            Li = winplane.Li(xi, Lbi)
                        except RuntimeWarning:
                            continue # no load should be applied to this segment
                        # set load for mesh element
                        print(f'Element: {bmid} --> {Lbi:4.2f} {xi:4.2f} {qinp} {qoutp} {Li}')
                        mlbeam = mlb_beam[bmid]
                        mlbeam.coordinate_system = lcoord_system
                        mlbeam.line = [qaxial[0], qinp[0], qoutp[0], qtorsion[0], 
                                       qaxial[1], qinp[1], qoutp[1], qtorsion[1], 
                                       Li[0], Li[1], lbload.title]
                #
                # Beam point load process
                for pbload in CBloads.point[bcname]:
                    #for pbload in rows:
                    label = pbload.load_name
                    print(f'Load Title: {label} - Point Load')
                    L1 = pbload.L0
                    # start loop beam steps
                    xi = 0                    
                    for step in cbeam.step:
                        bmid = step._mesh
                        beam = Mbeams[bmid]
                        Lbs = beam.L
                        xi += Lbs
                        #
                        mlbeam = mlb_beam[bmid]
                        mlbeam.coordinate_system = lcoord_system
                        #
                        if xi < L1: # no load for this beam step
                            continue
                        else:
                            L2 = xi - L1
                            Li = (Lbs - L2)
                            self._beam_pload(beam,
                                             pbload,
                                             mlb_node,
                                             mlbeam,
                                             Li=Li)
                            break
            #
            # Nodal load process
            for CPname, CPloads in Clb_item._node.items():
                point = CPoints[CPname]
                for pload in CPloads.load:
                    label = pload.load_name
                    print(f'Load Title: {label} - Nodal Load')
                    #
                    try:
                        node_name = Mnodes.get_point_name(point[:3])
                        node = Mnodes[node_name]
                        print(f'---> Load: {CPname} Point: {node_name} Node: {node.name}')
                        mlb_node[node.name].load = [*pload[:6], pload[7]]
                    except IOError: # check if point load on a beam
                        for beam_name, beam in Mbeams.items():
                            # TODO: avoid to loop all beams elements
                            if beam.intersect_point(point):
                                self._beam_pload(beam, pload,
                                                 mlb_node, mlb_beam, 
                                                 point=point)
                                break
        #
        Clwave = Cloads.metocean()
        Mlwave = mload.metocean()
        # Metocean
        for name, condition in Clwave.items():
            Mlwave[name] = condition
        #
        print('--> end mesh loading')
    #
    #
    def _beam_pload(self, beam, pload,
                    node_load, mlbeam,
                    Li=None, point=None):
        """ """
        n1, n2 = beam.nodes
        if point:
            L1 = n1.distance(point)
            L2 = n2.distance(point)
        else:
            L1 = Li
            L2 = beam.L - Li
        #
        load = (pload.fx, pload.fy, pload.fz,
                pload.mx, pload.my, pload.mz,
                pload.title)
        #1 / 0 # pload
        # TODO : adjust tolerancea
        if math.isclose(L1, 0, rel_tol=1e-09, abs_tol=0.0):
            print(f'Node {n1.name} Load')
            1 / 0 # FIXME: check load coord_system
            node_load[n1.name].load = load
        elif math.isclose(L2, 0, rel_tol=1e-09, abs_tol=0.0):
            print(f'Node {n2.name} Load')
            1 / 0 # FIXME: check load coord_system
            node_load[n2.name].load = load
        else:
            print(f'Beam {beam.name} Point load L1: {L1:4.2f}')
            #beam_load[beam.name].point = [L1, *load]
            mlbeam.point = [L1, *load]
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
    def Li(self, x:float, Lb:float):
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
#