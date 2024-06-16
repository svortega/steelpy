#
# Copyright (c) 2009 fem2ufo
#
from __future__ import annotations
# Python stdlib imports
#import math
from typing import NamedTuple
from dataclasses import dataclass
#
# package imports
#
import numpy as np
from numpy.matlib import repmat
#import matplotlib.pyplot as plt
#
from steelpy.utils.math.operations import linstep #, linspace, trnsload
from steelpy.utils.dataframe.main import DBframework
#
import xarray as xr
#import matplotlib.pyplot as plt
#
#

#
@dataclass
class BeamMorisonWave:
    __slots__ = ['_beam', 'surface', 'rho', '_up', 
                 '_data', '_type', 'uvector', 'nelev']
    def __init__(self, beam, rho: float,
                 nelev: int, up: str = 'y'):
        """
        rho : : Seawater density (1025)
        """
        self._beam = beam
        self.rho = rho
        self.nelev = nelev
        self._up = up
        # TODO : fix dircos
        dircos = self._beam.dircosines
        if up in ['y']:
            dircos = [dircos[0], dircos[2], dircos[1]]
        self.uvector = UnitVector(*dircos)
    #
    def coordinates(self):
        """get coordinates along beam"""
        n1, n2 = self._beam.nodes
        #sect = self._beam.section.geometry
        #steps = linstep(d=sect.Dh, L=self._beam.L,
        #                steps=self.nelev)
        steps = self.steps()
        coord = [n1[:3]]
        for step in steps[1:]:
            coord.append(self._beam.find_coordinate(node_distance=step))
        coord = list(map(list, zip(*coord)))
        return coord
    #
    def steps(self):
        """beam steps"""
        sect = self._beam.section.geometry
        steps = linstep(d=sect.Dh, L=self._beam.L,
                        steps=self.nelev)
        return steps
    #
    def elevations(self):
        """ Elevation Range"""
        # FIXME: 3D eleveation 
        n1, n2 = self._beam.nodes
        #Tb = self._beam.T3D()
        #nglobal = [*n1[:3], 0, 0, 0,
        #          *n2[:3], 0, 0, 0]
        #nlocal = trnsload(nglobal, r_matrix=Tb)
        #
        n1 = getattr(n1, self._up)
        n2 = getattr(n2, self._up)
        zmax = np.maximum(n1, n2)
        zmin = np.minimum(n1, n2)
        #xx = linspace(start=zmin, stop=zmax, num=nelev+1)
        #yy =  np.linspace(zmin, zmax, nelev+1)
        return np.linspace(zmin, zmax, self.nelev+1)
    #
    def dz(self):
        """ """
        elev = self.elevations(nelev=self.nelev)
        return np.diff(elev)
    #
    def Z(self):
        """ """
        elev = self.elevations(nelev=self.nelev)
        dz = np.diff(elev)
        # locating the middle point of each element
        Z = elev[:-1] + dz
        return Z

    #
    def Dh(self, mg):
        """Diametre hydrodynamic"""
        # TODO: section naming, element direction
        section = self._beam.section.geometry
        Dh = section.Dh
        #mg = self.MG(Z)
        Dh += 2 * mg
        At = np.pi * np.power(Dh, 2) / 4
        return Dh, At
    #
    def local_kin(self, kin, Vc):
        """ Kinematics local to the beam member

        beam : beam element
        kin : kinematics
        Vc : current
        """
        #
        uvector = self.uvector
        # uvector = [0.447, 0.525, 0.724]
        # print(f'Unit Vector [{uvector[0]}, {uvector[1]}, {uvector[2]}]')
        # print('')
        #
        # Components of velocity local to the member
        #
        Un = Vc + kin['ux'] - uvector[0] * (uvector[0] * (Vc + kin['ux']) + uvector[1] * kin['uz'])
        Vn = kin['uz'] - uvector[1] * (uvector[0] * (Vc + kin['ux']) + uvector[1] * kin['uz'])
        Wn = - uvector[0] * (uvector[0] * (Vc + kin['ux']) + uvector[1] * kin['uz'])
        #
        # print('')
        # print('========================================')
        # print('Components of velocity local to the member [N/m]')
        # print(f'Un ={np.max(Un): 1.4e}, Vn={np.max(Vn): 1.4e}, Wn={np.max(Wn): 1.4e}')
        # print(f'Un ={np.min(Un): 1.4e}, Vn={np.min(Vn): 1.4e}, Wn={np.min(Wn): 1.4e}')
        #
        # Water velocity normal to the cylinder axis
        #
        vn = Vc + np.sqrt(np.power(kin['ux'], 2) + np.power(kin['uz'], 2)
                          - np.power(uvector[0] * kin['ux'] + uvector[1] * kin['uz'], 2))
        #
        # print('')
        # print('Water velocity normal to the cilinder axis [N/m]')
        # print(f'vnmax ={np.max(vn): 1.4e}, vnmin={np.min(vn): 1.4e}')
        #
        # Components of acceleration local to the member
        #
        Anx = kin['ax'] - uvector[0] * (uvector[0] * kin['ax'] + uvector[1] * kin['az'])
        Anz = kin['az'] - uvector[1] * (uvector[0] * kin['ax'] + uvector[1] * kin['az'])
        Any = - uvector[0] * (uvector[0] * kin['ax'] + uvector[1] * kin['az'])
        #
        #
        # print('')
        # print('Components of acceleration local to the member [N/m]')
        # print(f'Anx ={np.max(Anx): 1.4e}, Any={np.max(Anz): 1.4e}, Anz={np.max(Any): 1.4e}')
        # print(f'Anx ={np.min(Anx): 1.4e}, Any={np.min(Anz): 1.4e}, Anz={np.min(Any): 1.4e}')
        # print('========================================')
        #
        return KinVel(Un, Vn, Wn, vn, self.rho), KinAcc(Anx, Anz, Any, self.rho)
    #
    def vn(self, kin, Vc):
        """ Absolute Water velocity normal to the cylinder axis
        
        Vc : current velocity
        """
        uvector = self.uvector
        vn = Vc + np.sqrt(np.power(kin['ux'], 2) + np.power(kin['uz'], 2)
                          - np.power(uvector.x * kin['ux'] + uvector.z * kin['uz'], 2))
        return vn
    #
    def Un(self, kin, Vc):
        """
        Instantaneous undisturbed velocity resolved normal
        to the member including both wave and current
        
        kin :
        Vc : current velocity
        """
        uvector = self.uvector
        #
        shape = kin['ax'].shape
        #Vc = permute1(Vc, order=shape[0])
        Vc = permute3(Vc, order=shape[0])
        #
        # Components of velocity local to the member
        comp0 = uvector.x * (Vc + kin['ux']) + uvector.z * kin['uz']
        #
        Un = Vc + kin['ux'] - uvector.x * comp0    # x
        Wn = - uvector.y * comp0                   # y
        Vn = kin['uz'] - uvector.z * comp0         # z
        # Water velocity normal to the cylinder axis
        vn = self.vn(kin, Vc)
        #
        return KinVel(Un, Vn, Wn, vn, self.rho)
    #
    def An(self, kin):
        """
        Instantaneous undisturbed acceleration resolved normal
        to the member
        
        kin : 
        """
        uvector = self.uvector
        comp0 = (uvector.x * kin['ax'] + uvector.z * kin['az'])
        # components of acceleration normal to the member in the x,y and z directions
        Anx = kin['ax'] - uvector.x * comp0
        Any = - uvector.y * comp0
        Anz = kin['az'] - uvector.z * comp0
        return KinAcc(Anx, Anz, Any, self.rho)
    #
    def dF(self, Dh, At, Cd, Cm,
           kinvel, kinacc,
           time, eta):
        """Components of the force per unit of cilinder length acting in
        the x, y and z dir are given by the generalized Morisson equation
        
        dF = Fn + Ft
        
        Dh : Hydrodynamic diametre
        At : Cross sectional Area
        Cd : Drag coefficient
        Cm : Mass coefficient
        
        """
        #
        dmx, dmy, dmz = kinacc.FIn(At, Cm)
        #
        ddx, ddy, ddz = kinvel.FDn(Dh, Cd)
        #
        fx = dmx + ddx
        fy = dmy + ddy
        fz = dmz + ddz
        #
        Fi = xr.Dataset(data_vars={'fx': fx,'fy': fy,'fz': fz})
        #
        #
        elev = self.elevations()
        coord = self.coordinates()
        steps = self.steps()        
        return BeamUnitForce(Fi, time, eta,
                             coord, steps, elev, 
                             self._beam, self._up)        
    #
    #
    def Fwave(self, Vc, MG, Cd, Cm,
              kinematics, eta: list, 
              time: list):
        """
        Wave force on a slender cilindrical element
        
        Vc : Current velocity
        MG : Marine Growth
        Cd : Drag Coefficient
        Cm : Inertia Coefficient
        WKF : Wave Kinematic Factor
        kinematics : Kinematic class
        nelev : number of elevations
        """
        #elev = self.elevations()
        # FIXME : element direction
        Dh, At = self.Dh(mg=MG)
        #dz = self.dz(nelev=nelev)
        #dz = np.diff(elev)
        #
        #shape = kinematics['ax'].shape
        #Vc = permute1(Vc, order=shape[0])
        #
        #coord = self.coordinates()
        #kinacc = {}
        #kinvel = {}
        #for step in range(len(wave_phase)):
        #    phase = wave_phase[step]
        #    print(f'----> step {step}, eta {eta[step]}, phase {phase}')
        #    item = kinematics.roll(x=step)
        #    #
        #    kinstep = item.interp(x=np.array(coord[0]),
        #                          y=np.array(coord[2]),
        #                          z=np.array(coord[1]),
        #                          method="linear",
        #                          assume_sorted=True,
        #                          kwargs={"fill_value": 0})
        #    #
        #    kinacc[phase] = self.An(kinstep)
        #    kinvel[phase] = self.Un(kinstep, Vc[step])
        #    #
        #
        #
        kinacc = self.An(kinematics)
        kinvel = self.Un(kinematics, Vc)
        #
        return self.dF(Dh, At, Cd, Cm,
                       kinvel, kinacc,
                       time=time, eta=eta)
        #return udl
    #
    #
#
#
class UnitVector(NamedTuple):
    """Unit vector """
    x: float
    y: float
    z: float
#
#
class BeamUnitForce(NamedTuple):
    """Components of the force per unit of cilinder lenght"""
    Fi: list
    time: list
    eta: list
    coordinates: list
    steps: list
    elevation: list
    beam: tuple
    up: str
    #
    #
    @property
    def df(self):
        """Dataframe of beam's partial linearly variable load
        
        [load_title, 'beam', beam_name, 'line',  qx0,qy0,qz0, qx1,qy1,qz1, L0,L1, comment(optional)]"""
        #
        #coords = self.qx.coords
        #rows = coords['x'].values
        #cols = coords['z'].values
        #wlength = coords['length'].values
        #
        #
        #Fx, Fy, OTM = self.span_loading()
        #
        #
        # FIXME: wave system to beam local system (is this fixed already?)
        #
        #qitem = self.qx.to_dataframe(name='qx').reset_index()
        #qy = self._get_line(qname='qx', qitem=qitem)
        #qitem = self.qy.to_dataframe(name='qy').reset_index()
        #qz = self._get_line(qname='qy', qitem=qitem)
        #qitem = self.qz.to_dataframe(name='qz').reset_index()
        #qx = self._get_line(qname='qz', qitem=qitem)
        #
        # TODO: L1 and L2 should not be zero for all cases
        #dftemp = []
        #for x, row in enumerate(rows):
        #    for idx, wstep in enumerate(wlength):
        #        for hstep, col in enumerate(cols):
        #            ldata = list(zip(qx[idx][hstep], qy[idx][hstep], qz[idx][hstep]))
        #            dftemp.append(['beam', self.beam_name, self.beam_number, 'line',
        #                           *ldata[0], *ldata[1], 0, 0,
        #                           float(Fx[x, hstep, idx].values),
        #                           float(OTM[x, hstep, idx].values), 
        #                           row, wstep, col])
        #
        dftemp = self.solve()
        #
        # setup df's columns
        header = ['element_type', 'element_name', 'element_id', 
                  'load_type',
                  'qx0', 'qy0', 'qz0', 'qx1', 'qy1', 'qz1',
                  'L0', 'L1', 'BS', 'OTM', 
                  'x', 'y', 'z']
        #
        df = DBframework()
        dfload = df.DataFrame(data=dftemp, columns=header, index=None)
        # print('--->')
        # 1 / 0
        return dfload
    #
    def solve2(self):
        """ """
        coords = self.qx.coords
        rows = coords['x'].values
        cols = coords['z'].values
        wlength = coords['length'].values
        #
        #
        Fx, Fy, OTM = self.span_loading()
        #
        #
        # FIXME: wave system to beam local system (is this fixed already?)
        #
        qitem = self.qx.to_dataframe(name='qx').reset_index()
        qy = self._get_line(qname='qx', qitem=qitem)
        qitem = self.qy.to_dataframe(name='qy').reset_index()
        qz = self._get_line(qname='qy', qitem=qitem)
        qitem = self.qz.to_dataframe(name='qz').reset_index()
        qx = self._get_line(qname='qz', qitem=qitem)
        
        # TODO : fix torsion
        qitem['qt'] = float(0.0)
        qt = self._get_line(qname='qt', qitem=qitem)
        #
        dftemp = []
        for x, row in enumerate(rows):
            for idx, wstep in enumerate(wlength):
                for hstep in range(len(qx[idx])):
                #for hstep, col in enumerate(cols):
                    ldata = list(zip(qx[idx][hstep],
                                     qy[idx][hstep],
                                     qz[idx][hstep],
                                     qt[idx][hstep]))
                    #
                    dftemp.append(['beam', self.beam.name, self.beam.number, 'line',
                                   *ldata[0], *ldata[1],                   # q0, q1
                                   ldata[2][0], self.beam.L - ldata[3][0], # L0, L1
                                   float(Fx[x, hstep, idx].values),
                                   float(OTM[x, hstep, idx].values), 
                                   row, wstep, cols[hstep]])
        return dftemp
    #
    def solve(self):
        """ """
        elev = self.elevation
        # TODO : beam direction
        coord = self.coordinates
        if self.up in ['y']:
            coord = [coord[0], coord[2], coord[1]]
        coord = UnitVector(*coord)        
        idx = [x for x in range(len(coord[0]))]
        bsteps = np.abs(np.concatenate(([0], np.diff(coord.z))))       
        #
        dftemp = []
        for step in range(len(self.eta)):
            #phase = self.wave_phase[step]
            #print(f'----> step {step}, eta {self.eta[step]}, phase {phase}')
            item = self.Fi.roll(x=step)
            Fstep = item.interp(x=np.array(coord.x),
                                y=np.array(coord.y),
                                z=np.array(coord.z),
                                method="linear",
                                assume_sorted=True,
                                kwargs={"fill_value": 0})
            #
            # FIXME : magic step to be fixed
            fy = Fstep['fx'].data[idx, idx, idx]
            fz = Fstep['fy'].data[idx, idx, idx]
            fx = Fstep['fz'].data[idx, idx, idx]
            # TODO : fix torsion
            ft = fz * 0
            #
            BS = np.sqrt(np.power((fx * bsteps), 2)
                         + np.power((fy * bsteps), 2))
            #
            OTM = BS * elev
            #
            dftemp.append(np.array([fx, fy, fz, ft, BS, OTM]))
            #
            #
            #for idx, x, in enumerate(coord[0]):
            #    y, z = coord[1][idx], coord[2][idx]
            #    print(x, y, z)
            #    ldata = item.interp(x=x, y=z, z=y, 
            #                        method='linear',
            #                        assume_sorted=True, 
            #                        kwargs={'fill_value': 0,
            #                                'bounds_error': True})
            #    #
            #    print(ldata)
            #    #
            #    #ldata = ldata.to_array().to_numpy()
            #    #print(ldata)
            #    #dftemp.append(['beam', self.beam.name, self.beam.number, 'line',
            #    #               *ldata[0], *ldata[1],                   # q0, q1
            #    #               ldata[2][0], self.beam.L - ldata[3][0]]) # L0, L1
            #    #               #float(Fx[x, hstep, idx].values),
            #    #               #float(OTM[x, hstep, idx].values), 
            #    #               #row, wstep, cols[hstep]])
        #        
        #
        # distance along beam 
        #LbeamX = np.array(coord[1])
        Lbeam = self.steps
        bl = self.beam.L
        #
        n1, n2 = self.beam.nodes
        n1 = getattr(n1, self.up)
        #
        if n1 != elev[0]:
            Lbeam = Lbeam[::-1]
        #
        #
        df = DBframework()
        qload = []
        for idx, step in enumerate(dftemp):
            q = step.T
            qtemp = [[self.eta[idx], self.time[idx], 
                      *q[x-1], *q[x], Lbeam[x-1], bl - Lbeam[x]]
                     for x in range(1, len(q))]
            #for x in range(1, len(q)):
            #    qtemp.append([*q[x-1], *q[x], Lbeam[x-1], bl - Lbeam[x]])              
            #qload.append(qtemp)
            #
            qload.append(df.DataFrame(data=qtemp,
                                      columns=['eta', 'time', 
                                               'qx0', 'qy0', 'qz0', 'qt0', 'BS0', 'OTM0',
                                               'qx1', 'qy1', 'qz1', 'qt1', 'BS1', 'OTM1',
                                               'L0', 'L1'],
                                      index=None))
        #
        #
        qload = df.concat(qload)
        qload['BS'] = qload['BS0'] + qload['BS1']
        qload['OTM'] = qload['OTM0'] + qload['OTM1']
        qload['element_type'] = 'beam'
        qload['element_name'] = self.beam.name
        qload['element_id'] = self.beam.number
        qload['type'] = 'line'
        #
        qload.drop(columns=['BS0', 'BS1', 'OTM0', 'OTM1'],
                   inplace=True)
        #
        return qload
    #
    def _get_line2(self, qname: str, qitem):
        """ """
        elev = self.elevation
        coords = qitem.coords
        rows = coords['x'].values
        cols = coords['col'].values
        wlength = coords['length'].values
        #
        1 / 0
        grpx = qitem.groupby('x')
        for keyx, item in grpx:
            grp2 = item.groupby('length')
            for keyw, wstep in grp2:
                keyw, wstep
                if -93.75 in wstep.coords['col']:
                    print('--')

        for row in rows:
            for wstep in wlength:
                for el in elev:
                    try:
                        cols = qitem.sel(x=row, col=el, length=wstep)
                        print('--')
                    except TypeError:
                        step = 0
        #
        1 / 0
        print('--')

    #
    #
    def _get_line(self, qname: str, qitem):
        """ """
        # FIXME : this is bad code
        elev = self.elevation
        dz = np.diff(elev)
        Lbeam = [float(item) * i for i, item in enumerate(dz)] + [self.beam.L]
        n1, n2 = self.beam.nodes
        n1 = getattr(n1, self.up)
        #n2 = getattr(n2, self.up)        
        #zmin = np.minimum(n1, n2)
        rev = False
        if n1 != elev[0]:
            rev = True
            Lbeam = Lbeam[::-1]
        #
        # data = qx.groupby(['x'])[name].agg(lambda x : x.tolist())
        qgrp = qitem.groupby(['x', 'length'])[['z', qname]]
        # TODO : optimize
        load_1 = []
        for key, item in qgrp:
            z, q = item.to_numpy().T
            #if rev:
            #    Lbeam = Lbeam[::-1]
            #    q = q[::-1]
            # changing z by elev
            qload = []
            for x in range(1, z.size):
                #try:
                #1 / float(q[x])
                if rev:
                    qload.append([q[x-1], q[x], Lbeam[x], Lbeam[x-1]])
                else:
                    qload.append([q[x-1], q[x], Lbeam[x-1], Lbeam[x]])
                #except ZeroDivisionError:
                #    #print('here')
                #    #pass
                #    if rev:
                #        qload.append([0.0, 0.0, Lbeam[x], Lbeam[x-1]])
                #    else:
                #        qload.append([0.0, 0.0, Lbeam[x-1], Lbeam[x]])
            #load_2 = []
            #step = None
            #for el in elev:
            #    idx = item.index[item['z'] == el].tolist()
            #    try:
            #        qn = np.round(item[qname][idx[0]], decimals=3)
            #        load_2.append([step, qn])
            #        step = qn
            #    except IndexError:
            #        step = 0
            # print(key, item)
            load_1.append(qload)
        #
        if n1 != elev[0]:
            load_1 = [item[::-1] for item in load_1]
        #1 / 0
        return load_1

    #
    def span_loading2(self):
        """ """
        qx = self.qx
        qz = self.qz
        dz = permute2(self.elevation, (qx.shape[0], qx.shape[2]), 1)
        Fx = qx.cumsum(dim='z') * dz
        Fz = qz.cumsum(dim='z') * dz
        #
        #Z = self.elevation[:-1] + self.dz
        Z = self.elevation
        Z = permute2(Z, (qx.shape[0], qx.shape[2]), 1)
        OTM = Fx * Z
        return Fx, Fz, OTM
    #
    def span_loading(self):
        """ """
        qx = self.Fi['fx']
        qz = self.Fi['fz']
        #
        blen = np.diff(self.coordinates[1])
        dz = permute5(blen, (qx.shape[0], qx.shape[1]), 1)
        Fx = qx.cumsum(dim='z') * dz
        Fz = qz.cumsum(dim='z') * dz
        #
        #Z = self.elevation[:-1] + self.dz
        Z = self.elevation
        Z = permute2(Z, (qx.shape[0], qx.shape[2]), 1)
        OTM = Fx * Z
        return Fx, Fz, OTM
#
#
class KinVel(NamedTuple):
    """

    Un : Kinematic components of velocity x
    Vn : Kinematic components of velocity y
    Wn : Kinematic components of velocity z
    vn : Fluid velocity normal to the cylinder axis
    rho : Sea water density (1025)
    """
    Un: list
    Vn: list
    Wn: list
    vn : list
    rho: float
    #
    def fdn(self, D, cd, UX, vn):
        """
        D  : Member diametre
        cd : Drag coefficient
        Ux : Instantaneus velocity resolved normal to the member
        Vn : Fluid velocity normal to the cylinder axis
        """
        # drag load per unit length
        Fdn = 0.5 * self.rho * cd * D * UX * vn
        return Fdn
    #
    def FDn(self, Dt:float, Cd:float):
        """
        Component of drag force per unit of cylinder length

        Dt : Diametre tubular
        Cd : Drag coefficient
        
        Return:
        FDn [x,y,z]
        """
        Dh = permute5(Dt, (self.Un.shape[0], self.Un.shape[1]), 1)
        cd = permute5(Cd, (self.Un.shape[0], self.Un.shape[1]), 1)
        #
        FDnx = self.fdn(Dh, cd, self.Un, self.vn)
        FDny = self.fdn(Dh, cd, self.Vn, self.vn)
        FDnz = self.fdn(Dh, cd, self.Wn, self.vn)
        #
        return FDnx, FDny, FDnz
#
#
class KinAcc(NamedTuple):
    """
    Anx : Kinematic components of velocity x
    Any : Kinematic components of velocity y
    Anz : Kinematic components of velocity z
    rho : Sea water density (1025)
    """
    Anx: list
    Any: list
    Anz: list
    rho: float
    #
    def fin(self, at, cm, An):
        """
        at : Cross sectional area
        cm : Inertia coeffient
        An : Instantaneus acceleration resolved normal to the member
        
        Retuns:
        Fin : inertia load per unit length
        """
        # inertia load per unit length
        Fin = self.rho * cm * at * An
        return Fin
    #
    def FIn(self, At:float, Cm:float):
        """
        Component of inertia force per unit of cylinder length normal to the member

        At : Area tubular
        Cm : Mass coefficient
        
        Returns
        FIn [x,y,z]
        """
        at = permute5(At, (self.Anx.shape[0], self.Anx.shape[1]), 1)
        cm = permute5(Cm, (self.Anx.shape[0], self.Anx.shape[1]), 1)
        #
        FInx = self.fin(at, cm, self.Anx)
        FIny = self.fin(at, cm, self.Any)
        FInz = self.fin(at, cm, self.Anz)
        return FInx, FIny, FInz
#
#
def permute5(A, order, axis:int=1):
    """ """
    A1 = np.transpose(A)
    A1 = repmat(A1, order[1], axis)
    A1 = np.expand_dims(A1, axis=0)
    A1 = np.transpose(A1)
    A1 = np.tile(A1, order[0])
    A1 = np.transpose(A1)
    return A1
#
def permute2(A, order, axis:int=1):
    """ """
    A1 = repmat(A, order[1], axis)
    A1 = np.transpose(A1)
    A1 = np.expand_dims(A1, axis=0)
    A1 = np.transpose(A1)
    A1 = np.tile(A1, order[0])
    return np.transpose(A1)
#
def permute3(A, order):
    """ """
    A1 = np.expand_dims(A, axis=0)
    A1 = np.transpose(A1)
    A1 = np.tile(A1, order)
    A1 = np.transpose(A1)
    return A1
#
def permute1(A, order):
    """ """
    A1 = np.transpose(A)
    A1 = np.expand_dims(A1, axis=0)
    A1 = np.transpose(A1)
    A1 = np.tile(A1, order)
    A1 = np.transpose(A1)
    return A1
#