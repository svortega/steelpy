#
# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
from array import array
from dataclasses import dataclass
#from typing import NamedTuple
#import math

# package imports
from steelpy.utils.math.operations import linstep
#
import matplotlib.pyplot as plt
import numpy as np
from numpy.matlib import repmat
import xarray as xr
from steelpy.utils.dataframe.main import DBframework
#
#

#
#
@dataclass
class KinematicResults:
    __slots__ = ['depth_points', 'surface', #'_theta',
                 '_data', '_type', '_stickup']
    
    def __init__(self, surface, kindata, 
                 #depth_points:int,
                 #theta:float,
                 stickup:float=0.10):
        """
        theta : wave angle in degrees
        """
        self.surface = surface
        #self.depth_points = depth_points
        #
        self._data = kindata
        self._type = 'regular'
        self._stickup = stickup # percentage
        #self._theta = theta
    #
    #
    def get_data(self, name: str, title: str):
        """ """
        #
        rows =  np.array([0.0])
        cols = self.z
        #cols =  np.arange(self.depth_points)
        #
        #etas = self.surface.eta
        xx =  self.surface.x
        item = self._data.groupby(['length'])[['elevation', name]]
        data = self.kindf(item, xx, zdepth=cols)
        #
        #
        ##lambda_x = self._data['x'] #.iloc[:self.surface_points]
        #data = self._data.groupby(['x'])[name].agg(lambda x : x.tolist())
        #lambda_x = data.index.values
        ##data = np.array(data.values)
        ##
        #data = self._data[name].values
        ## [lenght, z]
        #surface_points = len(self.surface.x)
        ##depth_points = len(self.depth_points)
        #data = data.reshape((surface_points, self.depth_points))
        # [z, lenght]
        data = np.transpose(data)
        # insert x axis to simulate a 3d wave [x, z, lenght]
        data = np.expand_dims(data, axis=0)
        ##
        return self.to_xarray(rows=rows,
                              cols=cols,
                              length=xx,
                              data=data,
                              name=title)
        #
        #return new_data
    #
    @property
    def Hw(self):
        """ """
        return self.surface.Hw
    #
    @property
    def Tw(self):
        """ """
        return self.surface.Tw    
    #
    @property
    def d(self):
        """ """
        return self.surface.d   
    #
    @property
    def finite_depth(self):
        """ """
        return self.surface.finite_depth     
    #
    #
    def __setitem__(self, name:str, parameters) -> None:
        """
        parameters = []
        """
        self._data[name] = parameters
    
    def __getitem__(self, name:str):
        """
        """
        return self._data[name]
    #
    def __str__(self):
        """ """
        output = "\n"
        output += "#\n"
        output += "# Velocity and acceleration profiles and Bernoulli checks\n"
        if self.finite_depth:
            header = "X"
            output += "#*******************************************************************************\n"
            output += "# y        u       v    dphi/dt   du/dt   dv/dt  du/dx   du/dy Bernoulli check\n"          
        else:
            header = "kX"
            output += "#*******************************************************************************\n"
            output += "# ky       u       v    dphi/dt   du/dt   dv/dt  du/dx   du/dy Bernoulli check\n"
        #
        output += "#*******************************************************************************\n"
        output += "# Note that increasing X and Phase here describes half of a wave for X >= 0.\n"
        #output += "# In a physical problem, where we might have a constant x, because\n"
        #output += "# the phase X = x - c * t, then as time t increases, X becomes increasingly\n"
        #output += "# negative and not positive as passing down the page here implies.\n"
        output += "#\n"       
        #
        data = self._data[['x', 'phase', 'z', 'u', 'v', 'dphidt', 'ut', 'vt', 'ux', 'uz', 'Bernoulli_check']]
        grpkin = data.groupby(['x', 'phase'])
        for key, items in grpkin:
            output += "# {:} = {: 8.4f}, Phase = {: 6.1f}\n".format(header, *key)
            
            for row in items.itertuples():
                #test = [f'{step: 1.3e}' for step in row[3:]]
                #output += (*test, sep=' ')
                for step in row[3:]:
                    output +=  f'{step: 1.3e} '
                output += "\n"
        #
        return output    
    #
    #
    #@property
    #def eta(self):
    #    """surface"""
    #    return self.surface
    #
    # -------------------------------------------------------
    #
    @property
    def z(self):
        """ 
        z  : finite water depth
        kz : infinite water depth
        """
        #
        grpz = self._data.groupby('length')
        gb_groups = grpz.groups
        key_list = list(gb_groups.keys())
        points = len(gb_groups[key_list[0]])
        #
        eta = self.surface.eta
        #points = self.depth_points
        stickup = self._stickup
        #
        crestmax = np.ceil(eta.max()) * (1 + stickup)
        crestmin = np.floor(eta.min())    
        step1 = int(np.ceil(points / 2))
        step2 = int(points - step1)
    
        return np.hstack([np.linspace(-self.d, crestmin, step1, endpoint=False),
                          np.linspace(crestmin, 0, step2, endpoint=False),
                          np.linspace(0, crestmax, step2, endpoint=False),
                          np.linspace(crestmax, 2*crestmax, step2)])        
    #    
    #
    @property
    def ux(self):
        """
        velocity horizontal
        
        Return:
        [Wave surface horizoltal coord, time,  Water depth]
        """
        return self.get_data('u', 'ux')
    #
    @property
    def uz(self):
        """ velocity vertical
        
        Return:
        [Wave surface vertical coord, time,  Water depth]
        """
        return self.get_data('v', 'uz')
    #
    @property
    def ax(self):
        """ acceleration horizontal"""
        return self.get_data('ut', 'ax')
    #
    @property
    def az(self):
        """ acceleration vertical """
        return self.get_data('vt', 'az')
    #
    @property
    def pressure(self):
        """ """
        return self.get_data('pressure', 'pressure')
    #
    @property
    def phi(self):
        """Wave potential function"""
        #return self._data[['phase','', 'dphidt']]
        return self.get_data('dphidt', 'phi')
    #
    # -------------------------------------------------------
    #
    def plot(self):
        """ """
        new_data = self._data.groupby(['length'])
        #xxx = self.kout
        try:
            zlev = new_data['elevation'].agg(lambda x : x.tolist())
        except:
            zlev = new_data['kz'].agg(lambda x : x.tolist())
        #
        velh = new_data['u'].agg(lambda x : x.tolist())
        velv = new_data['v'].agg(lambda x : x.tolist())
        acch = new_data['ut'].agg(lambda x : x.tolist())
        accv = new_data['vt'].agg(lambda x : x.tolist())
        #velh = self.ux
        #velv = self.uz
        #acch = self.ax
        #accv = self.az
        fig, axs = plt.subplots(2)
        
        for key, z in zlev.items():
            axs[0].plot(velh[key], z, color='blue', linewidth=1.0)
            axs[0].plot(velv[key], z, color='red', linewidth=1.0)
            #
            axs[1].plot(acch[key], z, color='blue', linewidth=1.0)
            axs[1].plot(accv[key], z, color='red', linewidth=1.0)            
        #axs[0].set_title('Velocities profiles over half a wave')
        axs[0].set_xlabel('Velocities $u$ and $v$ (m/sec)', fontsize = 8)
        axs[0].set_ylabel('$z$ (m)', fontsize = 8)
        #
        #axs[1].set_title('Acceleration profiles over half a wave')
        axs[1].set_xlabel('Accelerations $u$ and $v$ (m/sec^2)', fontsize = 8)
        axs[1].set_ylabel('$z$ (m)', fontsize = 8)        
        #
        # common axis labels
        fig.suptitle('Velocities and Acceleration profiles over half a wave',
                     fontsize = 10)
        #fig.supxlabel('fig.supxlabel')
        #fig.supylabel('$z$ (m)', fontsize = 10)        
        plt.show()
        #print('--')
    #
    def plot_vectorfield(self):
        """ """
        #
        new_data = self._data.groupby(['length'])
        try:
            zlev = new_data['elevation'].agg(lambda x : x.tolist())
        except:
            zlev = new_data['kz'].agg(lambda x : x.tolist())
        #
        velh = new_data['u'].agg(lambda x : x.tolist())
        velv = new_data['v'].agg(lambda x : x.tolist())        
        #velh = self.ux
        #velv = self.uz        
        #xx = self.x
        for key, z in zlev.items():
            Y = z
            X = [key for x in z]
            U = velh[key]
            V = velv[key]
            plt.quiver(X, Y, U, V, color='b',
                       scale_units='xy', scale=1, pivot='mid')
        # Vector origin location
        #X = [0]
        #Y = [0]
        #
        # Directional vectors
        #U = [2] 
        #V = [1] 
        
        # Creating plot
        #plt.quiver(X, Y, U, V, color='b', units='xy', scale=1)
        plt.title('Single Vector')
        
        # x-lim and y-lim
        #plt.xlim(-2, 5)
        #plt.ylim(-2, 2.5)
        
        # Show plot with grid
        plt.grid()
        plt.show()
        #print('---')
    #
    #
    def plot_contour(self):
        """ """
        #from matplotlib import rcParams
        #new_data = self._data.groupby(['length'])
        #try:
        #    zlev = new_data['elevation'].agg(lambda x : x.tolist())
        #except:
        #    zlev = new_data['kz'].agg(lambda x : x.tolist())
        #
        #lenght = new_data['length'].agg(lambda x : x.tolist())
        #velh = new_data['u'].agg(lambda x : x.tolist())
        #velv = new_data['v'].agg(lambda x : x.tolist())
        #
        #
        #contour_data =  self._data[['length', 'v', 'elevation']]
        #contour_data = contour_data.rename(columns={'length': 'x', 'v': 'y', 'elevation': 'z'})
        #Z = contour_data.pivot_table(index='x', columns='y', values='z').T.values
        #
        #X_unique = np.sort(contour_data.x.unique())
        #Y_unique = np.sort(contour_data.y.unique())
        #X, Y = np.meshgrid(X_unique, Y_unique)
        #
        #fig, ax = plt.subplots(1, 1)
        #
        # plots filled contour plot 
        #ax.contourf(lenght.iloc[0], velh.iloc[0], zlev.iloc[0])
        #ax.contourf(self._data['length'], self._data['v'], self._data['elevation'])
        #ax.contour(X, Y, Z)
        #
        #ax.set_title('Filled Contour Plot') 
        #ax.set_xlabel('feature_x') 
        #ax.set_ylabel('feature_y')
        #
        #
        vx = self.ux
        #vx.plot.contour()
        #
        # Initialize plot objects
        #rcParams['figure.figsize'] = 5, 5 # sets plot size
        #fig = plt.figure()
        #ax = fig.add_subplot(111)
        
        # Generate a contour plot
        #cp = ax.contour(X, Y, Z)        
        vx.plot()
        plt.show()         
        #
        #1 / 0
        
    #
    #
    # -------------------------------------------------------
    #
    def get_kin(self, elev: list, krf: list):
        """
        elev : Elevations
        krf : Kinematic reduction factor of horizontal velocity and acceleration
        """
        rows =  np.array([0.0])
        xx =  self.surface.x
        #items = self._data.groupby(['x'])[['z', 'u', 'ut']]
        items =  ['u', 'ut', 'v', 'vt']
        title = ['ux', 'ax', 'uz', 'az']
        #
        kdf = dict()
        for i, name in enumerate(items):
            item = self._data.groupby(['x'])[['z', name]]
            data = self.kindf(item, xx, elev)
            # krf on horizontal vel and acc
            if name in ['u', 'ut']:
                data *= krf     
            data = np.transpose(data)
            # insert x axis to simulate a 3d wave [x, z, lenght]
            data = np.expand_dims(data, axis=0)
            newname = title[i]
            kdf[newname] = self.to_xarray(rows=rows,
                                          cols=elev,
                                          length=xx,
                                          data=data,
                                          name=newname)
        #print('-->')
        #return kdf
        return xr.Dataset(data_vars=kdf)
    #
    def kindf(self, items, xx, zdepth):
        """ """
        dfkin = np.zeros((xx.size, zdepth.size))
        for i, step in enumerate(xx):
            #new = []
            item = items.get_group(name=(step, ))
            for j, point in enumerate(zdepth):
                dfkin[i, j] = np.interp(point,
                                        item.iloc[:, 0],
                                        item.iloc[:, 1],
                                        right=0)
            #
            #dfkin.append(new)
        #
        return dfkin
    #    
    def to_xarray(self, rows, cols, length, data, name: str):
        """ """
        #row_meshgrid, col_meshgrid = np.meshgrid(rows, cols, indexing='ij')
        return xr.DataArray(data=data,
                            dims=['x', 'z', 'length'],
                            coords=[rows, cols, length], 
                            #coords={
                            #    'row': (["x","z"], row_meshgrid),
                            #    'col': (["x","z"], col_meshgrid),
                            #    #'x': rows,
                            #    #'z': cols,                          
                            #    'length' : length},
                            name=name)
        #return ds
    #
    # -------------------------------------------------------
    #     
    #
    def get_kin2(self, beam, nelev: int):
        """ """
        xx =  self.surface.x
        #
        n1, n2 = beam.nodes
        sect = beam.section.geometry
        steps = linstep(d=sect.Dh, L=beam.L, steps=nelev)
        coord = [n1[:3]]
        for step in steps[1:]:
            coord.append(beam.find_coordinate(node_distance=step))
        #
        coord = list(map(list, zip(*coord)))
        #
        df_mulidx = self._data[['phase', 'x', 'z', 'u', 'ut', 'v', 'vt']].copy()
        df_mulidx.insert(1, 'y', float(0.0))
        #df_mulidx.loc[:, 'y'] = float(0.0)
        df_mulidx.rename(columns={'u': 'ux', 'ut': 'ax',
                                  'v': 'uz', 'vt': 'az'},
                         inplace=True)
        #
        dataset = {}
        dfgrp = df_mulidx.groupby('phase') #,'y','z'
        for key, item in dfgrp:
            #print(key)
            data = item.drop(columns=['phase'])
            data.set_index(['x', 'y', 'z'], inplace=True)
            #dataset[key] = data.to_xarray()
            data = data.to_xarray()
        #
        #for key, item in dataset.items():
            #print(f'phase : {key}')
            # TODO: coord in beam local system
            dataset[key] = data.interp(x=xr.DataArray(coord[0], dims="coordinates"),
                                       y=xr.DataArray(coord[2], dims="coordinates"),
                                       z=xr.DataArray(coord[1], dims="coordinates"),
                                       method="linear",
                                       kwargs={"fill_value": 0})
            #print(dataset[key])
        #
        #1 / 0
        return dataset
    #
    def get_kin3(self, elev: list, krf: list):
        """ """
        #df_mulidx = self._data[['phase', 'x', 'z', 'u', 'ut', 'v', 'vt']].copy()
        #
        items =  ['u', 'ut', 'v', 'vt']
        title = ['ux', 'ax', 'uz', 'az']
        #
        rows = np.array([0.0])
        #cols = self._data['x']
        #
        kdf = dict()
        for i, name in enumerate(items):
            dgroup = self._data.groupby('x')[['z', name]]
            data = []
            cols = list(dgroup.groups.keys())
            for key, grp in dgroup:
                data.append(np.interp(x=elev,
                                      xp=grp['z'],
                                      fp=grp[name],
                                      right=0))
                #cols.append(key)
            #
            # Expand wave
            cols2 = [item * -1 for item in cols[1:]]
            cols = np.hstack((list(reversed(cols2)), cols))
            data = np.vstack((list(reversed(data[1:])), data))
            #
            #data = self._data[['x', 'z', name]].copy()
            #data.insert(1, 'y', float(0.0))
            #cols = self._data['y']
            #
            #data.set_index(['x', 'y', 'z'], inplace=True)
            #data = data.to_xarray()
            #
            # krf on horizontal vel and acc
            if name in ['u', 'ut']:
                data *= krf
            #
            #data = np.transpose(data)
            # insert x axis to simulate a 3d wave [x, y, z]
            rows = np.array([cols[0], 0.0, cols[-1]])
            data = np.array([data, data, data])
            #data = np.expand_dims(data, axis=0)
            #data2 = np.append(data, data, axis=0)
            #cols = np.expand_dims(cols, axis=0)
            #rows = np.zeros((cols.shape))
            #
            #
            newname = title[i]
            kdf[newname] = xr.DataArray(data=data,
                                        coords=[rows, cols, elev],
                                        #coords={'row': (["x","y"], rows),
                                        #        'col': (["x","y"], cols),
                                        #        'z' : elev},
                                        dims=['y', 'x', 'z'],
                                        name=name)
        #
        ds = xr.Dataset(data_vars=kdf)
        #
        1 / 0
    #
    def get_kin4(self, elev: list, krf: list):
        """
        """
        kdf = self._get_xarray(elev, krf)
        # FIXME : rotate data according to wave direction
        #theta = self._theta
        return kdf
    #
    def _get_xarray(self, elev: list, krf: list):
        """ """
        #theta = self._theta
        #ctheta = np.cos(theta)
        #stheta = np.sin(theta)
        #
        items = ['u', 'ut', 'v', 'vt']
        title = ['ux', 'ax', 'uz', 'az']
        kdf = dict()
        for i, name in enumerate(items):
            dgroup = self._data.groupby('length')[['elevation', name]]
            data = []
            cols = list(dgroup.groups.keys())
            for key, grp in dgroup:
                data.append(np.interp(x=elev,
                                      xp=grp['elevation'],
                                      fp=grp[name],
                                      right=0))
            #
            # Expand wave
            cols2 = [item * -1 for item in cols[1:]]
            cols = np.hstack((list(reversed(cols2)), cols))
            data = np.vstack((list(reversed(data[1:])), data))
            #
            # krf on horizontal vel and acc
            if name in ['u', 'ut']:
                data *= krf
                # rotation according to theta angle (clockwise direction)
                #data = np.array([(ctheta + stheta) * data, # x' =  x*cos(theta) + y*sin(theta)
                #                 (ctheta - stheta) * data, # y' = -x*sin(theta) + y*cos(theta)
                #                 data])
            #else:
            #    data = np.array([data, data, data])
            #
            #data = np.transpose(data)
            # insert x axis to simulate a 3d wave [x, y, z]
            rows = np.array([cols[0], 0.0, cols[-1]])
            data = np.array([data, data, data])
            #data = np.expand_dims(data, axis=0)
            #data2 = np.append(data, data, axis=0)
            #cols = np.expand_dims(cols, axis=0)
            #rows = np.zeros((cols.shape))
            #
            #
            newname = title[i]
            kdf[newname] = xr.DataArray(data=data,
                                        coords=[rows, cols, elev],
                                        #coords={'row': (['y','x','z'], rows),
                                        #        'col': (['y','x','z'], cols),
                                        #        'z' : elev},
                                        dims=['y', 'x', 'z'],
                                        name=name)
        return xr.Dataset(data_vars=kdf)
#
#
#
def get_kinematic(n: int, z: list, B: list, Tanh: list, d: float,
                  surface: list, depth_points: int, 
                  is_finite: bool, g: exec = 9.80665):
    """
    n : order - Number of Fourier components or order of Stokes or cnoidal theory
    z : Solution vector
    B : Fourier coefficients
    TanH :
    d : water depth
    surface x :
    depth_points z :
    is_finite:
    g : 
    """
    g = g  # m/s^2
    pi = np.pi
    kd = z[1]
    c = z[4] / np.sqrt(z[1])
    ce = z[5] / np.sqrt(z[1])
    R = 1 + z[9] / z[1]
    #
    etas = surface['eta'].to_numpy()
    xx = surface['length'].to_numpy()
    phase = surface['phase'].to_numpy()
    time = surface['time'].to_numpy()
    #
    #npt = len(etas)
    #X = xx * kd / d # reset to dimensionless units
    #eta = [etas[j] * kd for j in range(npt)]
    #npoints = np.arange(points + 1) / points
    #npoints2 = repmat(npoints, m=npt, n=1)
    #
    depth_steps = np.arange(depth_points + 1) / depth_points
    #
    #y2 =  (1 + zdepth / d)
    #y2 = repmat(y, m=npt, n=1)
    #
    #
    #y4 = np.zeros((len(etas), len(zdepth)))
    #for i, eta in enumerate(etas):
    #    for j, point in enumerate(zdepth):
    #        if point > eta:
    #            if np.isclose(point, eta, atol=0.01):
    #                y4[i, j-1] = 1 + eta / d
    #            else:
    #                y4[i, j] = 1 + eta / d
    #            break
    #        else:
    #            y4[i, j] = 1 + point / d
    #
    if is_finite:
        # y = i * (1 + eta[j] / kd) / points
        #y = [[i / points * (1 + eta / d) for i in range(points + 1)]
        #     for eta in etas]
        #
        y = np.array([[point * (1 + eta / d) for point in depth_steps]
                      for eta in etas])
        #
        #output = [[Point(X[j], kd*((i*y[j])-1), kd, Tanh,B, n, ce, c, R, z, is_finite)
        #          for i in range(points)] for j in range(npt)]
        #
        #
        Y =  kd * (y - 1)
        output = pointkin(d, xx, Y, kd, Tanh, B, n, ce, c, R, z, is_finite, g)
        #
        #output = [Point(X[j], kd * (y[j][i] - 1), kd, Tanh, B, n, ce, c, R, z, is_finite)
        #          for j in range(npt) for i in range(points + 1)]
    else:
        # y = -pi + i / points * (eta[j] + pi)
        #
        #y = [[- pi + i / points * (etas / d + pi) for i in range(points + 1)]
        #      for eta in etas]
        #
        y = np.array([[- pi + point * (eta / d + pi) for point in depth_steps]
                      for eta in etas])
        #
        #
        output = pointkin(d, xx, y, kd, Tanh, B, n, ce, c, R, z, is_finite, g)
        #
        #output = [Point(X[j], y[j][i], kd, Tanh, B, n, ce, c, R, z, is_finite)
        #          for j in range(npt) for i in range(points + 1)]
    #
    # -----------------------------------------------------------
    # Dataframe setup
    #
    # [[xx[j], phase[j], time[j]]
    #outeq = [[xx[j], phase[j]]
    #         for j in range(npt) for i in range(points + 1)]
    #outeq = list(zip(*outeq))
    #dfkin = {'x': outeq[0], 'phase': outeq[1]} #'t': outeq[2]
    #
    #output = list(zip(*output))
    #output = [[row * factors[x] for row in col]
    #          for x, col in enumerate(output)]
    #dfkin.update({item: output[x] for x, item in enumerate(header)})
    #
    # df data format
    dfkin = {'x': repmat(xx, depth_steps.size, 1).flatten('F'),
             'phase': repmat(phase, depth_steps.size, 1).flatten('F'),
             'time': repmat(time, depth_steps.size, 1).flatten('F')}
    #dfkin.update({item: output[x].flatten('F') * factors[x]
    #              for x, item in enumerate(header)})
    dfkin.update(output)
    #
    #kindf(kin=dfkin, etas=etas, xx=xx, zdepth=zdepth)
    #
    df = DBframework()
    return df.DataFrame(dfkin)
    #return dfkin
#
#
def repmat2(A, n, axis:int):
    """
    """
    A1 = np.expand_dims(A, axis)
    A1 = np.transpose(A1)
    A1 = np.tile(A1, n)
    A1 = np.transpose(A1)
    return np.moveaxis(A1, 0, 1)
#
def permute2(A, order, axis:int=1):
    """ """
    A1 = repmat(A, order[0], axis)
    A1 = np.expand_dims(A1, axis=0)
    A1 = np.transpose(A1)
    A1 = np.tile(A1, order[1])
    return A1
    #return np.transpose(A1)
#
def permute(A, order, axis:int=1):
    """ """
    return np.transpose(repmat(A, order, axis))
#
#
def pointkin(d, x, Y, kd, Tanh, B, n, ce, c, R, z, Is_finite,
             g: float = 9.80665):
    """ """
    #
    X = x * kd / d # reset to dimensionless units
    #
    npoints = np.arange(n + 1)
    Xj = np.multiply.outer(X, npoints).T
    Yj = np.multiply.outer(Y, npoints).T
    #
    if Is_finite:
        coshdelta = np.cosh(Yj)
        sinhdelta = np.sinh(Yj)
        Tanh = permute2(Tanh, order=(sinhdelta.shape[1],
                                     sinhdelta.shape[2]))
        CH = coshdelta + sinhdelta * Tanh
        SH = sinhdelta + coshdelta * Tanh
    else:
        CH = np.exp(Yj)
        SH = np.exp(Yj)
    #
    yy = 1.0 + (Y / kd).T
    #Cos  =  np.cos(Xj)
    Cos = repmat2(np.cos(Xj), CH.shape[1], axis=0)
    Sin = repmat2(np.sin(Xj), CH.shape[1], axis=0)
    #
    B = permute2(B, order=(CH.shape[1], CH.shape[2]))
    npnt = permute2(npoints, order=(CH.shape[1], CH.shape[2]))
    #
    phi = np.sum(B * CH * Sin, axis=0)
    psi = np.sum(B * SH * Cos, axis=0)
    #
    u = np.sum(npnt * B * CH * Cos, axis=0)
    v = np.sum(npnt * B * SH * Sin, axis=0)
    #
    ux = np.sum(- npnt * npnt * B * CH * Sin, axis=0)
    vx = np.sum(npnt * npnt * B * SH * Cos, axis=0)
    #
    if Is_finite:
        header = ['z', 'u', 'v', 'dphidt', 'ut', 'vt', 'ux', 'uz', 'pressure', 'Bernoulli_check']
        factors = np.array([d, np.sqrt(g * d), np.sqrt(g * d), g * d, g, g,
                            np.sqrt(g / d), np.sqrt(g / d), g * d, 1])
        #
        # All PHI, PSI, u, v, ux and vx are dimensionless w.r.t. g & k.
        #Now convert to dimensionless w.r.t. d.
        #
        phi /= np.power(kd, 1.5)
        psi /= np.power(kd, 1.5)
        #
        u /= np.power(kd, 0.5)
        v /= np.power(kd, 0.5)
        ux *= np.power(kd, 0.5)
        vx *= np.power(kd, 0.5)
        #
        u = ce + u
        phi = ce * X + phi
        psi = ce * yy + psi
        dphidt = -c * u
        #
        ut = -c * ux
        vt = -c * vx
        uy = vx
        vy = -ux
        #
        dudt = ut + u*ux + v*uy
        dvdt = vt + u*vx + v*vy
        Pressure = (R - yy) - 0.5 * ((u - c) * (u - c) + v*v)
        Bernoulli_check = dphidt + Pressure + yy + 0.5*(u*u + v*v) - (R - 0.5*c*c)
    
    else:
        header = ['kz', 'u', 'v', 'dphidt', 'ut', 'vt', 'ux', 'uz', 'pressure', 'Bernoulli_check']
        factors = np.array([1 / kd, np.sqrt(g / kd), np.sqrt(g / kd), g / kd, g, g,
                            np.sqrt(g * kd), np.sqrt(g * kd), g / kd, 1])
        #
        u = z[5] + u
        phi = z[5] * X + phi
        dphidt = -z[4] * u
        ut = -z[4] * ux
        vt = -z[4] * vx
        uy = vx
        vy = -ux
        dudt = ut + u*ux + v*uy
        dvdt = vt + u*vx + v*vy
        Pressure = (z[9] - Y.T) - 0.5 * ((u - z[4]) * (u - z[4]) + v*v)
        Bernoulli_check = dphidt + Pressure + Y.T + 0.5*(u*u + v*v) - (z[9] - 0.5*z[4]*z[4])    
    #
    kinout = [yy, u, v, dphidt, ut, vt, ux, uy, Pressure, Bernoulli_check]
    kinout = np.array([item * factors[x] for x, item in enumerate(kinout)])
    # adjust water column coordinates
    kinout[0] -= d
    # name variables
    kinout = {item: kinout[x].flatten('F')
              for x, item in enumerate(header)}
    #print('-->')
    return kinout
#
#
def kindfX(kin, etas, xx, zdepth):
    """ """
    df = DBframework()
    dfkin = df.DataFrame(kin)
    grpkin = dfkin.groupby(['x'])
    #
    #y4 = np.zeros((len(etas), len(zdepth)))
    u = []
    for i, step in enumerate(xx):
        eta = etas[i]
        items = grpkin.get_group(name=step)
        for j, point in enumerate(zdepth):
            if point > eta:
                if np.isclose(point, eta, atol=0.01):
                    y4[i, j-1] = 1 + eta / d
                else:
                    y4[i, j] = 1 + eta / d
                break
            else:
                y4[i, j] = 1 + point / d    
    #
    dfkin
#
# -----------------------------------------------------
#  Velocities, accelerations, and pressure at a point
#
def Point(X, Y, kd, Tanh, B, n, ce, c, R, z, Is_finite):
    """ """
    #u = v = ux = vx = phi = psi = 0.
    psi = 0.0
    phi = 0.0
    vx = 0.0
    ux = 0.0
    v = 0.0
    u = 0.0
    y = 1. + Y/kd
    
    for j in range (1, n+1):
        Cos  = np.cos(j*X)
        Sin  = np.sin(j*X)
        if Is_finite:
            coshdelta = np.cosh(j*Y)
            sinhdelta = np.sinh(j*Y)
            C = coshdelta + sinhdelta*Tanh[j]
            S = sinhdelta + coshdelta*Tanh[j]
        else:
            #elif Is_deep:
            C = np.exp(j*Y)
            S = np.exp(j*Y)
        #
        phi += B[j] * C * Sin
        psi += B[j] * S * Cos
        u += j * B[j] * C * Cos
        v += j * B[j] * S * Sin
        ux += - j * j * B[j] * C * Sin
        vx += j * j * B[j] * S * Cos
    
    if Is_finite:
        # All PHI, PSI, u, v, ux and vx are dimensionless w.r.t. g & k.
        #Now convert to dimensionless w.r.t. d.
        phi /= pow(kd,1.5)
        psi /= pow(kd,1.5)
        u /= pow(kd,0.5)
        v /= pow(kd,0.5)
        ux *= pow(kd,0.5)
        vx *= pow(kd,0.5)
        u = ce + u
        phi = ce * X + phi
        psi = ce * y + psi
        dphidt = -c * u

        ut = -c * ux
        vt = -c * vx
        uy = vx
        vy = -ux
        dudt = ut + u*ux + v*uy
        dvdt = vt + u*vx + v*vy
        Pressure = R - y - 0.5 * ((u-c)*(u-c)+v*v)
        Bernoulli_check = dphidt + Pressure + y + 0.5*(u*u+v*v)-(R-0.5*c*c)
        #print("\n%f %f %f %f %f", R, y, 0.5*((u-c)*(u-c)+v*v),Pressure,Bernoulli_check)
    else:
    #elif Is_deep:
        u = z[5] + u
        phi = z[5] * X + phi
        dphidt = -z[4] * u
        ut = -z[4] * ux
        vt = -z[4] * vx
        uy = vx
        vy = -ux
        dudt = ut + u*ux + v*uy
        dvdt = vt + u*vx + v*vy
        Pressure = z[9] - Y - 0.5 * ((u-z[4])*(u-z[4])+v*v)
        Bernoulli_check = dphidt + Pressure + Y + 0.5*(u*u+v*v)-(z[9]-0.5*z[4]*z[4])
    #
    return y, u, v, dphidt, ut, vt, ux, uy, Pressure, Bernoulli_check
#
#
def get_kinematicX(n, z, Y, B, Tanh, nprofiles, points, is_finite):
    """ """
    pi = math.pi
    kd = z[1]
    c=z[4]/math.sqrt(z[1])
    ce=z[5]/math.sqrt(z[1])
    R=1+z[9]/z[1]    
    # Surface - print out Velocity and acceleration profiles plus check of Bernoulli
    #print("# %s\n", Title)
    #print("%s\n", Method)
    print("# Velocity and acceleration profiles and Bernoulli checks")
    if is_finite:
        print("# All quantities are dimensionless with respect to g and/or d")
    else:
        print("# All quantities are dimensionless with respect to g and/or k")
    #
    print("#*******************************************************************************")
    if is_finite:
        print("# y        u       v    dphi/dt   du/dt   dv/dt  du/dx   du/dy Bernoulli check")
        print("# -     -------------   -------  ------   -----  ------------- ---------------")
        print("# d        sqrt(gd)       gd        g       g       sqrt(g/d)        gd       ")
    else:
        print("# ky       u       v    dphi/dt   du/dt   dv/dt  du/dx   du/dy Bernoulli check")
        print("#       -------------   -------  ------   -----  ------------- ---------------")
        print("#         sqrt(g/k)       g/k       g       g       sqrt(gk)        g/k       ")
    print("#*******************************************************************************")
    print("# Note that increasing X/d and 'Phase' here describes half of a wave for")
    print("# X/d >= 0. In a physical problem, where we might have a constant x, because")
    print("# the phase X = x - c * t, then as time t increases, X becomes increasingly")
    print("# negative and not positive as passing down the page here implies.")
    print("#")
    #
    npt = number_steps(nprofiles)
    xx =  array('f', [0 for i in range(npt)])
    yy =  [] # array('f', [0 for i in range(npt)])
    eta = array('f', [0 for i in range(npt)])
    for j in range(npt):
        #X = 0.5 * L * (j / nprofiles)
        xx[j] = pi * (j / npt)
        eta[j] = Surface(xx[j], Y, n)
        yy.append([])
        print("# X/d = {: 8.4f}, Phase = {: 6.1f}".format(xx[j]/kd, xx[j] * 180 / pi ))
        for i in range(points):
            if is_finite:
                y = i * (1 + eta[j] / kd) / points
                yy[j].append((1 + eta[j] / kd) * i/points)
                (yyy, u, v, dphidt, ut, vt, ux, uy, 
                 Pressure, Bernoulli_check) = Point(xx[j], kd*(y-1), kd, Tanh,
                                                    B, n, ce, c, R, z, is_finite)                
            else:
                y = -pi + i / points * (eta[j] + pi)
                (yyy, u, v, dphidt, ut, vt, ux, uy, 
                 Pressure, Bernoulli_check) = Point(xx[j], y, kd, Tanh, 
                                                    B, n, ce, c, R, z, is_finite)

            print("{: 7.4f} {: 7.4f} {: 7.4f} {: 7.4f} {: 7.4f} {: 7.4f} {: 7.4f} {: 7.4f} {:7.0e}"
                  .format(y, u, v, dphidt, ut, vt, ux, uy, Bernoulli_check))
    print("#*******************************************************************************")
    #return xx, eta    
#