#
# Copyright (c) 2009-2022 steelpy
#
# Python stdlib imports
from array import array
from dataclasses import dataclass
import re
import math
from typing import NamedTuple, Tuple, Union, List, Dict

# package imports
from ..current.main_current import MeanCurrent
from ..operations.inout import (get_etas, get_kinematic, 
                               get_surface, get_kinematicX,
                               get_Height)
from steelpy.process.units.main import Units
#from steelpy.process.dataframe.dframe import DataFrame #, SeriesItem
from steelpy.process.dataframe.tb_deleted.dfmem import DataFrameBasic
#
#
#
#@dataclass
class WaveItem:
    __slots__ = ['H', 'Tw', 'd', 'title', 'Lw',
                 'order', 'nstep', 'niter', 'accuracy',
                 'current', 'c_type', 'method', 'finite_depth',
                 '_wave_length', '_z', '_Y', '_B', '_Tanh', '_Highest']
    
    def __init__(self, H:float, d:float, title:str, 
                 Tw:Union[float,None]=None, Lw:Union[float,None] = None, 
                 infinite_depth:bool=False,
                 current:float = 0.0, c_type:int = 1,
                 order:int=5, nstep:int=2,
                 niter:int=40, accuracy:float=1e-6):
        """ """
        self.H = H
        if Lw:
            self.Lw = Lw
            self.Tw = None
        else:
            self.Tw = Tw
            self.Lw = None
        self.d = d
        self.title: str
        # current 
        self.current = current
        self.c_type = c_type        
        #
        self.order = order
        self.nstep = nstep
        self.niter = niter
        self.accuracy = accuracy
        #
        self.finite_depth = True
        if infinite_depth:
            self.finite_depth = False
        #
        #self.wave_length:float
        #
    #
    def get_parameters(self, g:float = 9.80665):
        """ 
        g: gravity # m/s^2
        
        returns:
        MaxH : H/d
        T : Dimensionless period / None
        L : Dimensionless wavelength / None
        c_type  : Current criterion, 1 for Euler, or 2 for Stokes
        current : Current magnitude
        order :  Number of Fourier components or order of Stokes or cnoidal theory
        nsteps : Number of height steps to reach H/d
        niter  : Maximum number of iterations for each step (20)
        crit   : Criterion for convergence (1e-6)
        Height : ??
        finite_depth : True/False
        """
        MaxH = self.H / self.d
        current = self.current / math.sqrt( g * self.d )
        if self.Lw:
            case = "wavelength"
            L = self.Lw / self.d
            Height = get_Height(MaxH, case, self.finite_depth, L=L)
            T = None
        else:
            case = 'period'
            T = self.Tw * math.sqrt( g /  self.d )
            Height = get_Height(MaxH, case, self.finite_depth, T=T)
            L = None
        #
        return [MaxH, case, T, L, self.c_type, current,  
                self.order, self.nstep, self.niter, self.accuracy,
                Height, self.finite_depth]
    #
    @property
    def L(self):
        """ Wave Length"""
        try:
            return self._wave_length * self.d
        except AttributeError:
            self.__call__()
            return self._wave_length * self.d
    #
    #
    def surface(self, surface_points:int = 36,
                step_size=None):
        """Free surface (wave) profile
        surface_points : Number of points on free surface (the program clusters them near crest)
        step_size: deg
        """
        self.__call__()
        n = self.order
        kd = self._z[1]
        #method = self.method
        #x, eta = print_elevation(n, z, Y, B, Tanh, surface_points, method)
        # n, z, Y, B, Tanh, nprofiles, points, is_finite
        #surface_points = 36
        #points = 20
        #x, eta = get_etas(n, self._z, self._Y, self._B, self._Tanh, 
        #                  surface_points, self.finite_deep)
        #
        x, eta, phase = get_surface(n, kd, self._Y, 
                                    surface_points, self.finite_depth)
        #
        #dataframe = {'x': x, 'eta': eta, 'phase': phase,
        #             'zeta': [(1+item) for item in eta]}
        #
        dataframe = {'x': [item*self.d for item in x],
                     'eta': [item*self.d for item in eta],
                     'phase': phase,
                     'zeta': [(1+item)*self.d for item in eta]}
        #eta = [self.d + item for item in eta]
        return SurfaceResults(dataframe, self.Tw, self.d, self.finite_depth)
    #
    def kinematics(self, surface=None, depth_points:int = 20,
                   surface_points:int = 36):
        """ """
        self.__call__()
        if not surface:
            surface = self.surface(surface_points)        
        n = self.order
        #method = self.method
        #
        #get_kinematicX(n, self._z, self._Y, self._B, self._Tanh,
        #               surface_points, depth_points, self.finite_depth)
        #
        #
        kout = get_kinematic(n, self._z, self._B, self._Tanh, self.d, 
                             surface['eta'], surface['x'],
                             depth_points, self.finite_depth)
        #
        dataframe = KinematicResults(kout, depth_points+1, self.finite_depth)
        #kd = self._z[1]
        dataframe['x'] = surface['x'] #[item * self.d for item in surface.x_kd]
        dataframe['phase'] = surface['phase'] #dataframe['x'] * kd * 180 / math.pi
        return dataframe
        #return KinematicResults (dataframe,
        #                        depth_points+1, self.finite_depth)
        #print ( 'end' )
    #
    #
#
#
#
class WaveRegModule:
    __slots__ = [ 'order', 'nsteps', 'max_iter', 'accuracy',
                  '_labels', '_cases', '_current', # 'c_type', 
                  'infinite_depth']

    def __init__(self, n:int=5, nstep:int=2,
                 number:int=40, accuracy:float=1e-6):
        """
        n      : Stokes order (5)
        nstep  : Number of height steps to reach H/d (2)
        number : Maximum number of iterations for each step (20)
        accuracy   : Criterion for convergence
        """
        self.order = n
        self.nsteps = nstep
        self.max_iter = number
        self.accuracy = accuracy
        #self.c_type = 1
        #self._current = 0 # m/s
        self.infinite_depth = False
        #
        self._labels: List = [ ]
        self._cases: List = [ ]
        self._current = MeanCurrent()
    #
    def __getitem__(self, case_name: Union[int, str]) -> Tuple:
        """
        case_name : Wave name
        """
        try:
            index = self._labels.index(case_name)
            return self._cases[index]
        except ValueError:
            raise IndexError('   *** wave {:} does not exist'.format(case_name))
    #
    #
    @property
    def infinite_water_depth(self):
        """ Water of infinite depth"""
        self.infinite_depth = True
    #
    @property
    def mean_current(self):
        """ """
        return self._current
#
#
class SurfaceResults(DataFrameBasic):
    
    def __init__(self, data:Union[None, Dict], 
                 d:float, t:float, finite_depth:bool):
        """
        zeta:
        x  : 
        d  : Water depth
        k  : wave number
        """
        super().__init__(data)
        #eta_kd :List
        #A : List
        #w : List
        #phase : List
        #eta: List
        #x_kd : List
        self.t = t
        self.d = d
        #order:int=1
        #kd:float
        self.finite_depth = finite_depth
        #
        #def k(self):
        #    """Wave number"""
        #    #k, niter = wavenum(self.w, self.d)
        #    kd = self.z[1]
        #    return k
    #
    #@property
    #def eta(self):
    #    """ """
    #    return [xx*self.d for xx in self.eta_kd]
    ##
    #@property
    #def zeta(self):
    #    """ """
    #    return [(1+xx)*self.d for xx in self.eta_kd]
    ##
    #@property
    #def phase(self):
    #    """ """
    #    return [xx * 180 / math.pi for xx in self._data['x']]
    ##
    #@property
    #def x(self):
    #    """ """
    #    return [self.d*item for item in self.x_kd]
    #
    #
    #
    def plot(self):
        """ """
        import matplotlib.pyplot as plt
        #
        x = [-1*item for item in reversed(self._data['x'][1:])]
        x.extend(self._data['x'])
        #
        y = [item for item in reversed(self._data['zeta'][1:])]
        y.extend(self._data['zeta'])
        #
        plt.plot(x, y)
        plt.title('Surface')
        plt.xlabel('$\lambda$ (m)')
        #if self.finite_depth:
        #    plt.ylabel('$\zeta$ (m)')
        #else:
        plt.ylabel('$\zeta$ (m)')
        plt.show()
    #
    def printout(self):
        """ """
        print("#")
        print("# Surface of wave - trough-crest-trough,")
        if self.finite_depth:
            header1 = "X"
            header2 = "theta"
            #print("# Non-dimensionalised with respect to depth")
        else:
            header1 = "kX"
            header2 = "k eta"
            #print("# Non-dimensionalised with respect to wavenumber")
        #
        npt = len(self._data['eta'])
        for j in range(npt):
            print("# {:} = {: 8.4f}, Phase = {: 6.1f} {:} = {: 8.4f}"
                  .format(header1, self._data['x'][j], 
                          self._data['phase'][j], 
                          header2, self._data['eta'][j]))
#
#
#
#class KinematicResults(NamedTuple):
#    """ """
#    dataframe:List
#    #x:List
#    #kd:float
#    depth_points:int
#    finite_depth:bool
#
# (DataFrameBasic)
class KinematicResults(DataFrameBasic):
    
    def __init__(self, data:Union[None, Dict], 
                 depth_points:int, finite_depth:bool):
        """
        """
        super().__init__(data)

        self.depth_points = depth_points
        self.finite_depth = finite_depth
    
    def get_data(self, name):
        """ """
        x = self._data['x']
        #coldata = [item*self.factors[name] for item in self.dataframe[name]]
        items = to_matrix(self._data[name], n=self.depth_points)
        new_data = {step: items[i] for i, step in enumerate(x)}
        return DataFrameBasic(new_data)
        #raise IOError(f'item {name:} not found')
    #
    @property
    def z(self):
        """ 
        z  : finite water depth
        kz : infinite water depth
        """
        #xxx = self.kout
        #new_data = self.get_data('kz')
        try:
            return self.get_data('z')
        except IOError:
            return self.get_data('kz')
        
    #
    @property
    def ux(self):
        """ """
        return self.get_data('u')
    #
    #
    @property
    def uz(self):
        """ """
        return self.get_data('v')
    #
    @property
    def ax(self):
        """ """
        return self.get_data('ut')
    #
    @property
    def az(self):
        """ """
        return self.get_data('vt')   
    #
    @property
    def pressure(self):
        """ """
        return self.get_data('pressure')    
    #
    def printout(self):
        """ """
        #pi = math.pi
        x = self._data['x']
        phase = self._data['phase']
        #z = self.z
        #xx = self.dataframe
        headers = list(self._data.keys())
        #headers.pop()
        #del headers[10:] # remove x and phase
        surface_points = len(x)
        index = surface_points * self.depth_points
        #items = len(headers)
        #
        print("#")
        print("# Velocity and acceleration profiles and Bernoulli checks")
        if self.finite_depth:
            header = "X/d"
            #print("# All quantities are dimensionless with respect to g and/or d")
            print("#*******************************************************************************")
            print("# y        u       v    dphi/dt   du/dt   dv/dt  du/dx   du/dy Bernoulli check")
            #print("# -     -------------   -------  ------   -----  ------------- ---------------")
            #print("# d        sqrt(gd)       gd        g       g       sqrt(g/d)        gd       ")            
        else:
            header = "kX"
            #print("# All quantities are dimensionless with respect to g and/or k")
            print("#*******************************************************************************")
            print("# ky       u       v    dphi/dt   du/dt   dv/dt  du/dx   du/dy Bernoulli check")
            #print("#       -------------   -------  ------   -----  ------------- ---------------")
            #print("#         sqrt(g/k)       g/k       g       g       sqrt(gk)        g/k       ")
        #
        print("#*******************************************************************************")
        print("# Note that increasing X/d and 'Phase' here describes half of a wave for")
        print("# X/d >= 0. In a physical problem, where we might have a constant x, because")
        print("# the phase X = x - c * t, then as time t increases, X becomes increasingly")
        print("# negative and not positive as passing down the page here implies.")
        print("#")        
        #
        step = 0
        j = 0
        for i in range(index):
            if step == i:
                print("# {:} = {: 8.4f}, Phase = {: 6.1f}"
                      .format(header, x[j], phase[j] ))
                      #.format(header, x[j], self.kd*x[j] * 180 / pi))
                step += self.depth_points
                j += 1
            test = [f'{self._data[key][i]: 1.3e}' for key in headers[:10]]
            print(*test, sep=' ')
        #print('--')
    #
    def plot(self):
        """ """
        import matplotlib.pyplot as plt
        #xxx = self.kout
        zlev = self.z  
        #
        velh = self.ux
        velv = self.uz
        acch = self.ax
        accv = self.az
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
        import matplotlib.pyplot as plt
        #
        velh = self.ux
        velv = self.uz        
        #xx = self.x
        zlev = self.z
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
def zeros(m, n=False, code='d'):
    """
    Create zero matrix
    """
    if n:
        new_matrix = [array(code, [0 for row in range(n)])
                      for col in range(m)]
    else:
        new_matrix = array(code, [0 for row in range(m)])
    return new_matrix
#
def to_matrix(l, n):
    return [l[i:i+n] for i in range(0, len(l), n)]
#
def get_wave_data(case_data):
    """ """
    if isinstance(case_data, (list, tuple)):
        data = get_data_list(data=case_data)
    elif isinstance(case_data, dict):
        data = get_data_dic(data=case_data)
    else:
        raise IOError('input data not valid')
    return data
#
def get_data_dic(data:Dict)->List[float]:
    """[case, load, theta, phi, rho, phase]"""
    new_data = [0,0,0,0]
    for key, item in data.items():
        if re.match(r"\b((wave)?h(eight)?(w)?)\b", key, re.IGNORECASE):
            new_data[0] = item.value
        elif re.match(r"\b((wave)?period|t(w)?|s)\b", key, re.IGNORECASE):
        #elif re.match(r"\b(t(w)?|s)\b", key, re.IGNORECASE):
            new_data[1] = item.value
        elif re.match(r"\b((water)?d(epth)?)\b", key, re.IGNORECASE):
            new_data[2] = item.value
        elif re.match(r"\b((wave)?l(ength)?)\b", key, re.IGNORECASE):
            new_data[3] = item.value
        #elif re.match(r"\b(rho)\b", key, re.IGNORECASE):
        #    new_data[4] = item.value        
        #elif re.match(r"\b(phase)\b", key, re.IGNORECASE):
        #    new_data[5] = item.value        
    return new_data
#
#
def get_data_list(data:List[Units], steps:int=4)->List[float]:
    """
    [Hw, Tw, d, Lw]
    """
    new_data = [0] * steps
    for x, item in enumerate(data):
        new_data[x] = item.value
    
    return new_data
    
#