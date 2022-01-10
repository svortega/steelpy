#
# Copyright (c) 2009-2022 steelpy
#
# Python stdlib imports
from array import array
from dataclasses import dataclass
import math
from typing import NamedTuple, Tuple, Union, List, Dict

# package imports
from steelpy.metocean.regular.stokes.Inout import (print_output, print_elevation,
                                                   print_flowfield, print_velo)


#
#
class StokesModule:
    __slots__ = [ 'order', 'nsteps', 'max_iter', 'accuracy',
                  '_labels', '_cases', 'c_type', '_current']

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
        self.c_type = 1
        self._current = 0 # m/s
        #
        self._labels: List = [ ]
        self._cases: List = [ ]
    #
    def __setitem__(self, case_name: int,
                    case_data: Union[List[float], Dict[str, float]]) -> None:
        """
        case_name : Wave name
         H : Wave height [unit length]
         T : Wave Period [second]
         d : Water depth LTH + Tide and Surge [unit length]
         Lw : Wave Length [unit length]
         Phase : Wave phase [degree]
         Order : ??
        """
        try:
            self._labels.index(case_name)
            raise Exception('    *** warning wave {:} already exist'
                            .format(case_name))
        except ValueError:
            self._labels.append(case_name)
            self._cases.append(WaveItem(H=case_data[0].value, T=case_data[1].value,
                                        d=case_data[2].value, title=case_name,
                                        order=self.order, nstep=self.nsteps, number=self.max_iter,
                                        accuracy=self.accuracy))

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
@dataclass
class WaveItem:
    H:float
    T:Union[float,None]
    d:float
    title: str
    #
    order: int = 5
    nstep: int = 2
    number: int = 40
    accuracy: float = 1e-6
    #
    Lw: Union[ float, None ] = None
    current:float = 0.0
    c_type: int = 1
    #
    def elevation(self):
        pass
    #
    def flowfield(self, nprofiles:int = 16, points:int = 20):
        """
        nprofiles : Number of velocity/acceration profiles over half a wavelength to print out,
                including x= 0 and x = lambda/2
        points : Number of vertical points in each profile, including points at bottom and surface.
        """
        z, Y, B, Tanh = self.__call__ ()
        n = self.order
        print_flowfield (n, z, Y, B, Tanh, nprofiles, points)
        print( 'end' )
    #
    def surface(self, surface_points:int = 50):
        """Free surface (wave) profile
        surface_points : Number of points on free surface (the program clusters them near crest)
        """
        z, Y, B, Tanh = self.__call__()
        n = self.order
        x, eta = print_elevation(n, z, Y, B, Tanh, surface_points)
        return [x, eta]
        #print ( 'end' )
    #
    def velocities(self, surface_points:int = 50):
        """ """
        z, Y, B, Tanh = self.__call__ ()
        n = self.order
        H = self.H
        print_velo(H, n, z, Y, B, Tanh, surface_points)
        print ( 'end' )
    #
    def __call__(self):
        """ Solver """
        z, Y, B, Tanh = mainStokes(h=self.H, t=self.T, d=self.d, L=self.Lw,
                                   current=self.current, c_type=self.c_type,
                                   n=self.order, nstep=self.nstep,
                                   number=self.number, accuracy=self.accuracy)
        return z, Y, B, Tanh
#
#
def mainStokes(h:float, t:Union[float,None], d:float, L:Union[float,None],
               current:float, c_type:int=1,
               n:int=5, nstep:int=2, number:int=40, accuracy:float=1e-6):
    """
    Stokes theory calculations
        Input
    -------
    d : water depth (m)
    h : wave height (m)
    t : wave period (s)
    c : current magnitud (m/s)

    Optional
    --------
    c_type : 1 - Eularian mean, 2 - depth integrated mean (1)
    n      : Stokes order (5)
    nstep  : Number of height steps to reach H/d (2)
    number : Maximum number of iterations for each step (20)
    crit   : Criterion for convergence
    """
    #
    # Read_data()
    # inital values
    g = 9.80665  # m/s^2
    maxh = float ( h ) / float ( d )
    T = float ( t ) * math.sqrt ( g / float ( d ) )
    # _Height = maxh / (T * T)
    # wave = 1
    current = float ( current ) / math.sqrt ( g * float ( d ) )
    c_type = int ( c_type )
    #
    Case = 'Period'
    if L:
        Case = "Wavelength"

    z = zeros ( 2 * n + 10 + 1 )
    Y = zeros ( n + 1 )
    # coeff = zeros(n + 1)
    e = zeros ( n + 1 )
    # monitor = stdout
    H = maxh
    pi = math.pi

    if Case == 'Wavelength':
        kd = 2.0 * pi / L
        kH = kd * H
        ckd, skd, ss, t, C, D, E = CDE ( kd )
    else: # Period
        # if period is specified, solve dispersion relation using secant method
        # Until February 2015 the bisection method was used for this.
        # I found that in an extreme case (large current) the bracketting
        # of the solution was not correct, and the program failed,
        # without printing out a proper error message.
        print ( "# Period has been specified." )
        print ( "# Now solving for L/d _iteratively, printing to check convergence" )
        omega = 2 * pi / T
        # Fenton & McKee for initial estimate
        kFM = (omega * omega *
               pow ( 1 / math.tanh ( pow ( omega, 1.5 ) ),
                     (2.0 / 3.0) ))
        kd1 = kFM
        kd2 = kFM * 1.01
        ckd, skd, ss, t, C, D, E = CDE ( kd2 )
        F2 = F ( kd2, H, T, current, c_type, C, n, D )

        for _iter in range ( 1, number + 1 ):
            ckd, skd, ss, t, C, D, E = CDE ( kd1 )
            F1 = F ( kd1, H, T, current, c_type, C, n, D )
            Fd = (F2 - F1) / (kd2 - kd1)
            delta = F1 / Fd
            kd2 = kd1
            kd1 = kd1 - delta
            print ( "{: 8.4f}".format ( 2 * pi / kd1 ) )
            if abs ( delta / kd1 ) < accuracy:
                break
            F2 = F1
            if _iter >= number:
                raise ValueError ( "Secant for solution of wavenumber has not converged" )
                # print("Contact John Fenton johndfenton@gmail.com")
        # Next _iter
        kd = kd1
        kH = kd * H
    #
    z[ 1 ] = kd
    z[ 2 ] = kH
    SU = 0.5 * kH / pow ( kd, 3 )
    print ( "# Stokes-Ursell no.: {:}".format ( SU ) )
    if SU > 0.5:
        print ( " > 1/2. Results are unreliable" )
    else:
        print ( " < 1/2, Stokes theory should be valid" )
    #
    e[ 1 ] = 0.5 * kH
    for i in range ( 2, n + 1 ):
        e[ i ] = e[ i - 1 ] * e[ 1 ]
    # Calculate coefficients
    Y, z, A, B = AB ( skd, ss, t, n, z, e, C, kd, ckd )
    z[ 7 ] = C[ 0 ] + e[ 2 ] * C[ 2 ] + e[ 4 ] * C[ 4 ]  # ubar
    z[ 8 ] = - e[ 2 ] * D[ 2 ] - e[ 4 ] * D[ 4 ]
    z[ 9 ] = 0.5 * C[ 0 ] * C[ 0 ] + e[ 2 ] * E[ 2 ] + e[ 4 ] * E[ 4 ]

    if c_type == 1:
        z[ 5 ] = current * math.sqrt ( kd )
        z[ 4 ] = z[ 7 ] + z[ 5 ]
        z[ 6 ] = z[ 4 ] + z[ 8 ] / kd - z[ 7 ]

    elif c_type == 2:
        z[ 6 ] = current * math.sqrt ( kd )
        z[ 4 ] = z[ 6 ] - z[ 8 ] / kd + z[ 7 ]
        z[ 5 ] = z[ 4 ] - z[ 7 ]
    #
    if Case == 'Wavelength':
        z[ 3 ] = 2 * pi / z[ 4 ]

    elif Case == 'Period':
        z[ 3 ] = T * math.sqrt ( kd )

    # Tanh = zeros(n + 1)
    # for i in range(1, n + 1):
    #    Tanh[i] = math.tanh(i * z[1])
    #
    Tanh = [ math.tanh ( i * z[ 1 ] ) for i in range ( n + 1 ) ]
    print_output (n, z, Y, B, c_type, current )
    return z, Y, B, Tanh
#
# 
# Main program
def StokesMain(d:float, h:float, t:float, current:float,
               c_type:int=1, n:int=5, nstep:int=2,
               number:int=40, accuracy:float=1e-6):
    """
    Stokes theory calculations
        Input
    -------
    d : water depth (m)
    h : wave height (m)
    t : wave period (s)
    c : current magnitud (m/s)

    Optional
    --------
    c_type : 1 - Eularian mean, 2 - depth integrated mean (1)
    n      : Stokes order (5)
    nstep  : Number of height steps to reach H/d (2)
    number : Maximum number of iterations for each step (20)
    crit   : Criterion for convergence
    """
    #
    # Read_data()
    # inital values
    g = 9.80665  # m/s^2
    maxh = float(h) / float(d)
    T = float(t) * math.sqrt(g / float(d))
    #_Height = maxh / (T * T)
    # wave = 1
    current = float(current) / math.sqrt(g * float(d))
    c_type = int(c_type)
    Case = 'Period'

    z = zeros(2 * n + 10 + 1)
    Y = zeros(n + 1)
    #coeff = zeros(n + 1)
    e = zeros(n + 1)
    # monitor = stdout
    H = maxh
    pi = math.pi

    if Case == 'Wavelength':
        kd = 2.0 * pi / L
        kH = kd * H
        ckd, skd, ss, t, C, D, E = CDE(kd)

    elif Case == 'Period':
        # if period is specified, solve dispersion relation using secant method
        # Until February 2015 the bisection method was used for this.
        # I found that in an extreme case (large current) the bracketting
        # of the solution was not correct, and the program failed,
        # without printing out a proper error message.
        print("# Period has been specified.")
        print("# Now solving for L/d _iteratively, printing to check convergence")
        omega = 2 * pi / T
        # Fenton & McKee for initial estimate
        kFM = (omega * omega *
               pow(1 / math.tanh(pow(omega, 1.5)),
                   (2.0 / 3.0)))
        kd1 = kFM
        kd2 = kFM * 1.01
        ckd, skd, ss, t, C, D, E = CDE(kd2)
        F2 = F(kd2, H, T, current, c_type, C, n, D)

        for _iter in range(1, number + 1):
            ckd, skd, ss, t, C, D, E = CDE(kd1)
            F1 = F(kd1, H, T, current, c_type, C, n, D)
            Fd = (F2 - F1) / (kd2 - kd1)
            delta = F1 / Fd
            kd2 = kd1
            kd1 = kd1 - delta
            print("{: 8.4f}".format(2 * pi / kd1))
            if abs(delta / kd1) < accuracy:
                break
            F2 = F1
            if _iter >= number:
                raise ValueError("Secant for solution of wavenumber has not converged")
                #print("Contact John Fenton johndfenton@gmail.com")
        # Next _iter
        kd = kd1
        kH = kd * H
    z[1] = kd
    z[2] = kH
    SU = 0.5 * kH / pow(kd, 3)
    print("# Stokes-Ursell no.: {:}".format(SU))
    if SU > 0.5:
        print(" > 1/2. Results are unreliable")
    else:
        print(" < 1/2, Stokes theory should be valid")
    #
    e[1] = 0.5 * kH
    for i in range(2, n + 1):
        e[i] = e[i - 1] * e[1]
    # Calculate coefficients
    Y, z, A, B = AB(skd, ss, t, n, z, e, C, kd, ckd)
    z[7] = C[0] + e[2] * C[2] + e[4] * C[4]  # ubar
    z[8] = - e[2] * D[2] - e[4] * D[4]
    z[9] = 0.5 * C[0] * C[0] + e[2] * E[2] + e[4] * E[4]

    if c_type == 1:
        z[5] = current * math.sqrt(kd)
        z[4] = z[7] + z[5]
        z[6] = z[4] + z[8] / kd - z[7]

    elif c_type == 2:
        z[6] = current * math.sqrt(kd)
        z[4] = z[6] - z[8] / kd + z[7]
        z[5] = z[4] - z[7]
    #
    if Case == 'Wavelength':
        z[3] = 2 * pi / z[4]

    elif Case == 'Period':
        z[3] = T * math.sqrt(kd)
    
    #Tanh = zeros(n + 1)
    #for i in range(1, n + 1):
    #    Tanh[i] = math.tanh(i * z[1])
    #
    Tanh = [math.tanh(i * z[1]) for i in range(n + 1)]

    #	Output results and picture of wave

    # Solution=fopen("Solution.res","w")
    # Elevation = fopen("Surface.res","w")
    # Flowfield = fopen("Flowfield.res","w")

    # Output()
    #method = ("Stokes method order {:}".format(n))
    L = print_output(n, z, Y, B, c_type, current)

    surface_points = 50
    print_elevation(n, z, Y, B, Tanh, surface_points)

    nprofiles = 16
    points = 20
    print_flowfield(n, z, Y, B, Tanh, nprofiles, points)

    # fflush(Nothing)
    print("Touch key to continue ")
    # Console.ReadKey(True).KeyChar
    print("Finished")
    # End Function


# End Class
#
# 
def F(kd:float, H:float, T:float, Current:float, Current_criterion:int,
      C:array, n:int, D:array):
    """
    Evaluates dispersion relation - used in iterative solution for wavelength
    """
    pi = math.pi
    kH = kd * H
    e2 = 0.5 * kH * 0.5 * kH
    rhs = Current * math.sqrt(kd) - 2.0 * pi / T / math.sqrt(kd)
    CC = zeros(3)
    CC[0] = C[0]
    CC[2] = 0.0
    CC[1] = CC[2]
    if n >= 3:
        CC[1] = e2 * C[2]
    #
    if n >= 5:
        CC[2] = e2 * e2 * C[4]
    
    if Current_criterion == 2:
        if n >= 3:
            CC[1] += e2 * D[2] / kd
        
        if n >= 5:
            CC[2] += e2 * e2 * D[4] / kd
    #
    CC[1] = CC[0] + CC[1]
    CC[2] = CC[1] + CC[2]
    # return rhs+CC[2]; // Should this be current criterion rather than 2?
    return rhs + CC[Current_criterion]


#
def CDE(kd:float):
    """
    Calculate coefficient arrays C[], D[] and E[]
    """
    skd = math.sinh(kd)
    ckd = 1.0 / math.tanh(kd)
    tkd = math.tanh(kd)
    S = 1.0 / math.cosh(2.0 * kd)
    ss = [0, S]
    for i in range(2, 8 + 1):
        ss.append(ss[i - 1] * S)
    #
    t = [0, 1.0 - S]
    for i in range(2, 8 + 1):
        t.append(t[i - 1] * (1.0 - S))
    #
    C = zeros(5)
    D = zeros(5)
    E = zeros(5)
    C[0] = pow(tkd, 0.5)
    C[2] = pow(tkd, 0.5) * (2 + 7 * ss[2]) / (4 * t[2])
    C[4] = pow(tkd, 0.5) * (4 + 32 * ss[1] - 116 * ss[2] - 400 * ss[3] - 71 * ss[4] + 146 * ss[5]) / (32 * t[5])
    D[2] = -pow(ckd, 0.5) / 2
    D[4] = pow(ckd, 0.5) * (2 + 4 * ss[1] + ss[2] + 2 * ss[3]) / (8 * t[3])
    E[2] = tkd * (2 + 2 * ss[1] + 5 * ss[2]) / (4 * t[2])
    E[4] = tkd * (8 + 12 * ss[1] - 152 * ss[2] - 308 * ss[3] - 42 * ss[4] + 77 * ss[5]) / (32 * t[5])
    return ckd, skd, ss, t, C, D, E
#
# 
def AB(skd:float, ss:List, t:List, n:int,
       z:array, e:array, C:array, kd:float, ckd:float):
    """
    Calculate coefficient arrays A[] and B[] and Fourier coefficients
    """
    BB = zeros(n + 1, n + 1)
    A = zeros(n + 1, n + 1)
    B = zeros(n + 1)
    Y = zeros(n + 1)
    #
    A[1][1] = 1 / skd
    A[2][2] = 3*ss[2] / (2*t[2])
    A[3][1] = (-4 - 20*ss[1] + 10*ss[2] - 13*ss[3]) / (8*skd*t[3])
    A[3][3] = (-2*ss[2] + 11*ss[3]) / (8*skd*t[3])
    A[4][2] = (12*ss[1] - 14*ss[2] - 264*ss[3] - 45*ss[4] - 13*ss[5]) / (24*t[5])
    A[4][4] = (10*ss[3] - 174*ss[4] + 291*ss[5] + 278*ss[6]) / (48*(3 + 2*ss[1])*t[5])
    A[5][1] = ((-1184 + 32*ss[1] + 13232*ss[2] + 21712*ss[3] + 20940*ss[4] 
                + 12554*ss[5] - 500*ss[6] - 3341*ss[7] - 670*ss[8]) 
               / (64*skd*(3 + 2*ss[1]) * (4 + ss[1])*t[6]))
    A[5][3] = ((4*ss[1] + 105*ss[2] + 198*ss[3] - 1376*ss[4] - 1302*ss[5] - 117*ss[6] + 58*ss[7]) 
               / (32*skd*(3 + 2*ss[1])*t[6]))
    A[5][5] = ((-6*ss[3] + 272*ss[4] - 1552*ss[5] + 852*ss[6] + 2029*ss[7] + 430*ss[8]) 
               / (64*skd*(3 + 2*ss[1])*(4 + ss[1])*t[6]))

    for i in range(1, n + 1):
        z[n + 10 + i] = 0.0
        jj = ((i + 1) % 2) + 1
        for j in range(jj, n, 2):
            z[n + 10 + i] += A[j][i] * e[j]
        #
        z[n + 10 + i] *= C[0] * math.cosh(i * kd)
        B[i] = z[n + 10 + i]
    #
    BB[1][1] = 1.0
    BB[2][2] = ckd*(1 + 2*ss[1]) / (2*t[1])
    BB[3][1] = -3*(1 + 3*ss[1] + 3*ss[2] + 2*ss[3]) / (8*t[3])
    BB[3][3] = -BB[3][1]
    BB[4][2] = (ckd*(6 - 26*ss[1] - 182*ss[2] - 204*ss[3] - 25*ss[4] + 26*ss[5]) 
                / (6*(3 + 2*ss[1])*t[4]))
    BB[4][4] = (ckd*(24 + 92*ss[1] + 122*ss[2] + 66*ss[3] + 67*ss[4] + 34*ss[5]) 
                / (24*(3 + 2*ss[1])*t[4]))
    BB[5][3] = (9*(132 + 17*ss[1] - 2216*ss[2] - 5897*ss[3] - 6292*ss[4] 
                   - 2687*ss[5] + 194*ss[6] + 467*ss[7] + 82*ss[8]) 
                / (128*(3 + 2*ss[1])*(4 + ss[1])*t[6]))
    BB[5][5] = (5*(300 + 1579*ss[1] + 3176*ss[2] + 2949*ss[3] + 1188*ss[4] 
                   + 675*ss[5] + 1326*ss[6] + 827*ss[7] + 130*ss[8]) 
                / (384*(3 + 2*ss[1])*(4 + ss[1])*t[6]))
    BB[5][1] = -(BB[5][3] + BB[5][5])

    for i in range(1, n + 1):
        Y[i] = 0.0
        j = ((i + 1) % 2) + 1
        while j <= n:
            Y[i] += BB[j][i] * e[j]
            j += 2
    #
    return Y, z, A, B
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
#
#
#
if __name__ == "__main__":
    #
    print("--------------------------------------------")
    print('        Stokes Method Wave Theory')
    print('                 INPUT')
    print("--------------------------------------------")
    # dw = input('water depth (m) = ')
    dw = 139
    # tw = input('wave period (s) = ')
    tw = 19.2
    # hw = input("wave height (m) = ")
    hw = 44.2

    print("current type : ")
    print("(1) Eularian-mean")
    print("(2) depth-integrated")
    # u_type = input("input [1] = ")
    u_type = 1
    if not u_type:
        u_type = 2
    # u = input('current magnitud (m/s) [0] = ')
    u = 0
    if not u:
        u = 0.0
    #
    StokesMain(dw, hw, tw, u, u_type)
    #
    if not u:
        u = 0
    # PhaseStep = input('Number of grid points [21] = ')
    # if not PhaseStep:
    #    PhaseStep = 21
    #
    print("")
    print("Kinematics along : ")
    # ans = input('(H)orizontal, (V)ertical, (S)urface [S]: ')
    # ans = findKtic(ans)


