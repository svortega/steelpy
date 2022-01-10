#
# Copyright (c) 2009-2022 steelpy
#
# Python stdlib imports
from array import array
import math

# package imports
#import matplotlib.pyplot as plt

#
#
def title_block(is_finite, c_type, current, z):
    """
    PRINT OUT TITLE BLOCK
    """
    #
    #  Highest wave - eqn (32) of Fenton (1990)
    #
    pi = math.pi
    L = 2 * pi / z[1]
    Highest = ((0.0077829 * L * L * L + 0.0095721 * L * L + 0.141063 * L)
               / (0.0093407 * L * L * L + 0.0317567 * L * L + 0.078834 * L + 1))
    #
    current_name = {1:"Euler", 2:"Stokes"}
    #
    if is_finite:
        #print("# Height/Depth: {:6.3f}".format(MaxH))
        # print("# {:}".format(title))
        print("# ---------------------------------------------------------")
        print("# Height/Depth: {:7.4f}".format(z[2] / z[1]))
        print('# (A height of{:3.0f}% of the max of H/d={:1.3f} for this length)'
              .format(z[2] / z[1] / Highest * 100.0, Highest))        
        print("# Length/Depth: {:7.4f}".format(2 * math.pi / z[1]))
        print("# Dimensionless Period T*sqrt(g/d): {:7.4f}"
              .format(z[3] / math.sqrt(z[1])))
    else:
        print("# Height/Length: {:6.3f}, {:3.4}% of the maximum of H/L = {:6.3f}"
                .format(z[2]/2/pi,(z[2]/2/pi)/Highest*100., Highest))
        print("# Dimensionless Period T*sqrt(g/L): {:7.2f}".format(z[3]/math.sqrt(2*pi)))
    #
    print("# Current criterion: {:},  Dimensionless value: {:1.3f}"
          .format(current_name[c_type], current))    
#
def output(n, z, Y, B, Tanh, is_finite):
    """ """
    pi = math.pi
    kd = z[1]
    L=2*pi/z[1]
    H=z[2]/z[1]
    T=z[3]/math.sqrt(z[1])
    c=z[4]/math.sqrt(z[1])
    ce=z[5]/math.sqrt(z[1])
    cs=z[6]/math.sqrt(z[1])
    ubar=z[7]/math.sqrt(z[1])
    Q=ubar-z[8]/pow(kd,1.5)
    R=1+z[9]/z[1]
    
    if is_finite:
        pulse=z[8]+z[1]*z[5];
        ke=0.5*(z[4]*pulse-z[5]*Q*pow(kd,1.5))
    
        # Calculate potential energy, not by computing the mean of 1/2 (eta-d)^2
        # but by exploiting orthogonality of the cosine functions to give the sum of 1/4 Y[i]^2
        pe = 0
        for i in range(1,n):
            pe += 0.25*pow(Y[i], 2)
    
        ub2=2.*z[9]-z[4]*z[4]
        
        sxx=4.*ke-3.*pe+ub2*z[1]+2.*z[5]*(z[7]*z[1]-z[8])
        f=z[4]*(3.*ke-2.*pe)+0.5*ub2*(pulse+z[4]*z[1])+z[4]*z[5]*(z[7]*z[1]-z[8])
        q=z[7]*z[1]-z[8]
        r=z[9]+z[1]
        s=sxx-2.*z[4]*pulse+(z[4]*z[4]+0.5*z[1])*z[1]
    #
    print("# ---------------------------------------------------------")
    print("# Integral quantities - notation from Fenton (1988)")
    if is_finite:
        print("# (1) Quantity, (2) symbol, solution non-dimensionalised by")
        print("# (3) g & wavenumber, and (4) g & mean depth")
        print("# Water depth                        (d) {: 7.4f}  {: 7.4f}".format(z[1], 1.))
    else:
        print("# (1) Quantity, (2) symbol, solution non-dimensionalised by")
        print("# (3) g & wavenumber")
    #
    print("# Wave length                   (lambda) {: 7.4f}  {: 7.4f}".format(2*pi, L))
    print("# Wave height                        (H) {: 7.4f}  {: 7.4f}".format(z[2], H))
    print("# Wave period                      (tau) {: 7.4f}  {: 7.4f}".format(z[3], T))
    print("# Wave speed                         (c) {: 7.4f}  {: 7.4f}".format(z[4], c))
    print("# Eulerian current                 (u1_) {: 7.4f}  {: 7.4f}".format(z[5], ce))
    print("# Stokes current                   (u2_) {: 7.4f}  {: 7.4f}".format(z[6], cs))
    print("# Mean fluid speed in frame of wave (U_) {: 7.4f}  {: 7.4f}".format(z[7], ubar))
    print("# Volume flux due to waves           (q) {: 7.4f}  {: 7.4f}".format(z[8], z[8]/pow(kd,1.5)))
    print("# Bernoulli constant                 (r) {: 7.4f}  {: 7.4f}".format(z[9], z[9]/kd))
    #
    if is_finite:
        print("# Volume flux                        (Q) {: 7.4f}  {: 7.4f}".format(Q*pow(kd,1.5), Q))
        print("# Bernoulli constant                 (R) {: 7.4f}  {: 7.4f}".format(R*kd, R))
        print("# Momentum flux                      (S) {: 7.4f}  {: 7.4f}".format(s, s/kd/kd ))
        print("# Impulse                            (I) {: 7.4f}  {: 7.4f}".format(pulse, pulse/pow(kd,1.5)))
        print("# Kinetic energy                     (T) {: 7.4f}  {: 7.4f}".format(ke, ke/kd/kd))
        print("# Potential energy                   (V) {: 7.4f}  {: 7.4f}".format(pe, pe/kd/kd))
        print("# Mean square of bed velocity     (ub2_) {: 7.4f}  {: 7.4f}".format(ub2, ub2/kd))
        print("# Radiation stress                 (Sxx) {: 7.4f}  {: 7.4f}".format(sxx, sxx/kd/kd))
        print("# Wave power                         (F) {: 7.4f}  {: 7.4f}".format(f, f/pow(kd,2.5)))
    #
    print("# ---------------------------------------------------------")
    print("# Dimensionless coefficients in Fourier series" )
    print("# Potential/Streamfn  Surface elevations" )
    print("#    j          B[j]          E[j]  j=1..n" )
    #print("# N, number of dimensionless Fourier coefficients - j, B[j], & E[j] below", n)
    for i in range (1, n+1):
        print("{:6.0f} {: 1.6e} {: 1.6e}".format(i, B[i], Y[i]))
        #print("%2d\t%15.7e\t%15.7e\n", i, B[i], Y[i]);
    print("" )
    #    
#
def get_etas(n, z, Y, B, Tanh, nprofiles, is_finite):
    """
    Surface - print out coordinates of points on surface for plotting 
    plus check of pressure on surface.
    """
    pi = math.pi
    surface_points = nprofiles
    kd = z[1]
    c=z[4]/math.sqrt(z[1])
    ce=z[5]/math.sqrt(z[1])
    R=1+z[9]/z[1]
    # Surface - print out coordinates of points on surface for plotting 
    # plus check of pressure on surface
    #print("# %s\n", Title);
    #print("%s\n", Method);
    print("# Surface of wave - trough-crest-trough,")
    print("# note quadratic point spacing clustered around crest")
    if is_finite:
        print("# Non-dimensionalised with respect to depth")
        print("#    X/d   eta/d   check of surface pressure")
        #print("# Dummy point to scale plot")
    else:
        print("# Non-dimensionalised with respect to wavenumber")
        print("#    kX    k eta   check of surface pressure")    
    #
    s_range = surface_points // 2
    for i in range(-s_range, s_range+1):
        #NB Quadratic point spacing, clustered near crest
        X = 4 * pi * (i/surface_points)**2
        X = math.copysign(X, i)
        eta = Surface(X, Y, n)
        (y, Pressure, Bernoulli_check, 
         u, v, dphidt, ut, vt, ux, uy) = Point(X, eta, kd, Tanh, 
                                               B, n, ce, c, R, z, is_finite)
        if is_finite:
            print("{:8.4f} {:7.4f} {:7.0e}".format(X/kd, 1+eta/kd, Pressure))
        else:
            print("{:8.4f} {:7.4f} {:7.0e}".format(X, eta, Pressure))
    print("")
    #
    npt = number_steps(nprofiles)
    xx =  array('f', [0 for i in range(npt)])
    eta = array('f', [0 for i in range(npt)])    
    return xx, eta
#
def get_surface(n, z, Y, nprofiles):
    """ """
    pi = math.pi
    kd = z[1]
    #
    npt = number_steps(nprofiles)
    xx =  array('f', [0 for i in range(npt)])
    eta = array('f', [0 for i in range(npt)])
    for j in range(npt):
        #X = 0.5 * L * (j / nprofiles)
        X = pi * j / nprofiles
        eta[j] = Surface(X, Y, n)/kd
        xx[j] = X/kd
        #print("# X/d = {: 8.4f}, Phase = {: 6.1f} theta/d = {: 8.4f}"
        #      .format(xx[j], kd*xx[j] * 180 / pi, eta[j]))
    #print("")   
    return xx, eta 
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
    for j in range(nprofiles):
        #X = 0.5 * L * (j / nprofiles)
        xx[j] = pi * (j / nprofiles)
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
def get_kinematic(n, z, B, Tanh, d, etas, xx, points, is_finite):
    """ """
    g = 9.80665  # m/s^2
    pi = math.pi
    kd = z[1]
    c=z[4]/math.sqrt(z[1])
    ce=z[5]/math.sqrt(z[1])
    R=1+z[9]/z[1]
    #
    npt = len(etas)
    X = [xx[j] * kd for j in range(npt)]
    #eta = [etas[j] * kd for j in range(npt)]
    #
    if is_finite:
        # y = i * (1 + eta[j] / kd) / points
        y = [[i/points * (1 + etas[j]) for i in range(points+1)]
             for j in range(npt)]
        #
        #output = [[Point(X[j], kd*((i*y[j])-1), kd, Tanh,B, n, ce, c, R, z, is_finite) 
        #          for i in range(points)] for j in range(npt)]
        #
        #output = []
        #for j in range(npt):
        #    output.append([])
        #    for i in range(points):
        #        #print(y[j][i])
        #        output[j].append(Point(X[j], kd*(y[j][i]-1), 
        #                               kd, Tanh,B, n, ce, c, R, z, is_finite))             
        #
        factors = [d, (g*d)**0.5, (g*d)**0.5, g*d, g, g, (g/d)**0.5, (g/d)**0.5, g*d, 1]
        output = [['z', 'u', 'v', 'dphidt', 'ut', 'vt', 'ux', 'uz', 'pressure', 'Bernoulli_check']]
        #
        output.extend([Point(X[j], kd*(y[j][i]-1), 
                             kd, Tanh,B, n, ce, c, R, z, is_finite) 
                       for j in range(npt) for i in range(points+1)])
    else:
        #
        y = [[- pi + i / points * (etas[j]* kd + pi) for i in range(points+1)] 
              for j in range(npt)]
        #output = [[Point(X[j], y[j], kd, Tanh,B, n, ce, c, R, z, is_finite) 
        #          for i in range(points)] for j in range(npt)]
        #
        factors = [1/kd, (g/kd)**0.5, (g/kd)**0.5, g/kd, g, g, (g*kd)**0.5, (g*kd)**0.5, g/kd, 1]
        output = [['kz', 'u', 'v', 'dphidt', 'ut', 'vt', 'ux', 'uz', 'pressure', 'Bernoulli_check']]
        output.extend([Point(X[j], y[j][i], 
                             kd, Tanh,B, n, ce, c, R, z, is_finite) 
                       for j in range(npt) for i in range(points+1)] )       
    #
    output = list(zip(*output))
    #for j in range(npt):
    #    print("# X/d = {: 8.4f}, Phase = {: 6.1f}".format(X[j]/kd, X[j] * 180 / pi ))
    #    for i in range(points):
    #        print(*[f'{output[j][i][x]: 7.4f}' for x in range(9)], sep=' ')
    #
    #
    #for j in range(npt):
    #    X = xx[j] * kd
    #    eta = etas[j] * kd
    #    print("# X/d = {: 8.4f}, Phase = {: 6.1f}".format(X/kd, X * 180 / pi ))
    #    for i in range(points):
    #        if is_finite:
    #            y = i * (1 + eta / kd) / points
    #            (u, v, dphidt, ut, vt, ux, uy, 
    #             Pressure, Bernoulli_check) = Point(X, kd*(y-1), kd, Tanh,
    #                                                   B, n, ce, c, R, z, is_finite)                
    #        else:
    #            y = -pi + i / points * (eta + pi)
    #            (u, v, dphidt, ut, vt, ux, uy, 
    #             Pressure, Bernoulli_check) = Point(X, y, kd, Tanh, 
    #                                                  B, n, ce, c, R, z, is_finite)
    #
    #        print("{: 7.4f} {: 7.4f} {: 7.4f} {: 7.4f} {: 7.4f} {: 7.4f} {: 7.4f} {: 7.4f} {:7.0e}"
    #              .format(y, u, v, dphidt, ut, vt, ux, uy, Bernoulli_check))   
    #print('--')
    #for x, item in enumerate(header):
    #    output[x].insert(0, item)
    return output, factors
    
#
def print_velo(H, n, z, Y, B, Tanh, surface_points, method):
    """
    Surface velocity
    """
    pi = math.pi
    L = 2 * pi / z[ 1 ]
    # Highest wave - eqn (32) of Fenton (1990)
    Highest = ((0.0077829 * L * L * L + 0.0095721 * L * L + 0.141063 * L)
               / (0.0093407 * L * L * L + 0.0317567 * L * L + 0.078834 * L + 1))
    kd = z[1]
    c = z[4] / math.sqrt(z[1])
    ce = z[5] / math.sqrt(z[1])
    R = 1 + z[9] / z[1]
    X = 0.0
    eta = surface(X, kd, Y, n)
    Velo = [0]
    points = surface_points
    for i in range(0, points):
        y = i * eta / points
        Pressure, Bernoulli_check, u, v, dphidt, ut, vt, ux, uy = point(X, y, kd, Tanh, B, n, ce, c, R)
        # Velo[i] = u
        Velo.append(u)
    #
    i = 1
    sum1 = 0
    while i <= points - 1:
        sum1 += Velo[i]
        i += 2
    # Loop
    i = 2
    sum2 = 0
    while i <= points - 2:
        sum2 += Velo[i]
        i += 2
    # Loop
    ucm = (Velo[0] + 4 * sum1 + 2 * sum2 + Velo[points]) / 3.0 / points
    #
    print("# ---------------------------------------------------------")
    print("{:} {:4.0f} {:1.4e} {:1.4e} {:1.4e} {:1.4e} {:1.4e}"
          .format(method, n, H, L, 0.5 * z[2] / pow(z[1], 3),
                  z[2] / z[1] / Highest * 100.0, ucm))
    #
#
#  Surface elevation
def Surface(x, Y, n):
    """
    Surface elevation
    """
    kEta = 0
    kEta += sum([Y[j] * math.cos(j * x) 
                for j in range(1, n)])
    kEta += 0.5 * Y[n] * math.cos(n * x)
    return kEta
#
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
        Cos  = math.cos(j*X)
        Sin  = math.sin(j*X)
        if Is_finite:
            coshdelta = math.cosh(j*Y)
            sinhdelta = math.sinh(j*Y)
            C = coshdelta + sinhdelta*Tanh[j]
            S = sinhdelta + coshdelta*Tanh[j]
        else:
            #elif Is_deep:
            C = math.exp(j*Y)
            S = math.exp(j*Y)
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
    return [y, u, v, dphidt, ut, vt, ux, uy, Pressure, Bernoulli_check]
#
def get_Height(MaxH, case, is_finite, L=None, T=None):
    """ """
    case = case.lower()
    if is_finite:
        if case == "wavelength":
            Height = MaxH/L
        else:
            Height = MaxH/(T*T)
    else:
        Height = -MaxH
    
    return Height
#
def number_steps(StpLgth):
    """
    """
    npt = max(StpLgth, 2)
    if npt % 2 == 0:
        npt += 1
    return int(npt)
#
#
def Read_data(Extra):
    """ """
    fgets(Title, 400, Input1)
    Title = ChangeCharacter(Title, Title.Length - 1, NullChar)
    if Title == "FINISH":
        return 0
    # End if
    print(monitor, "# %s", Title)

    fscanf(Input1, "%""lf", MaxH)
    fgets(dummy, 100, Input1)
    # End if
    print(monitor, Lf & "# Height/Depth:%6.3f", MaxH)
    fscanf(Input1, "%s", [Case])
    fgets(dummy, 100, Input1)

    if Case == "Wavelength":
        fscanf(Input1, "%""lf", L)
        fgets(dummy, 100, Input1)
    # End if
    print("  Length/Depth:%7.2f", L)
    Height = MaxH / L
    # End if
    if Case == "Period":
        fscanf(Input1, "%""lf", T)
        fgets(dummy, 100, Input1)
        # End if
        Height = MaxH / (T * T)
        print("  Dimensionless Period T*sqrt(g/d):%7.2f", T)
        # End if
        fscanf(Input1, "%""d", Current_criterion)
        fgets(dummy, 100, Input1)
        # End if
        fscanf(Input1, "%""lf", Current)
        fgets(dummy, 100, Input1)
    # End if
    if Current_criterion == 1:
        Currentname = Current1
    # End if
    if Current_criterion == 2:
        Currentname = Current2
        # End if
        print("# Current criterion: %s,  Dimensionless value:%6.3lf", Currentname, Current)

        fscanf(Input1, "%""d", n)
        fgets(dummy, 100, Input1)
        # End if
        fscanf(Input1, "%""d", nstep)
        fgets(dummy, 100, Input1)
    #
    #  if wavelength is known at this stage the program calculates the highest wave,
    #  sees how close this wave is and allocates the number of steps automatically.
    if Case == "Wavelength":
        Highest = (0.0077829 * L * L * L + 0.0095721 * L * L + 0.141063 * L) / (
            0.0093407 * L * L * L + 0.0317567 * L * L + 0.078834 * L + 1)
        Console.Write(Lf & "{0:f}", MaxH / Highest)
        # nstep = 1+20*pow(MaxH/Highest,2);
    #
    if Extra == "n":
        return 1
    #
    #  Convergence criteria
    Input2 = fopen(Convergence_file, "r")
    fgets(dummy, 400, Input2)
    fscanf(Input2, "%d", number)
    fgets(dummy, 400, Input2)
    fscanf(Input2, "%le", crit)
    fgets(dummy, 400, Input2)
    fclose(Input2)

    #  Number of data points to present results for
    Input2 = fopen(Points_file, "r")
    fgets(dummy, 400, Input2)
    #  Number of points on surface profile (clustered quadratically near crest)
    fscanf(Input2, "%d", Surface_points)
    fgets(dummy, 400, Input2)
    #  Number of vertical profiles
    fscanf(Input2, "%d", Nprofiles)
    fgets(dummy, 400, Input2)
    #  Number of points in each profile
    fscanf(Input2, "%d", Points)
    fgets(dummy, 400, Input2)

    fclose(Input2)

    return 1
    # End Function
#
#

