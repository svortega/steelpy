# Copyright (c) 2019 steelpy


# Python stdlib imports
import math
from typing import Dict, List, Tuple

# package imports

#
# ************************************************
#                     PROGRAM RINGS
#              STRESSES IN CIRCULAR RINGS
#                R.J.ROARK & W.C.YOUNG
# ************************************************
#                    Log History
# Original Fortran Program (5th Ed) 27/12/10 - SVO
#
#
# ************************************************
#
#
#from buckingham import FindUnitLength, FindUnitForce
#import SectionPropertyModule as secProp

#
#
#
#-------------------------------------------------
#                 Roark Section
#-------------------------------------------------
#
def shear_factor(_c, _c1, _R, _e):
    '''
    Radial/horizontal web shear stress factor calculation
    Ref Roark 7 ed chapter 9 : Shear stress due to the radial 
                               shear force V
    
    tr : thickness of the section normal to the plane of curvature
         at the radial position r
    
    '''
    #
    #global _D, _Tw, _Bft, _Tft, _Bfb, _Tfb
    #
    #  
    # total area of the web 
    _Aw = (_D - _Tft - _Tfb) * _Tw       
    # Area of Flange
    _Af = (_Bft*_Tft) + (_Bfb*_Tfb)
    
    # Average Shear
    _tau_average = 1.0/_Aw
    print(' ')
    print('shear average : {: 1.4E}'.format(_tau_average))
    
    #
    # --------------------------------------------
    # Shear stress due to the radial shear force V
    _rn = _R -_e             # radial position to neautral axis
    _b = _R - _c - _Tft      # radial position to inner extreme fiber
    _tr = _Tw   # width of cross section where stresses are calculated
    
    # fraction of the web dep
    _dwo = (_c1 - _Tft - _e) / 3.0 
    _dwi = (-_c + _Tfb + _e) / 3.0
    _CoorZ = [(_c1 - _Tft), 2 * _dwo + _e, _dwo + _e, _e, 0,
              -_e, _dwi -_e, 2*_dwi -_e, (-_c + _Tfb)]
    
    _tau_radial = []
    for i in range(len(_CoorZ)):
        _r = _R + _CoorZ[i]      # radial position to point of interest
        _cr =  (_r - _b) / 2.0  # distance from centroide to extreme fibre
        _r1 = _cr + _b           # point of cross section from radial centre 
        _Ar = (_tr * math.log(_r /_b) )* _r  # area of part portion below r
        #_Ar = (_tr * math.log((_r1 / _cr + 1) / (_r1 / _cr - 1))) * _r
        _Qr = _Ar * _r1
        # roark Ed 7 equ (9.1-4)
        _tau_radial.append( (_rn / (_tr*_Aw*_e*_r**2)) * (_R*_Ar - _Qr))
    #
    #print(' ')
    #for i in range(len(_tau_radial)):
    #    print('tau : {: 1.4E} {: 1.4E}'
    #          .format(_tau_radial[i], _CoorZ[i]))
    #print(' ')
    print('shear radial  : {: 1.4E} '.format(max(_tau_radial)))
    
    _tau_y = max(_tau_average, max(_tau_radial))
    #print(' ')
    print('-----------------------------')
    print('Max Shear (tau) : {: 1.4E}'. format(_tau_y))
    
    return _tau_y
#
#-------------------------------------------------
#
#
def factors2(_E, _G, _nu, _A, _R, _I, _F, _e, _D):
    ''' 
    The Hoop-stress deformation factor : 
        alpha = I/AR^2
        
    The Transverse (radial) shear deformation factor : 
        beta = FEI/GAR^2  [thin rings]
        beta = 2F(1+v)e/R [thick rings]
            
    G is the shear modulus of elasticity
    F is the shape factor for the cross section [see Sec 8.10 (7ed)]
        
    Note that these constants are unity if not correction for hoop
    stress or shear stress is necessary or desired for use with thin
    rings
    
    '''
    if _R/_D < 0.60:
        raise RuntimeError('R/d [{:}] < 0.6 not applicable'.format(_R/_D ))
        #sys.exit('R/d < 0.6 not applicable')
    # thin ring
    if _R/_D > 8.0:
        # Hoop-stress deformation factor
        _alpha = _I / (_A * _R**2)     
        # Transverse (radial) shear deformation factor
        _beta = (_F * _E * _I) / (_G * _A * _R**2)
        _ring_status = 'thin ring'
    else: # thick ring
        # Hoop-stress deformation factor
        _alpha = _e / _R  
        # Transverse (radial) shear deformation factor
        _beta = (2 * _F * (1 + _nu) * _e) / _R   
        _ring_status = 'thick ring'
    # k constans
    k1 = 1 - _alpha + _beta
    k2 = 1 - _alpha
    #
    print('========= 0 ============')
    print(' k Factors')
    print ('R/d = [{:}] {:}'.format(_R/_D, _ring_status ))
    print('k1 =',k1, 'k2 =',k2)
    print('')
    return k1, k2
#
#
#
def solution(amom:float, ahop:float, ashr:float, 
             tm:float, tt:float, tv, istep:int, ist:int, crad:float, 
             sx:float, cx:float, j:int, 
             xmom:List[float], xhop:List[float], xshr:List[float]):
    '''
    General formulas for moment, hoop load and radial shear:
        M = MA - NAR(1-u) + VARz + LTM
        N = NAu + VAz + LTM
        V = -NAz + VAu + LTv
        
    Where LTM, LTN and LTv are loads terms given below for
    several types of loads.
    
    This funtion calculates the forces around the ring.
    only n/2 + 1 points are calculated around the ring
    for symmetric load cases 
    '''
    if abs(1000000.0 / math.pi - tv) < 0.0001:
        xmom[j] = amom
        xhop[j] = ahop
        xshr[j] = ashr
    else:
        xmom[j] = amom - ahop * crad * (1 - cx) + ashr * crad * sx + tm
        xhop[j] = ahop * cx + ashr * sx + tt
        xshr[j] = -ahop * sx + ashr * cx + tv
        

    if istep != ist:
        if j != 0 and j != ist - 1:
            k = istep - j 
            xmom[k] = xmom[j]
            xhop[k] = xhop[j]
            xshr[k] = -xshr[j]
    
    return xmom, xhop, xshr
#
#
#
def rotadd(istep:int, phase:float, 
           xmom:List[float], xhop:List[float], xshr:List[float]):
    """
    """
    _xloc = [0 for i in range(istep)]
    _xm = [0 for i in range(istep)]
    _xh = [0 for i in range(istep)]
    _xs = [0 for i in range(istep)]
    
    k = int(phase) * istep // 360
    #print('====> k :',k, phase, istep)
    step = 360.0 / float(istep)
    
    for i in range(istep):
        l = k - i 
        j = istep - l
        if l <= 0:
            j = abs(l) 
        #else:
        #    j = istep - l 
        #
        kk = i 
        akk = float(kk)
        _xloc[i] = akk * step
        _xm[i] += xmom[j]
        _xh[i] += xhop[j]
        _xs[i] += xshr[j]
    return _xloc, _xm, _xh, _xs
#
#
#
#
def stress2(smac, area, ceni, ceno, fki, fko, _shearFactor, istep, xloc, xm, xh, xs):
    '''
    Calculate stress for curved beams
    '''
    _cstmax = 0
    _step = 360.0 / float(istep)
    
    _bstri = [] # range(istep)    # Stresses Bending Inner
    _bstro = [] # range(istep)    # Stresses Bending Outer
    _hstr = [] # range(istep)     # Hoop Stresses
    _sstr = [] #range(istep)     # Shear Stresses
    _cstri = [] # range(istep)    # Stresses Bending Inner + Hoop
    _cstro = [] # range(istep)    # Stresses Bending Outer + Hoop
    
    for i in range(istep):
        # Circunferencial normal stress due tu pure bending M:
        #     Sigma = My/Aer or kiMc/Ix
        _bstri.append((fki  * xm[i] * ceni) / float(smac))
        _bstro.append(- ( fko * xm[i] * ceno)  / float(smac))
        
        # Circunferencial normal stress due to hoop tension N:
        #     Sigma = N/A
        _hstr.append(xh[i] / float(area))
        
        # Shear stress due to the radial shear force V:
        # Tau = 
        _sstr.append(xs[i] * _shearFactor)
        
        # Radial stress due to moment M and normal force N
        # Sigma = (My/Aer + N/A) * k
        _cstri.append(_bstri[i] + _hstr[i] * fki)
        _cstro.append(_bstro[i] + _hstr[i] * fko)
        
        k = i
        ak = float(k)
        xloc[i] = ak * _step
        
        x = abs(_cstri[i])
        y = abs(_cstro[i])
        
        if x > y:
            _temp = _cstri[i]
        else:
            _temp = _cstro[i]
        
        a = abs(_cstmax)
        b = abs(_temp)
        
        if a <= b:
            _cstmax = _temp
    
    
    _shmax = min(_sstr)
    if  max(_sstr) > abs(_shmax):
        _shmax = max(_sstr)
    
    #_CSTmax = max()
    
    #print('********************')
    #print ('_cstmax = {} _shmax = {}'.format(_cstmax, _shmax))
    return _cstmax, _shmax, _bstri, _bstro, _hstr, _sstr, _cstri, _cstro
#
#
#
#
#
#-------------------------------------------------
#                 Roark Cases
#-------------------------------------------------
#
# Diametrical opposed forces
def ring_1(crad, k1, k2, apld, phase, istep=24):
    '''
    Circular Ring Case 1 : 
    Diametri Point Loads
    
    Parameters
    ----------
    W = Load (force)
    
    Returns
    ----------
    M
    N
    V
    
    Notes
    ----------
    Formulas for stress and strain, 7th edition
    R.J.Roark & W.C.Young 
    
    Examples
    ----------
    '''
    
    global xmom, xhop, xshr
    
    xmom = [0 for i in range(istep)]
    xhop = [0 for i in range(istep)]
    xshr = [0 for i in range(istep)]
    
    print (' ')
    print ('roark case 1 diametri point loads')
    print ('')
    print('Applied load: ',apld)
    print('Phase: ',phase)
    print ('')
    #
    
    Z = 2 * math.pi
    amom = k2 * apld * crad / math.pi
    ahop = 0
    ashr = 0
    ist = (istep // 2) + 1
    
    for j in range(ist):
        i = j
        d = float(i) / float(istep)
        x = 2 * d * math.pi
        sx = math.sin(x)
        cx = math.cos(x)
        tm = (-apld * crad * sx) / 2.0
        tt = (-apld * sx) / 2.0
        tv = (-apld * cx) / 2.0
        
        solution(amom, ahop, ashr, tm, tt, tv, istep, ist, 
                 crad, sx, cx, j, i)
    
    xloc, xm, xh, xs = rotadd(istep, phase, xmom, xhop, xshr)
    
    return xloc, xm, xh, xs
    #
#
# Off centre opposed force
def ring_2(crad, k1, k2, apld, phase, theta, istep=24):
    '''
    Circular Ring Case 2 
    Point Loads at theta
    
    Parameters
    ----------
    W = Load (force)
    
    Returns
    ----------
    M
    N
    V
    
    Notes
    ----------
    Formulas for stress and strain, 7th edition
    R.J.Roark & W.C.Young 
    
    Examples
    ----------
    
    '''
    
    global xmom, xhop, xshr
    
    xmom = [0 for i in range(istep)]
    xhop = [0 for i in range(istep)]
    xshr = [0 for i in range(istep)]
    
    print (' ')
    print ('roark case 2 point loads at theta')
    print ('')
    print('Applied load: ',apld)
    print('Phase: ',phase)
    print('Theta: ',theta)
    print ('')
    
    Z = 2 * math.pi
    #a = theta * math.pi / 180.0
    a = math.radians(theta)
    sth = math.sin(a)
    cth = math.cos(a)
    #b = a / math.pi
    #c = sth / math.pi
    #amom = -apld * crad * ((1 - b) * (1 - cth) - c * (1 - cth * k2))
    amom = ((-apld * crad / math.pi) 
            * ((math.pi - a) * (1 - cth) -  sth * ( k2 - cth)))
    #ahop = -apld * (1 - b + c * cth)
    ahop = (-apld / math.pi)* (math.pi - a + sth * cth)
    ashr = 0
    ist = (istep // 2) + 1
    
    for j in range(ist):
        i = j
        d = float(i) / float(istep)
        x = 2 * d * math.pi
        sx = math.sin(x)
        cx = math.cos(x)
        
        if a > x:
            tm = 0
            tt = 0
            tv = 0
        
        else:
            tm = -apld * crad * (cth - cx)
            tt = apld * cx
            tv = -apld * sx
            
        solution(amom, ahop, ashr, tm, tt, tv, istep, ist, 
                 crad, sx, cx, j, i)
    
    xloc, xm, xh, xs = rotadd(istep, phase, xmom, xhop, xshr)
    
    return xloc, xm, xh, xs
    #
#
# Opposite moments
def ring_3(crad, k1, k2, apld, phase, theta, istep=24):
    '''
    Circular Ring Case 3 
    Moment Applied at theta
    
    Parameters
    ----------
    W = Load (force)
    
    Returns
    ----------
    M
    N
    V
    
    Notes
    ----------
    Formulas for stress and strain, 7th edition
    R.J.Roark & W.C.Young 
    
    Examples
    ---------- 
    ''' 
    
    global xmom, xhop, xshr
    
    xmom = [0 for i in range(istep)]
    xhop = [0 for i in range(istep)]
    xshr = [0 for i in range(istep)]
    
    print (' ')
    print ('roark case 3 moment applied at theta')
    print ("")
    print('Applied load: ',apld)
    print('Phase: ',phase)
    print('Theta: ',theta)
    print ("")
    
    a = theta * math.pi / 180.0
    sth = math.sin(a)
    cth = math.cos(a)
    #b = (2 * sth / (math.pi * fk1))
    #amom = -apld * (1 - (a / math.pi) - b)
    amom = ((-apld / math.pi) 
            * (math.pi - a - (2 * sth * k2 / k1)))
    #ahop = (apld / crad) * b
    ahop = ((apld / (crad * math.pi))
            * (2 * sth * k2 / k1))
    ashr = 0
    ist = (istep // 2) + 1
    
    for j in range(ist):
        i = j
        d = float(i) / float(istep)
        x = 2 * d * math.pi
        sx = math.sin(x)
        cx = math.cos(x)
        
        if (a > x):
            tm = 0
            tt = 0
            tv = 0
        
        else:
            tm = apld
            tt = 0
            tv = 0
        
        solution(amom, ahop, ashr, tm, tt, tv, istep, ist, 
                 crad, sx, cx, j, i)
    
    xloc, xm, xh, xs = rotadd(istep, phase, xmom, xhop, xshr)
    
    return xloc, xm, xh, xs
    #
#
# Two parallel forces with reaction
def ring_4(crad, k1, k2, apld, phase, theta, istep=24):
    '''
    Circular Ring Case 4 
    Vertical Point Loads at theta - Vertical 
    
    Parameters
    ----------
    W = Load (force)
    
    Returns
    ----------
    M
    N
    V
    
    Notes
    ----------
    Formulas for stress and strain, 7th edition
    R.J.Roark & W.C.Young 
    
    Examples
    ---------- 
    ''' 
    
    global xmom, xhop, xshr
    
    xmom = [0 for i in range(istep)]
    xhop = [0 for i in range(istep)]
    xshr = [0 for i in range(istep)]
    
    print (' ')
    print ('roark case 4 vertical point loads -vertical')
    print ("")
    print('Applied load: ',apld)
    print('Phase: ',phase)
    print('Theta: ',theta)
    print ("")
    
    a = theta * math.pi / 180.0
    sth = math.sin(a)
    cth = math.cos(a)
    #b = a / math.pi
    #c = sth * sth * fk4
    #p = 1.0 / math.pi
    #amom = -apld * crad * (p * (1 + cth + c) - sth * (1 - b))
    amom = ((-apld * crad / math.pi) 
            * (sth * (sth - math.pi + a) + k2 * (1 + cth)))
    #ahop = -apld * (c / math.pi)
    ahop = -apld * sth**2 / math.pi
    ashr = 0
    ist = (istep // 2) + 1
    
    for j in range(ist):
        i = j 
        d = float(i) / float(istep)
        x = 2 * d * math.pi
        sx = math.sin(x)
        cx = math.cos(x)
        ax = (x - 0.009)
        
        if (a > ax):
            tm = 0
            tt = 0
            tv = 0
        
        else:
            tm = apld * crad * (sx - sth)
            tt = apld * sx
            tv = apld * cx
        
        solution(amom, ahop, ashr, tm, tt, tv, istep, ist, 
                 crad, sx, cx, j, i, )
   
    xloc, xm, xh, xs = rotadd(istep, phase, xmom, xhop, xshr)
    
    return xloc, xm, xh, xs
    #
#
# Two radial forces with reaction
def ring_5(crad, k1, k2, apld, phase, theta, istep=24):
    '''
    Circular Ring Case 5 
    Normal Point Loads Applied at theta
    
    Parameters
    ----------
    W = Load (force)
    
    Returns
    ----------
    M
    N
    V
    
    Notes
    ----------
    Formulas for stress and strain, 7th edition
    R.J.Roark & W.C.Young 
    
    Examples
    ---------- 
    ''' 
    
    global xmom, xhop, xshr
    
    xmom = [0 for i in range(istep)]
    xhop = [0 for i in range(istep)]
    xshr = [0 for i in range(istep)]
    
    print (" ")
    print ('roark case 5 normal piont loads applied at theta')
    print (" ")
    print('Applied load: ',apld)
    print('Phase: ',phase)
    print('Theta: ',theta)
    print (" ")
    
    a = theta * math.pi / 180.0
    sth = math.sin(a)
    cth = math.cos(a)
    #b = a / math.pi
    #p = 1.0 / math.pi
    #amom = -apld * crad * (sth * (1 - b) - p * (1 + cth))
    amom = ((-apld * crad / math.pi) 
            * (sth * (math.pi - a) - k2 * (1 + cth)))
    #ahop = -apld * (sth * (1 - b))
    ahop = (-apld / math.pi) * sth*(math.pi - a)
    ashr = 0
    ist = (istep // 2) + 1
    
    for j in range(ist):
        i = j
        d = float(i) / float(istep)
        x = 2 * d * math.pi
        sx = math.sin(x)
        cx = math.cos(x)
        ax = (x - 0.009)
        
        if a > x:
            tm = 0
            tt = 0
            tv = 0
        
        else:
            tm = -apld * crad * math.sin(x - a)
            tt = -apld * math.sin(x - a)
            tv = -apld * math.cos(x - a)
        
        solution(amom, ahop, ashr, tm, tt, tv, istep, ist, 
                 crad, sx, cx, j, i)
    
    xloc, xm, xh, xs = rotadd(istep, phase, xmom, xhop, xshr)
    
    return xloc, xm, xh, xs
    #
#
# Two tangential forces with their reactions
def ring_6(crad, k1, k2, apld, phase, theta, istep=24):
    '''
    Circular Ring Case 6 
    Tangential Loads at theta
    
    Parameters
    ----------
    W = Load (force)
    
    Returns
    ----------
    M
    N
    V
    
    Notes
    ----------
    Formulas for stress and strain, 7th edition
    R.J.Roark & W.C.Young 
    
    Examples
    ''' 
    
    global xmom, xhop, xshr
    
    xmom = [0 for i in range(istep)]
    xhop = [0 for i in range(istep)]
    xshr = [0 for i in range(istep)]
    
    print (' ')
    print ('roark case 6 tangential loads at theta')
    print (" ")
    print('Applied load: ',apld)
    print('Phase: ',phase)
    print('Theta: ',theta)
    print (" ")
    
    a = theta * math.pi / 180.0
    sth = math.sin(a)
    cth = math.cos(a)
    b = a / math.pi
    c = sth / math.pi
    #amom = -apld * crad * (c * (1 + fk4) - (1 - b) * (1 - cth))
    amom = ((-apld * crad / math.pi) 
            * (sth * (1 + k2) - (math.pi - a) * (1 - cth)))
    #ahop = -apld * (c * fk4 + (1 - b) * cth)
    ahop = (-apld / math.pi) * (sth + (math.pi - a) * cth)
    ashr = 0
    ist = (istep // 2) + 1
    
    for j in range(ist):
        i = j
        d = float(i) / float(istep)
        x = 2 * d * math.pi
        sx = math.sin(x)
        cx = math.cos(x)
        #ax = (x - 0.009)
        
        if a > x :
            tm = 0
            tt = 0
            tv = 0
        
        else:
            tm = -apld * crad * (1 - math.cos(x - a))
            tt = apld * math.cos(x - a)
            tv = -apld * math.sin(x - a)
        
        solution(amom, ahop, ashr, tm, tt, tv, istep, ist, 
                 crad, sx, cx, j, i)
    
    xloc, xm, xh, xs = rotadd(istep, phase, xmom, xhop, xshr)
    
    return xloc, xm, xh, xs
    #
#
# Equal radial loads equally spaced
# THIS ROARK CASE IS NOT IMPLEMENTED
def ring_7(crad, k1, k2, apld, phase, theta, istep=24):
    #
    print ('***ERROR***THIS ROARK CASE IS NOT IMPLEMENTED')
    #
#
# Uniform unit load (at the base theta > pi/2)
def ring_8(crad, k1, k2, apld, phase, theta, istep=24):
    '''
    Circular Ring Case 8 
    UDL Applied at Base theta > pi/2
    
    Parameters
    ----------
    W = Load (force)
    
    Returns
    ----------
    M
    N
    V
    
    Notes
    ----------
    Formulas for stress and strain, 7th edition
    R.J.Roark & W.C.Young 
    
    Examples
    ''' 
    
    global xmom, xhop, xshr
    
    xmom = [0 for i in range(istep)]
    xhop = [0 for i in range(istep)]
    xshr = [0 for i in range(istep)]
    
    print (" ")
    print ('roark case 8 udl applied at base theta > pi/2')
    print ("")
    print('Applied load: ',apld)
    print('Phase: ',phase)
    print('Theta: ',theta)
    print ("")
    
    a = theta * math.pi / 180.0
    sth = math.sin(a)
    cth = math.cos(a)
    #c = sth * sth
    #d1 = sth * sth * sth
    #e = (sth + 0.75 * sth * cth + d1 * fk4 / 3 + 0.25 * a + a * c / 2)
    #p = 1.0 / math.pi
    #amom = apld * (crad * crad) * (0.25 + 0.5 * c - p * e)
    amom = ((apld * crad**2 / (2. * math.pi)) 
            * (math.pi * (sth**2 - 0.50) - ((sth * cth - a) / 2.) 
               - (sth**2 * (a + 2. * sth / 3.)) 
               - k2 * (2 * sth + sth * cth - math.pi + a)))
    #ahop = -apld * crad * (d1 * fk4 / (3 * math.pi))
    ahop = -apld * crad * sth**3 / (3.0 * math.pi)
    ashr = 0
    ist = (istep // 2) + 1
    
    for j in range(ist):
        i = j
        d = float(i) / float(istep)
        x = 2 * d * math.pi
        sx = math.sin(x)
        cx = math.cos(x)
        
        if a > x:
            tm = 0
            tt = 0
            tv = 0
        
        else:
            # Bug : Undefined variable CRAC replaced with the correct variable
            #       crad in the following statement by djd on 25/01/90
            tm = (-apld * crad**2 / 2.0) * (sx - sth)** 2
            tt = -apld * crad * sx * (sx - sth)
            tv = -apld * crad * cx * (sx - sth)
        
        solution(amom, ahop, ashr, tm, tt, tv, istep, ist, 
                 crad, sx, cx, j, i)
    
    xloc, xm, xh, xs = rotadd(istep, phase, xmom, xhop, xshr)
    
    return xloc, xm, xh, xs
    #
#
# Linear varying unit load (at the base theta > pi/2)
def ring_9(crad, k1, k2, apld, phase, theta, istep=24):
    '''
    Circular Ring Case 9 
    TDL at Base theta > pi/2
    
    Parameters
    ----------
    W = Load (force)
    
    Returns
    ----------
    M
    N
    V
    
    Notes
    ----------
    Formulas for stress and strain, 7th edition
    R.J.Roark & W.C.Young 
    
    Examples
    ''' 
    
    global xmom, xhop, xshr
    
    xmom = [0 for i in range(istep)]
    xhop = [0 for i in range(istep)]
    xshr = [0 for i in range(istep)]
    
    print ("")
    print ('roark case 9 tdl at base theta > pi/2')
    print ("")
    print('Applied load: ',apld)
    print('Phase: ',phase)
    print('Theta: ',theta)
    print ("")
    
    a = theta * math.pi / 180.0
    sth = math.sin(a)
    cth = math.cos(a)
    #c = sth * sth * sth
    #d1 = apld * crad * crad
    #e = a / math.pi
    #p = 1 / math.pi
    #f = p * (11 * sth * cth / 36.0 + (1 + cth) / (9.0 * sth) + 0.5 * sth)
    #g = c * k2 / (12.0 * math.pi)
    #amom = d1 * ((1 - e) * (0.25 + (sth * sth) / 6.0) - f - g)
    amom = ((apld * crad**2 / (36. * math.pi * sth))
            * ((math.pi - a) * (6. * sth**3 - 9. * sth) - (3 * sth**4)
               + 8 + (8 * cth) - (5 * sth**2 * cth) 
               - (6 * k2 * (3 * sth * (sth - math.pi + a)
                             + (sth**2 * cth) + 2 + (2 * cth)))))
    #ahop = -apld * crad * g
    ahop = -apld * crad * sth**3 / (12. * math.pi)
    ashr = 0
    ist = (istep // 2) + 1
    
    for j in range(ist):
        i = j
        d = float(i) / float(istep)
        x = 2 * d * math.pi
        sx = math.sin(x)
        cx = math.cos(x)
        
        if a > x:
            tm = 0
            tt = 0
            tv = 0
        
        else:
            #tm = (d1 / 6) * ((sx * sx * sx) / sth + 3 * sx * sth - sth * sth - 3 * sx * sx)
            tm = (apld * crad**2 / (6. * sth)) * (sx - sth)** 3
            #tt = ((apld * crad * sx) / sth) * ((sx * sx) / sth + sth - 2 * sx)
            tt = ((apld * crad * sx) / sth**2) * (sx - sth)** 2
            #tv = ((apld * crad * cx) / 2) * ((sx * sx) / sth + sth - 2 * sx)
            tv = ((apld * crad * cx) / (2. * sth)) * (sx - sth)** 2
        
        solution(amom, ahop, ashr, tm, tt, tv, istep, ist, 
                 crad, sx, cx, j, i)
    
    xloc, xm, xh, xs = rotadd(istep, phase, xmom, xhop, xshr)
    
    return xloc, xm, xh, xs
    #
#
# Uniform compressive unit load (from 180 to theta)
def ring_10(crad, k1, k2, apld, phase, theta, istep=24):
    '''
    Circular Ring Case 10 
    UDL from 180 to theta
    
    Parameters
    ----------
    W = Load (force)
    
    Returns
    ----------
    M
    N
    V
    
    Notes
    ----------
    Formulas for stress and strain, 7th edition
    R.J.Roark & W.C.Young 
    
    Examples
    ''' 
    
    global xmom, xhop, xshr
    
    xmom = [0 for i in range(istep)]
    xhop = [0 for i in range(istep)]
    xshr = [0 for i in range(istep)]
    
    print (" ")
    print ('roark case 10 u.d.l. from 180 to theta')
    print ("")
    print('Applied load: ',apld)
    print('Phase: ',phase)
    print('Theta: ',theta)
    print ("")
    
    a = theta * math.pi / 180.0
    sth = math.sin(a)
    cth = math.cos(a)
    #b = (sth**3) * k2 / 3.0
    #c = sth * sth / 2.0
    #d = a * cth
    #e = apld * crad
    #p = 1.0 / math.pi
    #amom = -e * crad * (cth + c - 0.75 + p * (sth - d + 0.75 * (a - sth * cth) - a * c - b))
    amom = ((-apld * crad**2 / (4. * math.pi))
            * ((math.pi - a)*(4 * cth + 2 * sth**2 - 1.0)
               + sth * (4 - (4 * sth**2 / 3.0) - cth)
               - 2 * k2 * (math.pi - a + sth * cth)))
    #ahop = -e * (cth + p * (sth - d - b))
    ahop = ((-apld * crad / math.pi) 
            * (math.pi * cth + sth - a * cth - sth**3 / 3.0))
    ashr = 0
    ist = (istep // 2) + 1
    
    for j in range(ist):
        i = j
        d = float(i) / float(istep)
        x = 2 * d * math.pi
        sx = math.sin(x)
        cx = math.cos(x)
        
        if a > x:
            tm = 0
            tt = 0
            tv = 0
        
        else:
            #tm = (-e * crad / 2.0) * (1 + cx * cx - sth * sth - 2 * cx * cth)
            tm = (-apld * crad**2 / 2.0) * (cth - cx)**2
            #tt = e * cx * (cth - cx)
            tt = (apld * crad * cx) * (cth - cx)
            #tv = -e * sx * (cth - cx)
            tv = (-apld * crad * sx) * (cth - cx)
        
        solution(amom, ahop, ashr, tm, tt, tv, istep, ist, 
                 crad, sx, cx, j, i)
    
    xloc, xm, xh, xs = rotadd(istep, phase, xmom, xhop, xshr)
    
    return xloc, xm, xh, xs
    #
#
# Linear varying unit load (from 180 to theta)
def ring_11(crad, k1, k2, apld, phase, theta, istep=24):
    '''
    Circular Ring Case 11 
    TDL from 180 to theta
    
    Parameters
    ----------
    W = Load (force)
    
    Returns
    ----------
    M
    N
    V
    
    Notes
    ----------
    Formulas for stress and strain, 7th edition
    R.J.Roark & W.C.Young 
    
    Examples
    ''' 
    
    global xmom, xhop, xshr
    
    xmom = [0 for i in range(istep)]
    xhop = [0 for i in range(istep)]
    xshr = [0 for i in range(istep)]
    
    print (' ')
    print ('roark case 11 tdl form 180 to theta')
    print ("")
    print('Applied load: ',apld)
    print('Phase: ',phase)
    print('Theta: ',theta)
    print ("")
    
    a = theta * math.pi / 180.0
    sth = math.sin(a)
    cth = math.cos(a)
    #b = a / math.pi
    #s2 = sth * sth
    #s3 = sth * sth * sth
    #c2 = cth * cth
    #c3 = cth * cth * cth
    #f = apld * crad * crad / (1 + cth)
    #g = 1 - b
    #h = (0.75 * cth - 0.5 * c2 - c3 / 3.0 - s2 * cth / 2.0 - 0.25)
    #q = (1.0 / math.pi) * (sth / 9.0 - (3 * sth * cth) / 4.0 + (11 * sth * c2) / 36.0)
    #r = k2 / (8 * math.pi)
    #r1 = r * (math.pi - a + sth * cth + (2 * s3 * cth) / 3.0)
    #amom = f * (g * h + q + r1)
    amom = ((-apld * crad**2 / (math.pi * (1 + cth)))
            * ((math.pi - a) * (((3 + 12 * cth**2 + 2 *cth + 4 * cth * sth**2) / 24.0)
                                - ((3 * sth**3 * cth - 3 * sth - 5 * sth**3) / 36.0)
                                + (5 * sth * cth / 8.0)
                                - k2 * (math.pi * cth / 2.0 - a * cth / 2.0 
                                         + sth**3 / 3.0 + sth * cth**2 / 2.0))))
    #w = (apld * crad) / (1 + cth)
    #p = 1 / math.pi
    #ahop = w * (p * (a / 4.0 + a * c2 / 2.0 - 3 * sth * cth / 4.0) - 0.25 - c2 / 2.0 + r1)
    ahop = ((-apld * crad / (math.pi * (1 + cth)))
            * ((math.pi - a) * (((1 + 4 * cth**2) / 8.0) + (5 * sth * cth / 8.0)
                                - (sth**3 * cth / 12.0))))
    ashr = 0
    ist = (istep // 2) + 1
    
    for j in range(ist):
        i = j
        d = float(i) / float(istep)
        x = 2 * d * math.pi
        sx = math.sin(x)
        cx = math.cos(x)
        
        if a > x:
            tm = 0
            tt = 0
            tv = 0
        
        else:
            #tm = f * ((cx**3) / 6.0 - c3 / 6.0 + cx * c2 / 2.0 - cx * cx * cth / 2.0)
            tm = (-apld * crad**2 / (6.0 * (1 + cth))) * (cth - cx)**3
            #tt = ((apld * crad * cx) / (2.0 * (1 + cth))) * ((cth - cx)**2)
            tt = (apld * crad * cx / (2.0 * (1 + cth))) * (cth - cx)**2
            #tv = ((apld * crad * sx) / (2.0 * (1 + cth))) * ((cth - cx)**2)
            tv = (-apld * crad * sx / (2.0 * (1 + cth))) * (cth - cx)**2
        
        solution(amom, ahop, ashr, tm, tt, tv, istep, ist, 
                 crad, sx, cx, j, i)
    
    xloc, xm, xh, xs = rotadd(istep, phase, xmom, xhop, xshr)
    
    return xloc, xm, xh, xs
    #
#
# Uniform radial unit load (from 180 to theta)
def ring_12(crad:float, k1:float, k2:float, apld:float, 
            phase:float, theta:float, istep:int=24):
    '''
    Circular Ring Case 12
    Normal UDL Base to theta
    
    Parameters
    ----------
    W = Load (force)
    
    Returns
    ----------
    M
    N
    V
    
    Notes
    ----------
    Formulas for stress and strain, 7th edition
    R.J.Roark & W.C.Young 
    
    Examples 
    ''' 
    xmom = [0 for i in range(istep)]
    xhop = [0 for i in range(istep)]
    xshr = [0 for i in range(istep)]
    
    print ('')
    print ('roark case 12 normal udl base to theta')
    print ('')
    print('Applied load: ',apld)
    print('Phase: ',phase)
    print('Theta: ',theta)
    print ('')
    
    ist = (istep // 2) + 1
    a = theta * math.pi / 180.0
    sth = math.sin(a)
    cth = math.cos(a)
    #p = 1.0 / math.pi
    #amom2 = (-apld * crad**2 * 
    #        (p * (a + 2 * sth - a * cth) - 1 + cth))
    amom = (-apld * crad**2 / math.pi
            * (sth + math.pi * cth - a * cth 
               - k2 * (math.pi - a - sth)))
    #ahop = -apld * crad * (p * (sth - a * cth) + cth)
    ahop = (-apld * crad / math.pi
            * (sth + math.pi * cth - a * cth))
    ashr = 0
    
    for j in range(ist):
        i = j
        d = i / istep
        x = 2 * d * math.pi
        sx = math.sin(x)
        cx = math.cos(x)
        
        if a > x:
            tm = 0
            tt = 0
            tv = 0
        else:
            tm = -apld * crad**2 * (1 - math.cos(x - a))
            tt = -apld * crad * (1 - math.cos(x - a))
            tv = -apld * crad * math.sin(x - a)
        
        xmom, xhop, xshr = solution(amom, ahop, ashr, tm, tt, tv, istep, ist, 
                                    crad, sx, cx, j, xmom, xhop, xshr)
    
    xloc, xm, xh, xs = rotadd(istep, phase, xmom, xhop, xshr)
    return xloc, xm, xh, xs
#
#
# Lineary varying radial unit load (from 180 to theta)
def ring_13(crad, k1, k2, apld, phase, theta, istep=24):
    '''
    Circular Ring Case 13
    Lineary varying radial unit load (from 180 to theta)
    
    Parameters
    ----------
    W = Load (force)
    
    Returns
    ----------
    M
    N
    V
    
    Notes
    ----------
    Formulas for stress and strain, 7th edition
    R.J.Roark & W.C.Young 
    
    Examples 
    ''' 
    global xmom, xhop, xshr
    
    xmom = [0 for i in range(istep)]
    xhop = [0 for i in range(istep)]
    xshr = [0 for i in range(istep)]
    
    print ('')
    print ('roark case 13 Lineary varying radial unit load')
    print ('')
    print('Applied load: ',apld)
    print('Phase: ',phase)
    print('Theta: ',theta)
    print ('')
    
    ist = (istep // 2) + 1
    a = theta * math.pi / 180.0
    sth = math.sin(a)
    cth = math.cos(a)
    
    amom = ((-apld * crad**2 / (math.pi * (math.pi - a)))
            * (2 + 2 * cth - sth * (math.pi - a) 
               + k2 * (1 + cth - (math.pi - a)**2 / 2.0)))
    
    ahop = ((-apld * crad / (math.pi * (math.pi - a)))
            * (2 + 2 * cth - sth * (math.pi - a)))
    
    ashr = 0
    
    for j in range(ist):
        i = j
        d = float(i) / float(istep)
        x = 2 * d * math.pi
        sx = math.sin(x)
        cx = math.cos(x)
        
        if a > x:
            tm = 0
            tt = 0
            tv = 0
        
        else:
            tm = ((-apld * crad **2 / (math.pi - a)) 
                  * (x - a - sx * cth + cx * sth))
            tt = ((-apld * crad / (math.pi - a)) 
                  * (x - a - sx * cth + cx * sth))
            tv = ((-apld * crad / (math.pi - a)) 
                  * (1 - cx * cth - sx * sth))
        
        solution(amom, ahop, ashr, tm, tt, tv, istep, ist, 
                 crad, sx, cx, j, i)
    
    
    xloc, xm, xh, xs = rotadd(istep, phase, xmom, xhop, xshr)
    
    return xloc, xm, xh, xs
    #
#
# Radial quadratic unit load (from 180 to theta)
def ring_14(crad, k1, k2, apld, phase, theta, istep=24):
    '''
    Circular Ring Case 14
    Radial quadratic unit load (from 180 to theta)
    
    Parameters
    ----------
    W = Load (force)
    
    Returns
    ----------
    M
    N
    V
    
    Notes
    ----------
    Formulas for stress and strain, 7th edition
    R.J.Roark & W.C.Young 
    
    Examples 
    ''' 
    global xmom, xhop, xshr
    
    xmom = [0 for i in range(istep)]
    xhop = [0 for i in range(istep)]
    xshr = [0 for i in range(istep)]
    
    print ('')
    print ('roark case 14 Radial quadratic unit load')
    print ('')
    print('Applied load: ',apld)
    print('Phase: ',phase)
    print('Theta: ',theta)
    print ('')
    
    ist = (istep // 2) + 1
    a = theta * math.pi / 180.0
    sth = math.sin(a)
    cth = math.cos(a)
    
    amom = ((-apld * crad**2 / (math.pi * (math.pi - a)**2))
            * (2 * (math.pi - a) * (2 - cth) - 6 * sth
               + k2 * (2 * (math.pi - a - sth) 
                        - (math.pi - a)**3 / 3.0)))
    
    ahop = ((-apld * crad / (math.pi * (math.pi - a)**2))
            * (2 * (math.pi - a) * (2 - cth) - 6 * sth))
    
    ashr = 0
    
    for j in range(ist):
        i = j
        d = float(i) / float(istep)
        x = 2 * d * math.pi
        sx = math.sin(x)
        cx = math.cos(x)
        
        if a > x:
            tm = 0
            tt = 0
            tv = 0
        
        else:
            tm = ((-apld * crad **2 / (math.pi - a)**2) 
                  * ((x - a)**2 - 2 + 2 * cx *cth 
                     + 2 * sx * sth))
            tt = ((-apld * crad / (math.pi - a)**2) 
                  * ((x - a)**2 - 2 + 2 * cx *cth
                     + 2 * sx * sth))
            tv = ((-2 * apld * crad / (math.pi - a)**2) 
                  * (x - a - sx * cth + cx * sth))
        
        solution(amom, ahop, ashr, tm, tt, tv, istep, ist, 
                 crad, sx, cx, j, i)
    
    
    xloc, xm, xh, xs = rotadd(istep, phase, xmom, xhop, xshr)
    
    return xloc, xm, xh, xs
    #
#
# Ring supported at the base and loaded by own weight
def ring_15(crad, k1, k2, apld, phase, theta, 
            _Area, _I, _beta, istep=24):
    '''
    Circular Ring Case 15
    Self Weight Supported at Base
    
    Parameters
    ----------
    W = Load (force)
    
    Returns
    ----------
    M
    N
    V
    
    Notes
    ----------
    Formulas for stress and strain, 7th edition
    R.J.Roark & W.C.Young 
    
    Examples 
    ''' 
    
    global xmom, xhop, xshr
    
    xmom = [0 for i in range(istep)]
    xhop = [0 for i in range(istep)]
    xshr = [0 for i in range(istep)]
    
    print (' ')
    print ('roark case 15 self weight supported at base')
    print ("")
    print('Applied load: ',apld)
    print('Phase: ',phase)
    print('Theta: ',theta)
    print ("")
    
    a = theta * math.pi / 180.0
    sth = math.sin(a)
    cth = math.cos(a)
    #amom = apld * (crad * crad) * k2 / 2.0
    _Kt = 1 + (_I / (_Area * crad**2))
    amom = ((apld * crad**2) 
            * (k2 - 0.50 - ((_Kt - 1) * _beta / k1)))
    #ahop = apld * crad * k2 / 2.0
    ahop = apld * crad * (0.50 + ((_Kt - 1) * k2 / k1))
    ashr = 0
    ist = (istep // 2) + 1
    
    for j in range(ist):
        i = j
        d = float(i) / float(istep)
        x = 2 * d * math.pi
        sx = math.sin(x)
        cx = math.cos(x)
        #tm = -apld * (crad * crad) * (x * sx + cx - 1)
        tm = (-apld * crad**2) * (x * sx + _Kt * (cx - 1))
        tt = -apld * crad * x * sx
        tv = -apld * crad * x * cx
        
        solution(amom, ahop, ashr, tm, tt, tv, istep, ist, 
                 crad, sx, cx, j, i)
    
    xloc, xm, xh, xs = rotadd(istep, phase, xmom, xhop, xshr)
    
    return xloc, xm, xh, xs
    #
#
# Unit axial segment of pipe filled with liquid
def ring_16(crad, k1, k2, apld, phase, theta, _rho, istep=24):
    '''
    Circular Ring Case 16
    Unit axial segment of pipe filled with liquid
    
    Parameters
    ----------
    W = Load (force)
    
    Returns
    ----------
    M
    N
    V
    
    Notes
    ----------
    Formulas for stress and strain, 7th edition
    R.J.Roark & W.C.Young 
    
    Examples 
    ''' 
    
    global xmom, xhop, xshr
    
    xmom = [0 for i in range(istep)]
    xhop = [0 for i in range(istep)]
    xshr = [0 for i in range(istep)]
    
    print (' ')
    print ('roark case 16 Unit axial segment of pipe filled with liquid')
    print ("")
    print('Applied load: ',apld)
    print('Phase: ',phase)
    print('Theta: ',theta)
    print ("")
    
    a = theta * math.pi / 180.0
    sth = math.sin(a)
    cth = math.cos(a)
    amom = ((_rho * crad**3) * (0.75 - k2 / 2.0))
    ahop = _rho * crad**2 * 0.75
    ashr = 0
    ist = (istep // 2) + 1
    
    for j in range(ist):
        i = j
        d = float(i) / float(istep)
        x = 2 * d * math.pi
        sx = math.sin(x)
        cx = math.cos(x)
        
        tm = _rho * crad**3 * (1 - cx - x * sx / 2.0)
        tt = _rho * crad**2 * (1 - cx - x * sx / 2.0)
        tv = _rho * crad**2 * (sx / 2.0 - x * cx / 2.0)
        
        solution(amom, ahop, ashr, tm, tt, tv, istep, ist, 
                 crad, sx, cx, j, i)
    
    xloc, xm, xh, xs = rotadd(istep, phase, xmom, xhop, xshr)
    
    return xloc, xm, xh, xs
    #
#
# Unit axial segment of pipe partly filled with liquid
def ring_17(crad, k1, k2, apld, phase, theta, _rho, istep=24):
    '''
    Circular Ring Case 17
    Unit axial segment of pipe partly filled with liquid
    
    Parameters
    ----------
    W = Load (force)
    
    Returns
    ----------
    M
    N
    V
    
    Notes
    ----------
    Formulas for stress and strain, 7th edition
    R.J.Roark & W.C.Young 
    
    Examples 
    ''' 
    
    global xmom, xhop, xshr
    
    xmom = [0 for i in range(istep)]
    xhop = [0 for i in range(istep)]
    xshr = [0 for i in range(istep)]
    
    print (' ')
    print ('roark case 17 Unit axial segment of pipe partly filled with liquid')
    print ("")
    print('Applied load: ',apld)
    print('Phase: ',phase)
    print('Theta: ',theta)
    print ("")
    
    a = theta * math.pi / 180.0
    sth = math.sin(a)
    cth = math.cos(a)
    amom = ((_rho * crad**3 / (4. * math.pi)) 
            * (2 * a * sth**2 + 3 * sth * cth - 3 * a + math.pi
               + 2 * math.pi * cth**2
               + 2 * k2 * (sth * cth - 2 * sth 
                            + (math.pi - a) * (1 - 2 * cth))))
    ahop = ((_rho * crad**2 / (4. * math.pi)) 
            * (3 * sth * cth + (math.pi - a) * (1 + 2 * cth**2)))
    ashr = 0
    ist = (istep // 2) + 1
    
    for j in range(ist):
        i = j
        d = float(i) / float(istep)
        x = 2 * d * math.pi
        sx = math.sin(x)
        cx = math.cos(x)
        
        if a > x:
            tm = 0
            tt = 0
            tv = 0
        
        else:
            tm = ((_rho * crad**3 / 2.0) 
                  * (2 * cth - sx * (x - a + sth * cth)
                     - cx * (1 + cth**2)))
            tt = ((_rho * crad**2 / 2.0) 
                  * (2 * cth - sx * (x - a + sth * cth)
                     - cx * (1 + cth**2)))
            tv = ((_rho * crad**2 / 2.0) 
                  * (sx * cth**2 - cx * (x - a + sth * cth)))
        
        solution(amom, ahop, ashr, tm, tt, tv, istep, ist, 
                 crad, sx, cx, j, i)
    
    xloc, xm, xh, xs = rotadd(istep, phase, xmom, xhop, xshr)
    
    return xloc, xm, xh, xs
    #
#
# Two antiparallel forces offset by angle phi
def ring_18(crad, k1, k2, apld, phase, theta, thi, istep=24):
    '''
    Circular Ring Case 18 
    Loads Applied at theta & phi + Shear
    
    Parameters
    ----------
    W = Load (force)
    
    Returns
    ----------
    M
    N
    V
    
    Notes
    ----------
    Formulas for stress and strain, 7th edition
    R.J.Roark & W.C.Young 
    
    Examples 
    ''' 
    #global xmom, xhop, xshr
    xmom = [0 for i in range(istep)]
    xhop = [0 for i in range(istep)]
    xshr = [0 for i in range(istep)]
    
    print (' ')
    print ('roark case 18 loads applied at theta & phi + Shear')
    print ("")
    print('Applied load: ',apld)
    print('Phase: ',phase)
    print('Theta: ',theta)
    print('Phi: ',thi)
    print ("")
    
    thi = thi * math.pi / 180.0
    a = theta * math.pi / 180.0
    sth = math.sin(a)
    cth = math.cos(a)
    sthi = math.sin(thi)
    cthi = math.cos(thi)
    #b = crad / (2.0 * math.pi)
    #c = sth * sth
    #d = sthi * sthi
    #amom = (-apld * b * (cth - cthi - (math.pi - a) *
    #                     sth + (math.pi - thi) * sthi + k2 * (c - d)))
    amom = ((apld * crad / (2.0 * math.pi)) 
            * (sthi**2 - sth**2 - (math.pi - thi) * sthi 
               + (math.pi - a) * sth - k2 * (cth - cthi)))
    #        
    #ahop = ((-apld / (2.0 * math.pi)) * (k2 * (c - d)))
    ahop = ((apld / (2.0 * math.pi)) * (sthi**2 - sth**2))
    #
    #ashr = ((apld / (2.0 * math.pi)) * 
    #        (a - thi + sth - sthi + k2 * 
    #         (sth * cth - sthi * cthi)))
    ashr = ((apld / (2.0 * math.pi)) 
            * (a - thi + sth - sthi + sth * cth - sthi * cthi))
    
    ist = (istep // 2) + 1
    
    for j in range(istep):
        i = j
        d = i / istep
        x = 2 * d * math.pi
        sx = math.sin(x)
        cx = math.cos(x)
        Z = 1
        z1 = 1
        ax = x - 0.009
        
        if thi > ax : 
            Z = 0
        
        if a > ax : 
            z1 = 0
        
        tm = (((-apld * crad / (2.0 * math.pi)) 
               * (sthi - sth) * (x - sx))
              + apld * crad * (sx - sth) * z1 
              - apld * crad * (sx - sthi) * Z)
        
        tt = (((apld / (2.0 * math.pi)) 
               * (sthi - sth) * sx) 
              + apld * sx * z1 
              - apld * sx * Z)
        
        tv = (((-apld / (2.0 * math.pi)) 
               * (sthi - sth) * (1 - cx))
              + apld * cx * z1 
              - apld * cx * Z)
        
        xmom, xhop, xshr = solution(amom, ahop, ashr, tm, tt, tv, istep, ist, 
                                    crad, sx, cx, j, xmom, xhop, xshr)
    
    xloc, xm, xh, xs = rotadd(istep, phase, xmom, xhop, xshr)
    return xloc, xm, xh, xs
#
#
# Moment on the rim
def ring_19(crad, k1, k2, apld, phase, theta, istep=24):
    '''
    Circular Ring Case 19 
    Ring Shear + Moment
    
    Parameters
    ----------
    W = Load (force)
    
    Returns
    ----------
    M
    N
    V
    
    Notes
    ----------
    Formulas for stress and strain, 7th edition
    R.J.Roark & W.C.Young 
    
    Examples 
    ''' 
    
    global xmom, xhop, xshr
    
    xmom = [0 for i in range(istep)]
    xhop = [0 for i in range(istep)]
    xshr = [0 for i in range(istep)]
    
    print (' ')
    print ('roark case 19 ring shear + moment')
    print ("")
    print('Applied load: ',apld)
    print('Phase: ',phase)
    print('Theta: ',theta)
    print ("")
    
    a = theta * math.pi / 180.0
    sth = math.sin(a)
    cth = math.cos(a)
    #p = 1.0 / math.pi
    #amom = (-apld * p / 2.0) * (math.pi - a - 2 * k2 * sth)
    amom = ((-apld / (2.0 * math.pi)) 
            * (math.pi - a - (2 * k2 * sth / k1)))
    #ahop = (apld * p / crad) * k2 * sth
    ahop = ((apld / (crad * math.pi)) 
            * (k2 * sth / k1))
    #ashr = (-apld * p / (2.0 * crad)) * (1 + 2 * k2 * cth)
    ashr = ((-apld / (2.0 * crad * math.pi)) 
            * (1 + (2 * k2 * cth / k1)))
    
    ist = (istep // 2) + 1
    
    for j in range(istep):
        i = j
        d = float(i) / float(istep)
        x = 2 * d * math.pi
        sx = math.sin(x)
        cx = math.cos(x)
        #z1 = 1
        
        #if a > x :
        #    Z = 0
        
        #tm = (-apld * p / 2.0) * (x - sx) + apld * z1
        tm = (((-apld / (2.0 * math.pi)) * (x - sx))
              + apld)
        #tt = apld * sx * p / (2.0 * crad)
        tt = apld * sx / (2.0 * math.pi * crad)
        #tv = (-apld * p / (2.0 * crad)) * (1 - cx)
        tv = (-apld / (2.0 * math.pi * crad)) * (1 - cx)
        
        solution(amom, ahop, ashr, tm, tt, tv, istep, ist, 
                 crad, sx, cx, j, i)
    
    xloc, xm, xh, xs = rotadd(istep, phase, xmom, xhop, xshr)
    
    return xloc, xm, xh, xs
    #
#
# Bulkhead or supporting ring in pipe
def ring_20(crad, k1, k2, apld, phase, theta, istep=24):
    '''
    Circular Ring Case 18 
    Point Support at Base + Tangential Shear
    
    Parameters
    ----------
    W = Load (force)
    
    Returns
    ----------
    M
    N
    V
    
    Notes
    ----------
    Formulas for stress and strain, 7th edition
    R.J.Roark & W.C.Young 
    
    Examples 
    ''' 
    global xmom, xhop, xshr
    
    xmom = [0 for i in range(istep)]
    xhop = [0 for i in range(istep)]
    xshr = [0 for i in range(istep)]
    
    print (' ')
    print ('roark case 20 point support at base + tangential shear')
    print ("")
    print('Applied load: ',apld)
    print('Phase: ',phase)
    print('Theta: ',theta)
    print ("")
    
    a = theta * math.pi / 180
    sth = math.sin(a)
    cth = math.cos(a)
    p = 1 / math.pi
    #amom = (apld * crad * p / 2.0) * (k2 - 0.5)
    amom = (apld * crad / (2.0 * math.pi)) * (k2 - 0.5)
    #ahop = (apld * p / 2.0) * (k2 + 0.5)
    ahop = (0.75 * apld / math.pi)
    ashr = 0
    ist = (istep // 2) + 1
    
    for j in range(ist):
        i = j
        d = float(i) / float(istep)
        x = 2 * d * math.pi
        sx = math.sin(x)
        cx = math.cos(x)
        #tm = (apld * crad * p) * (1 - cx - x * sx * 0.5)
        tm = (apld * crad / math.pi) * (1 - cx - (x * sx / 2.0))
        #tt = (-apld * p / 2.0) * x * sx
        tt = (-apld / (2.0 * math.pi))* x * sx
        #tv = (apld * p / 2.0) * (sx - x * cx)
        tv = (apld / (2.0 * math.pi)) * (sx - x * cx)
        
        solution(amom, ahop, ashr, tm, tt, tv, istep, ist, 
                 crad, sx, cx, j, i)
        
    xloc, xm, xh, xs = rotadd(istep, phase, xmom, xhop, xshr)
    
    return xloc, xm, xh, xs
    #
#
#
def ring_19a(crad, fk4, apld, phase, theta, istep=24):
    '''
    Circular ring case 19 
    self weight supported by tangential shear
    Formulas for stress and strain, 5th edition
    R.J.Roark & W.C.Young  
    ''' 
    global xmom, xhop, xshr
    
    xmom = [0 for i in range(istep)]
    xhop = [0 for i in range(istep)]
    xshr = [0 for i in range(istep)]
    
    print (' ')
    print ('roark case 19 self weight supported by tangential shear ')
    print ("")
    print('Applied load: ',apld)
    print('Phase: ',phase)
    print('Theta: ',theta)
    print ("")
    
    a = theta * math.pi / 180.0
    sth = math.sin(a)
    cth = math.cos(a)
    p = 1 / math.pi
    amom = (apld * (crad * crad) / 2.0) * (1 - fk4)
    ahop = (-apld * crad / 2.0) * (1 + fk4)
    ashr = 0
    ist = (istep // 2) + 1
    
    for j in range(ist):
        i = j
        d = float(i) / float(istep)
        x = 2 * d * math.pi
        sx = math.sin(x)
        cx = math.cos(x)
        tm = -apld * crad * crad * (1 - cx)
        tt = 0
        tv = -apld * crad * sx
        
        solution(amom, ahop, ashr, tm, tt, tv, istep, ist, 
                 crad, sx, cx, j, i)
    
    xloc, xm, xh, xs = rotadd(istep, phase, xmom, xhop, xshr)
    
    return xloc, xm, xh, xs
    #
#
#
def ring_25(crad, fk4, apld, phase, theta, istep=24):
    '''
    Circular Ring Case 25 (14) 
    Point Loads at theta + Shear
    Formulas for stress and strain, 4th edition
    R.J.Roark & W.C.Young  
    ''' 
    
    global xmom, xhop, xshr
    
    xmom = [0 for i in range(istep)]
    xhop = [0 for i in range(istep)]
    xshr = [0 for i in range(istep)]
    
    print (' ')
    print ('roark case 25 (14) point loads at theta + shear')
    print ("")
    print('Applied load: ',apld)
    print('Phase: ',phase)
    print('Theta: ',theta)
    print ("")
    
    # input applied load multiplied by 2 as roark 4th edition gives w/2
    apld = 2 * apld
    a = theta * math.pi / 180.0
    sth = math.sin(a)
    cth = math.cos(a)
    ist = (istep // 2) + 1
    
    for j in range(ist):
        i = j
        d = float(i) / float(istep)
        x = 2 * d * math.pi
        sx = math.sin(x)
        cx = math.cos(x)
        
        if x <= a :
            
            xmom[j] = (apld * crad * (0.23868 * cx - sth / 2.0 + 0.15915 * 
                                      (x * sx + a * sth + cth - cx * cth * cth)))
            
            xhop[j]  = (apld * (0.15915 * (x * sx - cx * cth * cth) - 0.07958 * cx))
            
            xshr[j]  = (apld * (0.15915 * (x * cx - sx / 2.0 + sx * cth * cth)))
                                
        else:
            xmom[j]  = (apld * crad * (0.23868 * cx - sx / 2.0 + 0.15915 * 
                                      (x * sx + a * sth + cth - cx * cth * cth)))
            
            xhop[j]  = (apld * (0.15915 * (x * sx - cx * cth * cth) - 0.07958 * cx - sx / 2.0))
            
            xshr[j]  = (apld * (0.15915 * (x * cx - sx / 2.0 + sx * cth * cth) - cx / 2.0))
        
        #if (j != 1 or j != ist):
        if j != 0 and j != (ist-1):
            #k = istep + 2 - j
            k = istep  - j 
            xmom[k] = xmom[j]
            xhop[k] = xhop[j]
            xshr[k] = -xshr[j]
    
    xloc, xm, xh, xs = rotadd(istep, phase, xmom, xhop, xshr)
    
    return xloc, xm, xh, xs
    #
#
#
def ring_81(crad, fk4, apld, phase, theta, thi, istep=24):
    '''
    Circular ring case 8 (15) 
    loads applied at theta & phi
    Formulas for stress and strain, 4th edition
    R.J.Roark & W.C.Young  
    ''' 
    
    global xmom, xhop, xshr
    
    xmom = [0 for i in range(istep)]
    xhop = [0 for i in range(istep)]
    xshr = [0 for i in range(istep)]
    
    print (' ')
    print ('roark case 8 (15) loads applied at theta & phi')
    print ("")
    print('Applied load: ',apld)
    print('Phase: ',phase)
    print('Theta: ',theta)
    print('Thi: ',thi)
    print ("")
    
    thi = thi * math.pi / 180.0
    a = theta * math.pi / 180.0
    sth = math.sin(a)
    cth = math.cos(a)
    sthi = math.sin(thi)
    cthi = math.cos(thi)
    ist = (istep // 2) + 1
    
    for j in range(ist):
        i = j
        d = float(i) / float(istep)
        x = 2 * d * math.pi
        sx = math.sin(x)
        cx = math.cos(x)
        
        if x <= a :
            xmom[j] = (apld * crad * (0.3183 * 
                                      (sthi * thi + cthi -sth * a - 
                                       cth - cx * sth * sth + 
                                       cx * sthi * sthi) -sthi + sth))
            
            xhop[j] = (apld * (0.3183 * cx * (sthi * sthi - sth * sth)))
            
            xshr[j] = (apld * (0.3183 * sx * (sth * sth - sthi * sthi)))
        
        elif x <= thi:
            xmom[j] = (apld * crad * (0.3183 * (sthi * thi + cthi - sth * 
                                                a- cth - cx * sth * sth +cx * 
                                                sthi * sthi) - sthi + sx))
            
            xhop[j] = (apld * (0.3183 * cx * (sthi * sthi - sth * sth) + sx))
            
            xshr[j] = (apld * (0.3183 * sx * (sth * sth - sthi * sthi) + cx))
        
        else:
            xmom[j] = (apld * crad * (0.3183 * (sthi * thi + cthi - sth * 
                                                a - cth - cx * sth * sth + 
                                                cx * sthi * sthi)))
            
            xhop[j] = (apld * (0.3183 * cx * (sthi * sthi - sth * sth)))
            
            xshr[j] = (apld * (0.3183 * sx * (sth * sth - sthi * sthi)))
        
        #if (j != 1 or j != ist):
        if j != 0 and j != (ist-1):
            #k = istep + 2 - j
            k = istep  - j 
            xmom[k] = xmom[j]
            xhop[k] = xhop[j]
            xshr[k] = -xshr[j]
    
    xloc, xm, xh, xs = rotadd(istep, phase, xmom, xhop, xshr)
    
    return xloc, xm, xh, xs
    #
#
#
def circular_ring(ring_case:int, ring_load:float, 
                  theta:float, phi:float, 
                  phase:float, istep:int,
                  R:float, k1:float, k2:float,):
    """
    theta:
    """
    #global k1, k2, _alpha, _beta
    #global fk1, fk2, fk3, fk4, fk5
    #global _R, _Area, _Ic
    #global smac,area,warea,ceni,ceno,fki,fko,crad
    #---------------------
    # start of main subroutineme loop - each load case is called 
    # in turn control is passed to roark routines 
    #
    #  insert extra roark cases here
    # k1 = 1, k2 = 1,
    #
    #
    print ('==> roark case number: ',ring_case)
    # Updated to the 7th edition
    # Case 1
    if ring_case == 1:
        _xloc, _xm, _xh, _xs =  ring_1(R, k1, k2, 
                                       ring_load, phase, istep)
    
    # Case 2
    elif ring_case == 2:
        _xloc, _xm, _xh, _xs = ring_2(R, k1, k2, 
                                      ring_load, phase, theta, istep)
    
    # Case 3
    elif ring_case == 3:
        _xloc, _xm, _xh, _xs = ring_3(R, k1, k2, 
                                      ring_load, phase, theta, istep)
    
    # Case 4
    elif ring_case == 4:
        _xloc, _xm, _xh, _xs = ring_4(R, k1, k2, 
                                      ring_load, phase, theta, istep)
    
    # Case 5
    elif ring_case == 5:
        _xloc, _xm, _xh, _xs = ring_5(R, k1, k2, 
                                      ring_load, phase, theta, istep)
    
    # Case 6
    elif ring_case == 6:
        _xloc, _xm, _xh, _xs = ring_6(R,  k1, k2, 
                                      ring_load, phase, theta, istep)
    
    # Case 7
    #***error***this roark case is not implemented
    elif ring_case == 7:
        _xloc, _xm, _xh, _xs = ring_7(R, k1, k2, 
                                      ring_load, phase, theta, istep)
    
    # Case 8
    elif ring_case == 8:
        _xloc, _xm, _xh, _xs = ring_8(R, k1, k2, 
                                      ring_load, phase, theta, istep)
     
    # Case 9
    elif ring_case == 9:
        _xloc, _xm, _xh, _xs = ring_9(R, k1, k2, 
                                      ring_load, phase, theta, istep)
    
    # Case 10
    elif ring_case == 10:
        _xloc, _xm, _xh, _xs = ring_10(R, k1, k2, 
                                       ring_load, phase, theta, istep)
    
    # Case 11
    elif ring_case == 11:
        _xloc, _xm, _xh, _xs = ring_11(R, k1, k2,  
                                       ring_load, phase, theta, istep)
    
    # Case 12
    elif ring_case == 12:
        _xloc, _xm, _xh, _xs = ring_12(R, k1, k2, 
                                       ring_load, phase, theta, istep)
    
    # check it
    # Case 13
    elif ring_case == 13:
        _xloc, _xm, _xh, _xs = ring_13(R, k1, k2, 
                                       ring_load, phase, theta, istep)
    
    # check it
    # Case 14
    elif ring_case == 14:
        _xloc, _xm, _xh, _xs = ring_14(R, k1, k2, 
                                       ring_load, phase, theta, istep)
    #
    # Case 15
    elif ring_case == 15:
        _xloc, _xm, _xh, _xs = ring_15(R, k1, k2, 
                                       ring_load, phase, 
                                       theta, _Area, _Ic, _beta, istep)
     
    # check it beta & gamma factors to be fixed
    # Case 16
    elif ring_case == 16:
        _rho = 0.01
        _xloc, _xm, _xh, _xs = ring_16(R, k1, k2, 
                                       ring_load, phase, 
                                       theta, _rho, istep)
    # check it beta & gamma factors to be fixed
    # Case 17
    elif ring_case == 17:
        _rho = 0.01
        _xloc, _xm, _xh, _xs = ring_17(R, k1, k2, 
                                       ring_load, phase, 
                                       theta, _rho, istep)
     
    #
    # Case 18
    elif ring_case == 18:
        _xloc, _xm, _xh, _xs = ring_18(R, k1, k2, 
                                       ring_load, phase, theta, phi, istep)
    
    # Case 19
    elif ring_case == 19:
        _xloc, _xm, _xh, _xs = ring_19(R, k1, k2, 
                                       ring_load, phase, theta, istep)
    
    # Case 20
    elif ring_case == 20:
        _xloc, _xm, _xh, _xs = ring_20(R, k1, k2, 
                                       ring_load, phase, theta, istep)
    
    #
    #------ ?? ------
    #
    # Case 14a
    elif ring_case == 144:
        _xloc, _xm, _xh, _xs = ring_25(R, k2, 
                                       ring_load, phase, theta, istep)
    
    # Case 15a
    elif ring_case == 155:
        _xloc, _xm, _xh, _xs = ring_81(R, fk4, 
                                       ring_load, phase, theta, phi, istep)
    # Case 19a
    elif ring_case == 199:
        _xloc, _xm, _xh, _xs = ring_19a(R, fk4, 
                                        ring_load, phase, theta, istep)
    
    #
    else:
        print('***ERROR***THIS ROARK CASE IS NOT IMPLEMENTED')
        print('This case will not be included in calculations')
    #
    #
    #
    #print ('LOCATION    MOMENT      THRUST      SHEAR')
    #print ('             N.mm         N           N')
    
    #for i in range (len(_xloc)):
    #    print("{:>5.1f}    {: 1.4e} {: 1.4e} {: 1.4e}"
    #          .format(_xloc[i], _xm[i], _xh[i], _xs[i]))
    #
    #
    print (' ')
    # end of main subroutineme loop
    return _xloc, _xm, _xh, _xs
    #
#
#
def print_results():
    """
    """
    #global wlh, wth, fow, fot, fiw, fit
    #global _orad, _fow, _fot, _fiw, _fit, _wlh, _wth
    
    print (" ")
    print ("        STRESSES IN CIRCULAR RINGS 5th edition")
    print ("                R.J.ROARK & W.C.YOUNG")
    print ("******************************************************")
    print (" ")
    print ('RINGS DATA ECHO')
    #print ('number of load cases: {}'.format(iload))
    #print ('number of points around ring: {}'.format(istep))
    print ('outside radius of ring: {}'.format(_orad))
    print ('outside flange: {} x {} thk'.format(_fow, _fot))
    print ('inside flange:  {} x {} thk'.format(_fiw, _fit))
    print ('web             {} x {} thk'.format(_wlh, _wth))
    print ("")
    
    _Type = 'I SECTION'
    _D = _wlh + _fot + _fit
    _Tw = _wth
    _Bft = _fow
    _Tft = _fot
    _Bfb = _fiw
    _Tfb = _fit
    _Units = 'mm'
    
    print(" ")
    print("_______________________________________________________________________________________")
    print(" ")
    print("                            SECTION PROPERTIES REPORT ")
    print("                                  [30/07/2011]")
    print(" ")
    print("                              TUBULAR              WEB(i)               FLANGE(i)")
    print("Section Type            Diameter   Thickness   Depth   Thickness    Width    Thickness")
    print(" ")
    print(".......................................................................................")
    print(" ")
    
    if _Type == 'TUBULAR':
        print("{:1.4E} {:1.4E}".format(_D, _Tw))
    
    elif _Type == 'I SECTION':
        print("{:23s} {:>19} {:1.4E} {:1.4E} {:1.4E} {:1.4E}"
              .format(_Type, "", _D, _Tw, _Bft, _Tft))
        print("{:>65} {:1.4E} {:1.4E}".format("",_Bfb, _Tfb))
    
    elif _Type == 'TEE':
        print("{:23s} {:>19} {:1.4E} {:1.4E} {:1.4E} {:1.4E}"
              .format(_Type, "", _D, _Tw, _Bft, _Tft))
    
    elif _Type == 'CHANNEL':
        print("{:23s} {:>19} {:1.4E} {:1.4E} {:1.4E} {:1.4E}"
              .format(_Type, "", _D, _Tw, _Bft, _Tft))
        print("{:>65} {:1.4E} {:1.4E}".format("",_Bft, _Tft))
    
    elif _Type == 'BOX':
        print("{:23s} {:>19} {:1.4E} {:1.4E} {:1.4E} {:1.4E}"
              .format(_Type, "", _D, _Tw, _Bft, _Tft))
        print("{:>43} {:1.4E} {:1.4E} {:1.4E} {:1.4E}"
              .format("",_D, _Tw, _Bft, _Tft))
    
    elif _Type == 'RECTANGULAR BAR':
        print("{:>20} {:1.4E} {:1.4E}".format("",_D, _Bft))
    
    elif _Type == 'ANGLE':
        print("{:23s} {:>19} {:1.4E} {:1.4E} {:1.4E} {:1.4E}"
              .format(_Type, "", _D, _Tw, _Bft, _Tft))
    #
    
    if _Units == 'mm':
        test2 = {'unitlength' : 'mm', 'unitmass' :"kg/m", 'mayoraxis' : "y", 'minoraxis' : "x"}
    #
    print(" ")
    print("_______________________________________________________________________________________")
    print(" ")
    print("                               SECTION DERIVED PROPERTIES                UNITS [  {:4s}]".format(_Units))
    print(" ")
    
    #
#
