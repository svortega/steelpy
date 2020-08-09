#
# *******************************************
#                 Bug History
#
# DNV Code development     SVO     01/03/2011
# Checked against Mathcad R1 A.Aal 06/12/2011
#
# *******************************************
#
from __future__ import print_function
import math
import datetime
#
#
#
def PlateStiffener (b, t, hw, tw, bf, tf, L, fyp, fyw = 0, fyf = 0, E = 210000):
    #
    # plate (flange) section
    _b = b
    _t = t
    # web section
    _hw = hw
    _tw = tw
    # flange section
    _bf = bf
    _tf = tf
    # Fy plate
    _Fyp = fyp
    # Fy web
    if fyw == 0: _Fyw = _Fyp
    else: _Fyw = fyw
    # fy flange
    if fyf == 0 : _Fyf = _Fyp
    else: _Fyf = fyf
    # section length
    _L = L
    # Elatic modulus
    _E = E
    #
    # Cross Sectional Area
    # area plate (flange)
    _Ap = _b * _t
    # area web
    _Aw = _hw * _tw
    # area flange
    _Af = _bf * _tf
    # cross sectional area of stiffener
    _As = _Aw + _Af
    # total area
    _A = _Ap + _Aw + _Af
    #
    #
    _ef = _tw / 2.0
    #
    #
    # Equivalent Yield Strength over
    # the cross-section
    _Fyeq = (_Ap * _Fyp + _Aw * _Fyw + _Af * _Fyf) / _A
    #
    # Distance from outer surface of attached
    # plating to elastic horizontal neutral axis
    #
    _Zp = (((0.50 * _b * _t**2) + _Aw * (_t + 0.50 * _hw) 
            + _Af * (_t + _hw + 0.50 * _tf)) / _A)
    #
    _Ztf = (t + _hw + _tf - _Zp)
    #
    #print('Z0 ==> ',_Zp, _Ztf )
    #
    # Moment of Inertia
    #
    _I = ((_b * _t**3 / 12.0) + (_Ap * (_Zp - _t / 2.0)**2) 
          + (_hw**3 * _tw / 12.0) + _Aw * (_Zp - _t - _hw / 2.0)**2 
          + (_bf * _tf**3 / 12.0) + _Af * (_t + _hw + _tf / 2.0 - _Zp)**2)
    #
    #print('Ixx =',_I)
    #
    # Section Modulus
    #
    _Wep = _I / _Zp
    #
    _Wes = _I / _Ztf
    #
    _W = min(_Wep, _Wes)
    #
    #print ('--->', _Wep, _Wes, _W)
    #
    # Radius of gyration
    #
    _r = math.sqrt(_I / _A)
    #
    #
    # Column Slenderness Ratio
    #
    _Lambda = (_L / (math.pi * _r) * math.sqrt(_Fyeq / _E))
    #
    #
    # Torsional Moment of Inertia
    _Iz = (_Af * _bf**2 / 12.0) + (_ef**2 * _Af / float(1 + (_Af / _Aw)))
    #
    #print('Iz ===>',_Iz)
    #
    print ('Af = ',_Af)
    print ('Aw = ',_Aw)
    print ('Iz = ',_Iz)
    print ('ef = ',_ef)
    print ('bf = ',_bf)
    # Plate Slenderness Ratio
    #
    try:
	_Beta = _b / _t * math.sqrt(_Fyp / _E)
    #
    except :
	_Beta = 0
    #
    #
    return _A, _As, _Zp, _I, _Iz, _Wes, _Wep, _Fyeq
#
#
def LoadCase(LoadDescription):
    #
    _LoadDescDummy = ''
    #
    if 'LONGITUDINAL' in LoadDescription:
	#
	_LoadDescDummy = 'LONGITUDINAL ' + str(_LoadDescDummy)
    #
    #
    if 'TRANSVERSE' in LoadDescription:
	#
	_LoadDescDummy = 'TRANSVERSE ' + str(_LoadDescDummy)
    #
    #
    if 'SHEAR' in LoadDescription:
	#
	_LoadDescDummy = 'SHEAR ' + str(_LoadDescDummy)
    #
    #
    if 'LATERAL-PRESSURE' in LoadDescription:
	#
	_LoadDescDummy = 'LATERAL-PRESSURE ' + str(_LoadDescDummy)
    #
    return _LoadDescDummy
    #
    #
#
#
def PlateCase(PlateDescription):
    #
    _PlateDescDummy = ''
    #
    if 'PLATE' in PlateDescription:
	#
	#
	if 'STIFFENER' in PlateDescription:
	    #
	    if 'GIRDER' in PlateDescription:
		#
		#print('Girder Supporting Stiffened Panel')
		_PlateDescDummy  = 'GIRDER STIFFENED PLATE'
	    #
	    else:
		#
		#print('Longitudinal Stiffened Pannel')
		_PlateDescDummy = 'STIFFENED PLATE'
	#
	#
	else:
	    #
	    #print('Unstiffened Plate')
	    _PlateDescDummy = 'UNSTIFFENED PLATE'
	    #
    #
    else:
	#
	self.msg = self.msg + '\n' + 'not defined yet'
    #
    #
    return _PlateDescDummy
    #
#
#
def StressDefinition(Sigma1Sd, Sigma2Sd, L, L1):
    #
    #
    if Sigma2Sd == 0 :
	#
	Sigma2Sd = Sigma1Sd
	_SigmaSd = Sigma1Sd
	_SigmaType = 'UNIFORM'
	#
    #
    elif Sigma1Sd == Sigma2Sd :
	#
	# _SigmaSd = Sigmay1Sd
	_SigmaSd = Sigma1Sd
	_SigmaType = 'UNIFORM'
	#
    #
    else:
	# Average SigmaSd
	_SigmaSd = (max( (Sigma2Sd 
                          + (Sigma1Sd - Sigma2Sd) 
                          * (L - L1) / L ) 
                         , 0.75 * Sigma1Sd ))
	#
	_SigmaType = 'VARYING'
	#
    #
    #
    return Sigma2Sd, _SigmaSd, _SigmaType 
    #
#
#
def StressSign(SigmaxSd, SigmaySd):
    #
    _Stress = 'TENSION'
    #
    if SigmaxSd and SigmaySd != 0:
	#
	_SignX = SigmaxSd / abs(SigmaxSd)
	_SignY = SigmaySd / abs(SigmaySd)
	#
	if _SignX == 1.0 and _SignY == 1.0:
	    #
	    _Stress = 'COMPRESSION'
	    #
    #
    else:
	#
	if SigmaxSd !=  0:
	    #
	    _SignX = SigmaxSd / abs(SigmaxSd)
	    #
	    if _SignX == 1.0:
		#
		_Stress = 'COMPRESSION'
		#
	    #
	#
	if SigmaySd !=  0:
	    #
	    _SignY = SigmaySd / abs(SigmaySd)
	    #
	    if _SignY == 1.0:
		#
		_Stress = 'COMPRESSION'
		#
	    #
	#
    #
    return _Stress
    #
#
#
#
#
class SectionProperty():
    #
    #
    def BeamSection (self):
	#
	_s = (self.s1 + self.s2)/2
	# Area
	self.Apf = _s * self.t
	# print ('Area :', self.Apf)
	#
	#
	# Cross sectional area of stiffener
	_Aw = self.hw*self.tw
	_Af = self.b * self.tf
	_As = _Aw + _Af
	# print ('Cross sectional area of stiffener :',_As)
	#
	# Total Area
	self.Asf = self.Apf + _Aw + _Af
	# print ('Total Area :', self.Asf )
	#
	#
	# Neutral Axis Location
	self.Ztf = ((self.Apf*((0.5*self.t + self.hw + self.tf) + 
                               (_Aw * (0.5*self.hw + self.tf) + 
                                0.5*self.tf*_Af))) / self.Asf)
	#
	# print ('Ztf:', self.Ztf )
	#
	_zpf = 0.5*self.t + self.hw + self.tf - self.Ztf
	# print ('_zpf:', _zpf )
	#
	# Moment of Inertia
	self.Ipf = (((self.Apf * self.t**2)/12) + 
                    self.Apf * (0.5*self.t + self.hw + self.tf - self.Ztf)**2)
	# print ('Ipf:', self.Ipf )
	#
	self.Iwf = (((self.tw * self.hw**3)/12) + 
                    _Aw * (0.5*self.hw + self.tf - self.Ztf)**2)
	# print ('Iwf :', self.Iwf )
	#
	self.Iff = (((self.b * self.tf**3)/12) + 
                    _Af * (0.5*self.tf - self.Ztf)**2)
	# print ('Iff :', self.Iff )
	#
	self.Is = self.Ipf + self.Iwf + self.Iff 
	#
	# print ('Moment of Inertia:', self.Is )
	#
	# 
	#
    #
    #
    def PlateStiffener (self):
	#
	# Cross Sectional Area
	#
	_Ap = _b * _t
	#
	_Ape = _be * _t
	#
	_Aw = _hw * _tw
	#
	_Af = _bf * _tf
	#
	_A = _Ap + _Aw + _Af
	#
	_Ae = _Ape + _Aw + _Af
	#
	# Equivalent Yield Strength over
	# the cross-section
	_Fyeq = (_Ap * _Fyp + _Aw * _Fyw + _Af * _Fyf) / _A
	#
	# Distance from outer surface of attached
	# plating to elastic horizontal neutral axis
	#
	_Z0 = (((0.50 * _b * _t**2) + _Aw * (_t + 0.50 * _hw) 
                + _Af * (_t + _hw + 0.50 * _tf)) / _A)
	#
	#
	_Zp = (((0.50 * _be * _t**2) + _Aw * (_t + 0.50 * _hw) 
                + _Af * (_t + _hw + 0.50 * _tf)) / _Ae)
	#
	# Moment of Inertia
	#
	_I = ((_b * _t**3 / 12.0) + (_Ap * (_Z0 - _t / 2.0)**2) 
              + (_hw**3 * _tw / 12.0) + _Aw * (_Z0 - _t - _hw / 2.0)**2 
              + (_bf * _tf**3 / 12.0) + _Af * (_t + _hw + _tf / 2.0 - Z0)**2)
	#
	#
	_Ie = ((_be * _t**3 / 12.0) + (_Ape * (_Zp - _t / 2.0)**2) 
               + (_hw**3 * _tw / 12.0) + _Aw * (_Zp - _t - _hw / 2.0)**2 
               + (_bf * _tf**3 / 12.0) + _Af * (_t + _hw + _tf / 2.0 - Zp)**2)
	#
	# Radius of gyration
	#
	_r = math.sqrt(_I / _A)
	#
	_re = math.sqrt(_Ie / _A)
	#
	# Column Slenderness Ratio
	#
	_Lambda = (_L / (math.pi * _r) * math.sqrt(_Fyeq / _E))
	#
	_Lambdae = (_L / (math.pi * _re) * math.sqrt(_Fyeq / _E))
	#
	# Plate Slenderness Ratio
	#
	_Beta = _b / _t * math.sqrt(_Fyp / _E)
	#
	#
    #
    #
#
#
#
class DNVRPC201():
    #
    #
    def __init__(self):
	#
	#self.ur = 0.0
	self.ur_combined = 0.0
	self.ur_combined_Plate = 0.0
	self.ur_combined_Girder = 0.0
	self.ur_lateral = 0.0
	self.ur_longcomp = 0.0
	self.ur_transComp = 0.0
	self.ur_biaxial = 0.0
	self.msg = ''
	self._Taucrl = 0.0
	#

    # Torsional Buckling Strength
    def fT(self,_LT, _Beta, _Iz):
	#
	# area web
	_Aw = self.hw * self.tw
	# area flange
	_Af = self.bf * self.tf
	#
	# For L- and T-stiffeners fET may be calculated as:
	if self.StiffenerType == 'L' or self.StiffenerType == 'T':
	    # (7.33)
	    #_Iz = (((1.0/12.0) * _Af * _b**2) 
	    #       + (_ef**2 * (_Af / (1.0 + _Af / _Aw))))
	    #
	    # (7.32)
	    _fET = ((_Beta * ((_Aw + (self.tf / self.tw)**2 * _Af) 
                              / (_Aw + 3 * _Af)) 
                     * self.G * (self.tw / self.hw)**2) 
                    + ((math.pi**2 * self.E * _Iz) 
                       / ((_Aw / 3.0 + _Af) * _LT**2)))
	    #
	#
	# For flatbar stiffeners fET may be calculated as:
	elif self.StiffenerType == 'F':
	    # (7.34)
	    _fET = ((_Beta + 2 * (self.hw / _LT)**2) 
                    * self.G * (self.tw / self.hw)**2)
	    #
	#
	# Generally fET may be calculated as:
	else:
	    _G = self.E / (2.0 * ( 1.0 + self.Poisson)) 
	    # (7.31)
	    _fET = ((_Beta * (_G * _It) / (Ipo)) 
                    + (math.pi**2 * (self.E * _hs**2 * Iz) 
                   / (_Ipo * _LT**2)))
	    #
	#
	# (7.30)
	_LambdaT = math.sqrt(self.fy / _fET)
	#
	# (7.29)
	_MuT = 0.35 * (_LambdaT - 0.60)
	#
	#
	# The torsional buckling strength may be calculated as:
	#
	# (7.28)
	if _LambdaT > 0.60:
	    #
	    _ft = (self.fy * (1 + _MuT + _LambdaT**2 
                              - math.sqrt((1 + _MuT
                                           + _LambdaT**2)**2 
                                          - 4 * _LambdaT**2))
                   / (2 * _LambdaT**2))
	    #
	    #print('==>',1 + _MuT + _LambdaT**2 )
	    #print(math.sqrt((1 + _MuT
	    #                 + _LambdaT**2)**2 
	    #                - 4 * _LambdaT**2) )
	    #print(2 * _LambdaT**2)
	    #
	#
	# (7.27)
	else:
	    _ft = self.fy
	#
	#
	#
	return _fET, _LambdaT, _MuT, _ft
	#
    #
    #
    def fk(self, _side, _fT, _fE, _ie, _Zt, _Zp):
	#
	# for check at stiffener side
	if _side == 'P':
	    _fr = self.fy
	#
	elif _side == 'S':
	    _fr = _fT
	#
	# for check at stiffener side if LambdaT <= 0.6,
	else:
	    self.msg = self.msg + '\n' + '**error**'
	    exit
	#
	# (7.23)
	_Lambda = math.sqrt(_fr / _fE)
	#
	#
	# for check at plate side
	# (7.25)
	if _side == 'P':
	    _Mu = ((0.34 + 0.08 * (_Zp / _ie)) 
                   * (_Lambda - 0.20))
	#
	# for check at stiffener side
	# (7.26)
	elif _side == 'S':
	    _Mu = ((0.34 + 0.08 * (_Zt / _ie)) 
                   * (_Lambda - 0.20))
	#
	else:
	    self.msg = self.msg + '\n' + '**error**'
	    exit
	#
	#
	#
	# The characteristic buckling strength for stiffeners 
	# may be found from:
	#
	# (7.22)
	# Plate
	if _Lambda > 0.20:
	    _fk = (_fr * (1 + _Mu + _Lambda**2 
                          - math.sqrt((1 + _Mu 
                                       + _Lambda**2)**2 
                                      - 4 * _Lambda**2))
                   / (2 * _Lambda**2))
	#
	# (7.21)
	else:
	    _fk = _fr
	#
	#
	#
	return _fr, _Lambda, _Mu, _fk
	#
    #
    #
    # eq (8.6)
    def TauCR(self, _LambdaCE):
	#
	if _LambdaCE > 1.0:
	    return (0.60*self.fy / _LambdaCE**2)
	#
	else:
	    return (0.60*self.fy)
	#
    # 
    #
    # eq (8.27)
    def fTG(self, _LambdaTG, _Mu):
	# 
	if _LambdaTG > 0.60 :
	    #
	    _fTG = (self.fyG * ((1.0 + _Mu + _LambdaTG**2 
                                 - math.sqrt((1.0 + _Mu + _LambdaTG**2)**2 
                                            - 4.0 * _LambdaTG**2)) 
                                / (2.0 * _LambdaTG**2)))
	    #
	#
	# _LambdaTG <= 0.60
	else:
	    #
	    _fTG = self.fyG
	#
	return _fTG
	#
    #
    #
    def frG(self, ftG, side):
	#
	if side == 'P':
	    return self.fyG
	#
	else:
	    return ftG
	#
    #
    #
    # Combined Unity Check 
    # equations (7.50) to (7.57) or (7.59) to (7.64)
    def CombinedUnityCheck(self, _Stiffener, _NSd, _L, _qSdp, _qSds, _NksRd, _NkpRd, _MpRd, _MstRd, _Ms1Rd, _Ms2Rd, _NE, _NRd,  _z = 0.0 , _u = 0):
	#
	# In the equations (7.50) to (7.57) or (7.59) to (7.62) u = 0 for
	# girders.
	#
	# Girders may be checked for shear forces similar to stiffeners
	# see sec. 7.8.
	#
	# For simplification, assuming:
	self.msg = self.msg + '\n' + 'z =' + str(_z )
	# where z is the distance from the neutral axis of the effective
	# section to the working point of axial force.
	#
	self.msg = self.msg + '\n' + 'u =' + str( _u)
	#
	# M1Sd for continuous stiffeners with equal spans
	#      and equal lateral pressure in all spans
	#      = absolute value of the actual largest support
	#      moment for continuous stiffeners with unequal spans
	#      and/or unequal lateral pressure in adjacent spans
	#
	# Lateral pressure on plate side:
	#print('')
	self.msg = self.msg + '\n' + 'Lateral pressure on plate side'
	self.msg = self.msg + '\n' + '----------'
	#
	# Maximum end support moment
	_M1Sd = abs(_qSdp * _L**2 / 12.0)
	self.msg = self.msg + '\n' + 'M1Sd =' +str(_M1Sd ) + ' ' + str(_L)
	#
	# Maximum mid-beam moment
	_M2Sd = _M1Sd / 2.0
	self.msg = self.msg + '\n' + 'M2Sd =' + str (_M2Sd )
	#
	# (7.50)
	_UR_p750 = ((_NSd / _NksRd) 
                    + ((_M1Sd - _NSd * _z) 
                       / (_Ms1Rd * (1.0 - _NSd / _NE))) 
                    + _u)
	#
	self.msg = self.msg + '\n' + 'URp eq7.50' + str( _UR_p750)
	#
	# (7.51)
	_UR_p751 = ((_NSd / _NkpRd) - (2 * _NSd / _NRd)
                    + ((_M1Sd - _NSd * _z) 
                       / (_MpRd * (1.0 - _NSd / _NE))) 
                    + _u)
	#
	self.msg = self.msg + '\n' + 'URp eq7.51' + str(_UR_p751)
	#
	# (7.52)
	_UR_p752 = ((_NSd / _NksRd) - (2 * _NSd / _NRd)
                    + ((_M2Sd + _NSd * _z) 
                       / (_MstRd * (1.0 - _NSd / _NE))) 
                    + _u)
	#
	self.msg = self.msg + '\n' + 'URp eq7.52' + str(_UR_p752)
	#
	# (7.53)
	_UR_p753 = ((_NSd / _NkpRd) 
                    + ((_M2Sd - _NSd * _z) 
                       / (_MpRd * (1.0 - _NSd / _NE))) 
                    + _u)
	#
	self.msg = self.msg + '\n' + 'URp eq7.53' + str(_UR_p753)
	#
	_URp = max(_UR_p750, _UR_p751, _UR_p752, _UR_p753)
	self.msg = self.msg + '\n' + 'URp =' + str(_URp)
	#
	#
	# Lateral pressure on stiffener/Girder side:
	#print('')
	self.msg = self.msg + '\n' + 'Lateral pressure on stiffener/Girder side:'
	#
	_M1Sd = abs(_qSds * _L**2 / 12.0)
	#
	#print('----------')
	self.msg = self.msg + '\n' + 'M1Sd =' + str( _M1Sd ) + ' ' + str( _L)
	#
	_M2Sd = _M1Sd / 2.0
	#
	self.msg = self.msg + '\n' + 'M2Sd =' + str( _M2Sd)
	#
	# Note: The above moments calculated are for the case of equal
	#       span with equal lateral pressure. If the stiffener's spans
	#       or lateral pressure are unequal, the above moments should
	#       be re-calculated to suit.
	#
	#
	# (7.54)
	_UR_s754 = ((_NSd / _NksRd) - (2 * _NSd / _NRd)
                    + ((_M1Sd + _NSd * _z) 
                       / (_MstRd * (1.0 - _NSd / _NE))) 
                    + _u)
	#
	self.msg = self.msg + '\n' + 'URp eq7.54' + str( _UR_s754)+ ' ' + str(_M1Sd + _NSd * _z) + ' ' + str((_MstRd * (1.0 - _NSd / _NE)))
	#
	# (7.55)
	_UR_s755 = ((_NSd / _NkpRd) 
                    + ((_M1Sd + _NSd * _z) 
                       / (_MpRd * (1.0 - _NSd / _NE))) 
                    + _u)
	#
	self.msg = self.msg + '\n' + 'URp eq7.55' + str( _UR_s755)
	#
	# (7.56)
	_UR_s756 = ((_NSd / _NksRd) 
                    + ((_M2Sd - _NSd * _z) 
                       / (_Ms2Rd * (1.0 - _NSd / _NE))) 
                    + _u)
	#
	self.msg = self.msg + '\n' + 'URp eq7.56' + str( _UR_s756)
	#
	# (7.57)
	_UR_s757 = ((_NSd / _NkpRd) - (2 * _NSd / _NRd)
                    + ((_M2Sd - _NSd * _z) 
                       / (_MpRd * (1.0 - _NSd / _NE))) 
                    + _u)
	#
	self.msg = self.msg + '\n' + 'URp eq7.57' +str(_UR_s757)
	#
	_URs = max(_UR_s754, _UR_s755, _UR_s756, _UR_s757)
	self.msg = self.msg + '\n' + 'URs =' + str( _URs)
	#
	#
	#
	#
	# 7.7.2 Simple supported stiffener (sniped stiffener)
	# ---------------------------------------------------
	#
	#print(' ')
	self.msg = self.msg + '\n' + 'Simple supported stiffener (sniped stiffener)'
	#
	# Simple supported stiffener (sniped stiffener):
	#
	# Lateral pressure on plate side:
	# (7.59)
	_UR_p759 = ((_NSd / _NksRd) - (2 * _NSd / _NRd) 
                    + ((abs(_qSdp * _L**2 / 8.0) + _NSd * _z) 
                       / (_MstRd * (1.0 - _NSd / _NE))) + _u)
	#
	#print('')
	self.msg = self.msg + '\n' + 'Lateral pressure on plate side :'
	#print('----------')
	self.msg = self.msg + '\n' + 'URp eq7.59 =' + str(_UR_p759)
	#
	# (7.60)
	_UR_p760 = ((_NSd / _NkpRd) 
                    + ((abs(_qSdp * _L**2 / 8.0) + _NSd * _z) 
                       / (_MpRd * (1.0 - _NSd / _NE))) + _u)
	#
	self.msg = self.msg + '\n' + 'URp eq7.60 =' + str( _UR_p760)
	#
	_URps = max(_UR_p759, _UR_p760)
	self.msg = self.msg + '\n' + 'URps =' + str( _URps)
	#
	#
	# Lateral pressure on stiffener side:
	#
	# (_qSd * _L**2 / 8.0) <= (_NSd * _Zc)
	#
	# Section at middle stiffener side
	# (7.61)
	_UR_s761 = ((_NSd / _NksRd) 
                    + ((abs(_qSds * _L**2 / 8.0) - _NSd * _z) 
                       / (_Ms2Rd * (1.0 - _NSd / _NE))) + _u)
	#
	# Section at middle plate side
	# (7.62)
	_UR_s762 = ((_NSd / _NkpRd) - (2 * _NSd / _NRd) 
                    + ((abs(_qSds * _L**2 / 8.0) - _NSd * _z) 
                       / (_MpRd * (1.0 - _NSd / _NE))) + _u)
	#
	#
	#print('')
	self.msg = self.msg + '\n' + 'Lateral pressure on stiffener side:'
	self.msg = self.msg + '\n' + '----------'
	self.msg = self.msg + '\n' + 'URs eq7.61 =' + str(_UR_s761)
	self.msg = self.msg + '\n' + 'URs eq7.62 =' + str( _UR_s762)
	#
	#
	# Section at middle stiffener side
	# (7.63)
	_UR_s763 = ((_NSd / _NksRd) - (2 * _NSd / _NRd) 
                    + (( _NSd * _z - abs(_qSds * _L**2 / 8.0)) 
                       / (_MstRd * (1.0 - _NSd / _NE))) + _u)
	#
	# Section at middle plate side
	# (7.64)
	_UR_s764 = ((_NSd / _NkpRd) 
                    + ((_NSd * _z - abs(_qSds * _L**2 / 8.0)) 
                       / (_MpRd * (1.0 - _NSd / _NE))) + _u)
	#
	#
	self.msg = self.msg + '\n' + '----------'
	self.msg = self.msg + '\n' + 'URs eq7.63 =' + str( _UR_s763)
	self.msg = self.msg + '\n' + 'URs eq7.64 =' + str( _UR_s764)
	#
	#
	if (_qSds * _L**2 / 8.0) < (_NSd * _z):
	    _URss1 = _UR_s763
	    _URss2 = _UR_s764
	#
	else:
	    _URss1 = _UR_s761
	    _URss2 = _UR_s762
	#
	self.msg = self.msg + '\n' + '----------'
	self.msg = self.msg + '\n' + 'URss1 =' + str( _URss1)
	self.msg = self.msg + '\n' + 'URss2 =' + str( _URss2)
	#
	_URss = max(_URss1, _URss2)
	#
	self.msg = self.msg + '\n' + 'URss =' + str( _URss)
	#
	if _Stiffener == 'C':
	    _URstiffP = _URp
	    _URstiffS = _URs
	#
	elif _Stiffener == 'S':
	    _URstiffP = _URps
	    _URstiffS = _URss
	#
	else:
	    self.msg = self.msg + '\n' + 'Check Stiffener'
	    exit
	#
	#print('')
	self.msg = self.msg + '\n' + '----------'
	self.msg = self.msg + '\n' + 'URstiffP =' + str( _URstiffP)
	self.msg = self.msg + '\n' + 'URstiffS =' + str( _URstiffS)
	#
	self.ur_combined = max(_URstiffP, _URstiffS)
	#
	#print('')
	self.msg = self.msg + '\n' + '----------'
	self.msg = self.msg + '\n' + 'UR total combined =' + str( self.ur_combined)
	#
	#
	#
    #
    #
    #
    # Table 3-1
    def PlateBucklingCheck(self):
	#
	# Reference table for buckling checks of plates
	#
	self.msg = self.msg + '\n' + 'ok'
	#
    #
    #
    # Section 5
    def LateralLoadedPlates(self):
	#
	#print(' ')
	self.msg = self.msg + '\n' + 'Lateral Loaded Plates'
	# For plates subjected to lateral pressure, either alone or in
	# combination with in-plane stresses, the stresses may be
	# checked by the following formulas:
	#
	# (5.4)
	SigmajSd = (math.sqrt(self.SigmaxSd**2 + self.SigmaySd**2 
                              - self.SigmaxSd * self.SigmaySd
                              + 3 * self.TauSd**2))
	#
	# (5.3)
	Psix = ((1 - (SigmajSd / self.fy)**2) 
                / math.sqrt(1 - (3.0/4.0) * (self.SigmaySd / self.fy)**2 
                            - 3 * (self.TauSd / self.fy)**2))
	#
	# (5.2)
	Psiy = ((1 - (SigmajSd / self.fy)**2) 
                / math.sqrt(1 - (3.0/4.0) * (self.SigmaxSd / self.fy)**2 
                            - 3 * (self.TauSd / self.fy)**2))
	#
	# (5.1 ) - PSd = design lateral pressure
	self.PSR = (4.0 * (self.fy / self.GammaM) * (self.t/ self.S)**2 
                    * (Psiy + (self.S / self.L)**2 * Psix))
	#
	self.ur_laterl = self.PSd/self.PSR
	self.msg = self.msg + '\n' + 'PSd ='  + str(self.PSd)  + 'PSa =' + str(self.PSR)
	self.msg = self.msg + '\n' + 'UR =' + str(self.PSd/self.PSR)
	#
	#
	# This formula for the design of a plate subjected to lateral
	# pressure is based on yield-line theory, and accounts for the
	# reduction of the moment resistance along the yield-line due
	# to applied in-plane stresses. The reduced resistance is
	# calculated based on von Mises equivalent stress. It is
	# emphasised that the formulation is based on a yield pattern
	# assuming yield lines along all four edges, and will give
	# uncertain results for cases where yield-lines can not be
	# developed along all edges. Furthermore, since the formula
	# does not take account of second-order effects, plates
	# subjected to compressive stresses shall also fulfil the
	# requirements of Chapter 6 and 7 whichever is relevant.
	#
    #
    # Section 6
    #def BucklingOfUnstiffenedPlates(self):
    #
    # 6.1 General
    # This section presents recommendations for calculating the
    # buckling resistance of unstiffened plates.
    # For plates that are part of a stiffened panel, the plate are
    # checked as part of the buckling checks according to Chapter 7.
    # Then additional check of the plate according to this section
    # is not required.
    # Buckling checks of unstiffened plates in compression shall
    # be made according to the effective width method. The
    # reduction in plate resistance for in-plane compressive forces
    # is expressed by a reduced (effective) width of the plate which
    # is multiplied by the design yield strength to obtain the design
    # resistance, see Figure 6-1 of DNV.
    #
    # 6.2 Buckling of unstiffened plates under
    #     longitudinally uniform compression
    #
    def BucklingOfUnstiffenedPlatesLongitudinalUniformCompression(self):
	#
	#print(' ')
	self.msg = self.msg + '\n' + '6.2 Buckling of unstiffened plates under'
	self.msg = self.msg + '\n' + '    longitudinally uniform compression'
	#
	# The design buckling resistance of an unstiffened plate under
	# longitudinal compression force may be calculated as:
	#
	# (6.3)
	_Lambdap = 0.525 * (self.S / self.t) * math.sqrt(self.fy / self.E)
	#
	# (6.2)
	if _Lambdap > 0.673 :
	    _Cx = ((_Lambdap - 0.22) / _Lambdap**2)
	    #
	#
	else: _Cx = 1.0
	#
	# (6.1)
	self.SigmaxRd = (_Cx * (self.fy / self.GammaM))
	#
	# in which
	# s = plate width
	# t = plate thickness
	# fcr = critical plate buckling strength
	# The resistance of the plate is satisfactory when:
	#
	self.ur_longcomp = self.SigmaxSd/self.SigmaxRd
	self.msg = self.msg + '\n' + 'SigmaxSd =' + str( self.SigmaxSd)
	self.msg = self.msg + '\n' + 'SigmaxRd =' + str(self.SigmaxRd)
	self.msg = self.msg + '\n' + 'UR =' + str( self.SigmaxSd/self.SigmaxRd)
	#
    #
    #
    # 6.3 Buckling of unstiffened plates with
    #     transverse compression
    #
    def BucklingOfUnstiffenedPlatesTransverseCompression(self):
	#
	#print(' ')
	self.msg = self.msg + '\n' + '6.3 Buckling of unstiffened plates with'
	self.msg = self.msg + '\n' + '    transverse compression'
	#
	# The design buckling resistance of a plate under transverse
	# compression force may be found from:
	#
	# (6.8)
	_Lambdac = 1.10 * (self.S / self.t) * math.sqrt(self.fy / self.E)
	self.msg = self.msg + '\n' + 'Lambdac =' + str( _Lambdac)
	#
	# (6.9)
	_Mu = 0.21 * (_Lambdac - 0.20)
	self.msg = self.msg + '\n' + 'Mu ='  + str(_Mu)
	#
	# (6.7)
	if _Lambdac >= 2.0:
	    _k = (1.0/(2.0 * _Lambdac**2)) + 0.070
	    #
	#
	elif _Lambdac <= 0.20:
	    _k = 1.0
	#
	# 0.20 > _Lambdac < 2.0
	else:
	    _k = ((1.0/(2.0 * _Lambdac**2)) 
                  * (1 + _Mu + _Lambdac**2 
                     - math.sqrt((1 + _Mu + _Lambdac**2 )**2 
                                 - 4 * _Lambdac**2)))
	#
	self.msg = self.msg + '\n' + 'k = '  +str( _k)
	#
	# t = plate thickness
	# l = plate length
	# s = plate width
	# The reduction factor due to lateral load kp
	# (6.11)
	_halpha = max((0.050 * (self.S/self.t) - 0.75) , 0)
	#
	# (6.10)
	if self.PSd <= (2 * (self.t/self.S)**2 * self.fy):
	    _kp = 1.0
	#
	else:
	    _kp = max(1.0 - _halpha * ((self.PSd/self.fy) - 2 * (self.t/self.S)**2),0.0)
	#
	# print ('kp =',_kp)
	#
	# (6.6)
	self.SigmayR = (((1.3 * self.t / self.L) 
                         * math.sqrt(self.E/self.fy) 
                         + _k * (1.0 - (1.3 * self.t / self.L) 
                                 * math.sqrt(self.E/self.fy))) 
                        * self.fy * _kp)
	#
	self.msg = self.msg + '\n' + 'SigmayR = ' + str( self.SigmayR)
	#
	# (6.5)
	self.SigmayRd = self.SigmayR / self.GammaM
	self.msg = self.msg + '\n' + 'SigmayRd= ' + str( self.SigmayRd)
	#
	# The resistance of the plate is satisfactory when:
	# (6.12)
	self.ur_transComp = self.SigmaySd/ self.SigmayRd
	self.msg = self.msg + '\n' + 'UR = '  + str(self.SigmaySd/ self.SigmayRd)
	#
	#
    #
    #
    #
    # 6.4 Buckling of unstiffened plate with shear
    #
    def BucklingOfUnstiffenedPlatesShear(self):
	#
	# print (' ')
	# print ('6.4 Buckling of unstiffened plate with shear')
	#
	# Shear buckling of a plate can be checked by
	#
	# (6.17)
	if self.L < self.S:
	    _kl = 5.34 * (self.S / self.L)**2 + 4.0
	#
	# L >= S
	else:
	    _kl = 5.34 + 4 * (self.S / self.L)**2
	#
	self.msg = self.msg + '\n' + 'kl ='   + str( _kl)
	#
	# (6.16)
	_Lambdaw = 0.795 * (self.S / self.t) * math.sqrt(self.fy / (self.E * _kl))
	#
	self.msg = self.msg + '\n' + 'Lambdaw ='  + str( _Lambdaw)
	#
	# (6.15)
	if _Lambdaw > 1.20:
	    _Ctau = 0.90 / _Lambdaw 
	#
	elif _Lambdaw <= 0.80:
	    _Ctau = 1.0
	#
	# 0.80 < Lambdaw <= 1.20
	else:
	    _Ctau = 1.0 - 0.625 * ( _Lambdaw - 0.80)
	#
	self.msg = self.msg + '\n' + 'Ctau ='  + str( _Ctau)
	#
	# (6.14)
	self.TauRd = (_Ctau/self.GammaM) * (self.fy / math.sqrt(3.0))
	#
	self.msg = self.msg + '\n' + 'Tau =' +  str( self.TauRd)
	#
	# (6.13)
	self.msg = self.msg + '\n' + 'URs = '   + str( self.TauSd / self.TauRd)
	#
    #
    #
    # 6.5 Buckling of unstiffened biaxially loaded
    #     plates with shear
    #
    def BucklingOfUnstiffenedPlatesBiaxiallyLoadedShear(self):
	#
	# print (' ')
	# print ('6.5 Buckling of unstiffened biaxially loaded')
	# print ('    plates with shear')
	#
	#
	# A plate subjected to biaxially loading with shear should fulfil
	# the following requirement:
	#
	# where if both SigmaxSd and SigmaySd is compression (positive) then
	if self.Stress == 'COMPRESSION':
	    #
	    if (self.S / self.t) > 1.20:
		_ci = 0
	    #
	    else: _ci = (1.0 - (self.S / (120 * self.t)))
	    #
	#
	# If either of SigmaxSd and SigmaySd or both is in tension (negative), then
	else: _ci = 1.0
	#
	# SigmaxRd is given by eq. (6.1) and SigmayRd is given by eq. (6.5).
	# In case of tension, apply fy/GammaM.
	# TauRd is given by eq. (6.19) in cases where SigmaySd is positive
	# (compression) and by eq. (6.14) in cases where SigmaySd is zero
	# or negative (in tension).
	#
	# Shear buckling of a plate can be checked by
	#
	# (6.17)
	if self.L < self.S:
	    _kl = 5.34 * (self.S / self.L)**2 + 4.0
	#
	# L >= S
	else:
	    _kl = 5.34 + 4 * (self.S / self.L)**2
	#
	print('kl =', _kl)
	#
	# (6.16)
	_Lambdaw = 0.795 * (self.S / self.t) * math.sqrt(self.fy / (self.E * _kl))
	#
	self.msg = self.msg + '\n' + 'Lambdaw =' + str( _Lambdaw)
	#
	# (6.20)
	if _Lambdaw > 1.25:
	    _Ctaue = 1.0 / _Lambdaw**2
	#
	elif _Lambdaw <= 0.80:
	    _Ctaue = 1.0
	#
	# 0.8 > _Lambdaw < 1.25
	else:
	    _Ctaue = 1.0 - 0.80*( _Lambdaw - 0.80)
	#
	# (6.19)
	self.TauRd = (_Ctaue / self.GammaM) * (self.fy / math.sqrt(3))
	#
	self.msg = self.msg + '\n' + 'Tau =' + str( self.TauRd)
	#
	# (6.18)
	self.URd = ((self.SigmaxSd / self.SigmaxRd)**2 
                    + (self.SigmaySd / self.SigmayRd)**2 
                    - _ci * (self.SigmaxSd / self.SigmaxRd)*(self.SigmaySd / self.SigmayRd)
                    + (self.TauSd / self.TauRd)**2)
	#
	self.ur_biaxial = self.URd
	self.msg = self.msg + '\n' + 'UR = ' + str( self.URd)
	#
    #
    #
    # 6.6 Buckling of unstiffened plates with varying
    #     longitudinal stress. Internal compression
    #     elements
    #
    def BucklingOfUnstiffenedPlatesVaryingLongStress(self):
	#
	# print (' ')
	# print ('6.6 Buckling of unstiffened plates with varying')
	# print ('    longitudinal stress. Internal compression')
	# print ('    elements')
	#
	#
	# The buckling resistance of an unstiffened plate with
	# varying longitudinal stress may be found from:
	#
	# s = plate width
	# Psi = Sigma2/Sigma1 Stress ratio. s1 is largest stress with
	# compressive stress taken as positive.
	# t = plate thickness
	# fcr = critical plate buckling strength
	#
	#_Psi = self.Sigma2 / self.Sigma1
	_Psi = self.Sigmax2Sd / self.Sigmax1Sd
	#
	_Epsilon = math.sqrt(235.0 / self.fy)
	#
	if _Psi <= 0:
	    #
	    if _Psi >= -2.0:
		_ksigma = 2.98 * (1.0 - _Psi)**2
	    #
	    # _Psi >= -1.0:
	    else:
		_ksigma = 7.81 - 6.29 * _Psi + 9.78 * _Psi**2
	#
	# _Psi > 0:
	else:
	    if _Psi <= 1.0:
		_ksigma = 8.20 / (1.05 - _Psi)
	    #
	    else:
		# print ('Psi > 1.0 ==> ksigma not available')
		pass
	    #
	#
	# where Lambdap is the plate slenderness given by:
	# (6.24)
	_Lambdap = (self.S / self.t) * 1.0 / (28.4 * _Epsilon * math.sqrt(_ksigma))
	# 
	# (6.23)
	if _Lambdap > 0.73 :
	    _Cx = (_Lambdap - 0.055 * (3 + _Psi)) / _Lambdap**2
	#
	# (6.22)
	else:
	    _Cx = 1.0
	#
	# (6.21)
	self.SigmaxRd = _Cx * self.fy / self.GammaM
	self.msg = self.msg + '\n' + 'SigmaxRd ='  + str( self.SigmaxRd)
	#
	# The resistance of the plate is satisfactory when:
	self.URd = self.SigmaxSd/self.SigmaxRd
	# print('URd =', self.URd)
	#
	# In order to perform cross sectional checks for members
	# subjected to plate buckling the local buckling effects 
	# can be accounted for by checking the resistance by using
	# the effective width according to Table 6-1.
	#
	#
    #
    #
    # 6.7 Buckling of outstand compression elements
    #
    def BucklingOfUnstiffenedPlatesOutstandCompElemnts(self):
	#
	self.msg = self.msg + '\n' + ' '
	self.msg = self.msg + '\n' + '6.7 Buckling of outstand compression elements'
	#
	# The buckling resistance of an outstand compression element
	# with varying or constant longitudinal stress may be found
	# from:
	#
	# in which
	# s = plate width
	# t = plate thickness
	# fcr = critical plate buckling strength
	#
	_Epsilon = math.sqrt(235.0/ self.fy)
	#
	self.msg = self.msg + '\n' + 'Epsilon =' + str( _Epsilon)
	#
	# For outstand with largest compression stress at free edge:
	if self.OutstandEdge == 'FREE':
	    #
	    if _Psi >= -1.0 and _Psi <= 1.0:
		_ksigma = 0.57 - 0.21 * _Psi + 0.07 * _Psi**2
	    #
	    else:
		#
		self.msg = self.msg + '\n' + 'ERROR ==> Psi < -1 or Psi > 1'
		pass
		#
	#
	# For outstand with largest compression stress at supported
	# edge:
	else:
	    #
	    if _Psi >= 0 and _Psi <= 1.0:
		#
		_ksigma = 0.578 / (0.34 + _Psi)
	    #
	    #
	    elif _Psi >= -1 and _Psi <= 0:
		#
		_ksigma = 1.70 - 5.0 * _Psi + 17.10 * _Psi**2
	    #
	    else:
		#
		self.msg = self.msg + '\n' + 'ERROR ==> Psi < -1 or Psi > 1'
		pass
	    #
	#
	self.msg = self.msg + '\n' + 'ksigma ='  + str( _ksigma)
	#
	# where Lambdap is the plate slenderness given by:
	# (6.29)
	_Lambdap = ((self.S / self.t) 
                    * (1.0 / (28.4 * _Epsilon * math.sqrt(_ksigma))))
	#
	# (6.28)
	if _Lambdap > 0.749:
	    #
	    _Cx = (_Lambdap - 0.188) / _Lambdap**2
	#
	# when _Lambdap <= 0.749
	else:
	    _Cx = 1.0
	#
	self.msg = self.msg + '\n' + 'Cx =' + str(_Cx)
	#
	# (6.26)
	self.SigmaxRd = _Cx * self.fy / self.GammaM
	#
	self.msg = self.msg + '\n' + 'SigmaxRd ='  + str(self.SigmaxRd)
	#
	# Cross sectional checks of members subjected to plate
	# buckling local buckling effects can be accounted for by
	# checking the resistance by using the effective width
	# according to Table 6-2 and Table 6-3 for outstand elements
	# with largest compression stress at free edge or supported
	# edge respectively.
	#
    #
    #
    # 6.8 Buckling of unstiffened plates with varying
    # transverse stress
    def BucklingOfUnstiffenedPlates(self):
	if self.LoadCase == 7:
	    #
	    # In case of linear varying transverse stress the capacity check
	    # can be done by use of the design stress value at a distance l1
	    # from the most stressed end of the plate, but not less than 0.75
	    # of maximum sy,Sd. The resistance sy,Rd should be calculated
	    # from eq. (6.5).
	    # l1 = minimum of 0.25 l and 0.5 s
	    self.msg = self.msg + '\n' + '6.8 not finish'
	    pass
	    #
	#
	#
	# 6.9 Buckling of unstiffened plate with
	#     longitudianal and transverse varying
	#      stress and with shear stress
	elif self.LoadCase == 8:
	    #
	    # The check of combined varying loads may be done according
	    # to eq. (6.18) with the resistance calculated according to eq.
	    # (6.21) and eq. (6.5) using the stress point defined in sec. 6.8.
	    self.msg = self.msg + '\n' + '6.9 not finish'
	    pass
	    #
	#
	#
	else:
	    self.msg = self.msg + '\n' + 'Case not defined'
	    pass
	#
	#
    #
    #
    # Section 7
    def BucklingOfStiffenedPlates(self):
	#
	self.msg = self.msg + '\n' + ''
	self.msg = self.msg + '\n' + 'Buckling of Stiffened Plates'
	#
	# 7.1 General
	# -----------
	#
	# This chapter deals with stiffened plate panels subjected to
	# axial stress in two directions, shear stress and lateral load.
	# There are different formulas for stiffeners being continuous
	# (or connected to frames with their full moment resistance)
	# and simple supported (sniped) stiffeners.
	#
	# An example of a stiffened plate panel is shown in Figure 3-1.
	# The stiffener cross section needs to fulfil requirements to
	# avoid local buckling given in Chapter 9.
	#
	# For shear lag effects see Commentary Chapter 10.
	#
	# The plate between stiffeners will normally be checked
	# implicitly by the stiffener check since plate buckling is
	# accounted for by the effective width method. However, in
	# cases where sy,Sd stress is the dominant stress it is necessary
	# to check the plate resistance according to eq. (7.19).
	#
	# For slender stiffened plates the load carrying resistance in the
	# direction transverse to the stiffener may be neglected. Then
	# sy,Sd stresses may be assumed to be carried solely by the
	# girder. In such cases the effective girder flange may be
	# determined by disregarding the stiffeners, and the stiffener
	# with plate may be checked by neglecting sy,Sd stresses
	# (method 2 in sec. 8.4). See also Commentary to 8 in Chapter
	# 10.
	#
	#
	self.msg = self.msg + '\n' + ''
	self.msg = self.msg + '\n' + '7.2 Forces in the idealised stiffened plate'
	self.msg = self.msg + '\n' + ''
	#
	# 7.2 Forces in the idealised stiffened plate
	# -------------------------------------------
	#
	# Stiffened plates subjected to combined forces, see Figure 7-1
	# should be designed to resist an equivalent axial force
	# according to eq. (7.1) and an equivalent lateral load
	# according to eq. (7.8).
	#
	# Assumption of tension field action implies that no (or
	# negligible) resistance of the plate against transverse
	# compression stresses (Sigmay) can be assumed. See also
	# Commentary Chapter 10.
	#
	#
	# Where
	# (7.7)
	if self.L < self.S:
	    #
	    _kl = 5.34 * (self.S / self.L)**2 + 4
	#
	# self.L >= self.S
	else:
	    #
	    _kl = 5.34 + 4*(self.S / self.L)**2
	#
	self.msg = self.msg + '\n' + 'kl ='  + str(_kl)
	#
	# Taucrl = critical shear stress for the plate panel 
	# between twostiffeners, according to eq. (7.6).
	# (7.4)
	self._Taucrl = _kl * 0.904 * self.E * (self.t / self.S)**2
	print ('Tau = ' + str(self._Taucrl))
	#
	self.msg = self.msg + '\n' + 'Taucrl : ' + str( self._Taucrl)
	#
	# (7.5)
	if self.L > self.LG:
	    #
	    _kg = 5.34 *(self.L / self.LG)**2 + 4
	#
	# self.L <= self.LG
	else:
	    #
	    _kg = 5.34 + 4*(self.L / self.LG)**2
	#
	self.msg = self.msg + '\n' + 'kg ='  + str( _kg)
	#
	# Taucrg = critical shear stress for the plate with the 
	# stiffeners removed, according to eq. (7.4).
	# (7.6)
	_Taucrg = _kg * 0.904 * self.E * (self.t / self.L)**2
	#
	self.msg = self.msg + '\n' + 'Taucrg : '  + str( _Taucrg)
	#
	# (7.2)
	if self.TauSd > (self._Taucrl / self.GammaM):
	    # tension field action is allowed
	    self.Tautf = self.TauSd - _Taucrg
	#
	# (7.3)
	else:
	    self.Tautf = 0
	#
	self.msg = self.msg + '\n' + 'Tautf : ' + str(self.Tautf)
	#
	# The equivalent axial force should be taken as:
	#
	# where
	# As = cross sectional area of stiffener
	# s = distance between stiffeners
	# t = plate thickness
	# SigmaxSd = axial stress in plate and stiffener
	#            with compressive stresses as positive
	#
	# (7.1)
	self.NSd = (self.SigmaxSd * (self.As + self.S * self.t) 
                    + self.Tautf * self.S * self.t)
	#
	self.msg = self.msg + '\n' + 'Equivalent axial force NSd : ' + str(self.NSd)
	#
	# Sigmay1Sd = larger design stress in the transverse direction,
	# with tensile stresses taken as negative
	#
	# Sigmay2Sd = smaller design stress in the transverse direction,
	# with tensile stresses taken as negative
	#
	_Psi = self.Sigmay2Sd / float(self.Sigmay1Sd)
	#
	self.msg = self.msg + '\n' + 'Psi =' + str(_Psi)
	#
	# Wes = section modulus for stiffener with effective plate
	# at flange tip
	# mc = 13.3 for continuous stiffeners or,
	# mc = 8.9 for simple supported stiffeners (sniped stiffeners)
	# Is = moment of inertia of stiffener with full plate width
	#
	# (7.12)
	_kc = (2 * (1 + math.sqrt(1 + (10.9 * self.Is) 
                                  / (self.t**3 * self.S))))
	#
	self.msg = self.msg + '\n' + 'kc : ' + str(_kc)
	#
	#
	#
	# 7.3 Effective plate width
	# -------------------------
	#
	self.msg = self.msg + '\n' + ''
	self.msg = self.msg + '\n' + '7.3 Effective plate width'
	#
	# Only applicable for continous stifener??
	#
	# (7.15)
	_Lambdap = 0.525 * (self.S / self.t) * math.sqrt(self.fy / self.E)
	#
	self.msg = self.msg + '\n' + 'Lambdap : ' + str( _Lambdap)
	#
	# The reduction factor due to stresses in the longitudinal
	# direction, Cxs, is
	#
	# (7.14)
	if _Lambdap > 0.673:
	    self.Cxs = (_Lambdap - 0.22)/ _Lambdap**2
	#
	else:
	    self.Cxs = 1.0
	#
	self.msg = self.msg + '\n' + 'Cxs : '  +str(self.Cxs)
	# 
	if (self.S / self.t) > 120:
	    _ci = 0
	#
	else:
	    _ci = 1 - self.S / (120 * self.t)
	#
	self.msg = self.msg + '\n' + 'ci : ' + str( _ci)
	#
	#
	# and the reduction factor for compression stresses in the
	# transverse direction, Cys, is found from:
	if self.Stress == 'COMPRESSION':
	    #
	    # (7.16)
	    self.Cys = math.sqrt(1 - (self.SigmaySd / self.SigmayR)**2 
                                 + _ci * ((self.SigmaxSd * self.SigmaySd)
                                          / (self.Cxs * self.fy * self.SigmayR)))
	    #
	    # SigmayR is calculated according to eq. (6.6).
	    # In case of linear varying stress, SigmaySd may be determined as
	    # described in sec. 6.8
	#
	# The reduction factor for tension stresses in the transverse
	# direction, Cys, is calculated as:
	else:
	    #
	    # (7.17)
	    self.Cys1 = ((1.0/2.0) * (math.sqrt(4 - 3 * (self.SigmaySd / self.fy)**2) 
                                      + (self.SigmaySd / self.fy)))
	    #
	    self.Cys = min(self.Cys1, 1.0)
	    #
	    # Tensile stresses are defined as negative.
	    #
	#
	self.msg = self.msg + '\n' + 'Cys =' + str(self.Cys)
	#
	#
	# The effective plate width for a continuous stiffener subjected
	# to longitudinal and transverse stress and shear is calculated
	# as:
	# (7.13)
	self.Se = self.S * self.Cxs * self.Cys
	#
	self.msg = self.msg + '\n' + 'Se =' + str(self.Se)
	#
	#
	#
	# The following resistance parameters are used in the
	# interaction equations for stiffeners:
	#
	# effective elastic section modulus on plate side,
	# see Figure 7-3.
	#_Wep = _Ie / _Zp
	#
	# effective elastic section modulus on stiffener 
	# side, see Figure 7-3.
	#_Wes = _Ie / _Zt
	#
	# effective area of stiffener and plate
	#_Ae = (_As + self.Se * self.t) 
	#
	# As = cross sectional area of stiffener
	# Se = effective width, see sec. 7.3
	#
	# Section Property of Stiffener with Efective Plate
	self.msg = self.msg + '\n' + ''
	self.msg = self.msg + '\n' + 'Section Property of Stiffener with Efective Plate'
	#
	_Ae, _Aec, _Zp, _Ie, _Iz, _Wes, _Wep, _fy = PlateStiffener(self.Se, self.t,
                                                                   self.hw, self.tw,
                                                                    self.bf, self.tf, 
                                                                    self.L, 
                                                                    self.fyp, self.fyS, self.fyS,
                                                                    self.ES)
	#
	#
	_Zt = (self.t + self.hw + self.tf - _Zp)
	self.msg = self.msg + '\n' + 'fy ='  + str(_fy)
	#self.msg = self.msg + '\n' + 'Flang Cross Sectional Area Af =' + str( _Af)
	#self.msg = self.msg + '\n' + 'Web Cross Sectional Area Aw =' + str( _Aw)
	self.msg = self.msg + '\n' + 'Total Cross Sectional Area Ae =' + str( _Ae)
	self.msg = self.msg + '\n' + 'Neautral Axis Loacation Zt =' + str( _Zt)
	self.msg = self.msg + '\n' + 'Neautral Axis Loacation Zp =' + str( _Zp)
	self.msg = self.msg + '\n' + 'Moment of Inertia Ie =' + str( _Ie)
	self.msg = self.msg + '\n' + 'Torsional Moment of Inertia Iz =' + str( _Iz)
	self.msg = self.msg + '\n' + 'Section Modulus = Wes' + str( _Wes)
	self.msg = self.msg + '\n' + 'Section Modulus = Wep' + str(_Wep)
	self.msg = self.msg + '\n' + ''
	#
	#
	#
	# (7.11)
	_C0 = ((_Wes * self.fy * self.mc) 
               / (_kc * self.E * self.t**2 * self.S))
	#
	self.msg = self.msg + '\n' + 'C0 ='  + str( _C0)
	#
	# p0 shall be applied in the direction of the external pressure
	# pSd. For situations where pSd is less than p0, the stiffener need
	# to be checked for p0 applied in both directions (i.e. at plate
	# side and stiffener side).
	#
	if self.SigmaySd > 0:
	    # (7.9)
	    if _Psi > -1.50:
		#
		_p0 = (0.60 + 0.40 * _Psi) * _C0 * self.Sigmay1Sd
	    #
	    # (7.10)
	    else: _p0 = 0

	else: _p0 = 0
	#
	self.msg = self.msg + '\n' + 'p0 =' + str( _p0)
	#
	# The equivalent lateral line load should be taken as:
	#
	# pSd = design lateral pressure
	# s = stiffener spacing
	#
	# (7.8)
	#
	# Equivalent lateral line load on plate side
	if self.PSd < 0:
	    #
	    if _p0 < abs(self.PSd):
		self.qsdp = 0
	    #
	    else:
		self.qsdp = (self.PSd + _p0) * self.S
	#
	else:
	    self.qsdp = (self.PSd + _p0) * self.S
	#
	self.msg = self.msg + '\n' + 'qsdp =' + str( self.qsdp)
	#
	#
	# Equivalent lateral line load on stiffener side
	if self.PSd < 0:
	    self.qsds = (_p0 - self.PSd) * self.S
	#
	else:
	    #
	    if _p0 < abs(self.PSd):
		self.qsds = 0
	    #
	    else:
		self.qsds = (_p0 - self.PSd) * self.S
	    #
	#
	self.msg = self.msg + '\n' + 'qsds =' + str( self.qsds)
	#
	#
	#
	# 7.4 Resistance of plate between stiffeners
	# ------------------------------------------
	#
	self.msg = self.msg + '\n' + ''
	self.msg = self.msg + '\n' + '7.4 Resistance of plate between stiffeners'
	#
	# The plate between stiffeners shall be checked for:
	#
	# (7.20)
	_ksp = math.sqrt(1.0 - 3.0 * (self.TauSd / self.fy)**2)
	#
	self.msg = self.msg + '\n' + 'ksp =' + str( _ksp)
	#
	# (7.19)
	_URs = (_ksp * self.SigmayRd)/ self.SigmaySd
	# SigmayRd is determined from eq. (6.5).
	#
	self.msg = self.msg + '\n' + 'URs ='  + str( _URs)
	#
	# (7.18)
	self.TauRd = self.fy / (math.sqrt(3.0) * self.GammaM)
	#
	self.msg = self.msg + '\n' + 'TauRd =' + str(self.TauRd )
	#
	_URt = self.TauSd / self.TauRd 
	#
	self.msg = self.msg + '\n' + 'URt =' + str( _URt)
	#
	_URs = self.SigmaySd / (_ksp * self.SigmayRd)
	#
	self.msg = self.msg + '\n' + 'URs ='  +  str( _URs)  + 'SigmaySd ='  + str( self.SigmaySd) + 'SigmayRd =' + str(self.SigmayRd)
	#
	#
	#
	# When this check and stiffener check according to sec. 7.7 is
	# carried out it is not necessary to check the plate between
	# stiffeners according to Chapter 6.
	# See also Commentary Chapter 10.
	#
	#
	# 7.5 Characteristic buckling strength of stiffeners
	# --------------------------------------------------
	#
	self.msg = self.msg + '\n' + ''
	self.msg = self.msg + '\n' + '7.5 Characteristic buckling strength of stiffeners'
	#
	# 7.5.1 General
	# -------------
	#
	# where
	#
	# (7.38)
	self.SigmajSd = (math.sqrt(self.SigmaxSd**2 + self.SigmaySd**2 
                                   - self.SigmaxSd * self.SigmaySd 
                                   + 3 * self.TauSd**2))
	#
	self.msg = self.msg + '\n' + 'SigmajSd ='  + str( self.SigmajSd)
	# where:
	#
	# (7.41)
	_c = 2.0 - (self.S / self.L)
	#
	# (7.42)
	_fEpx = 3.62 * self.E * (self.t / self.S)**2
	#
	# (7.43)
	_fEpy = 0.90 * self.E * (self.t / self.S)**2
	#
	# (7.44)
	_fEptau = 5.0 * self.E * (self.t / self.S)**2
	#
	self.msg = self.msg + '\n' + 'fEptau =' + str( _fEptau)
	# 
	# SigmaxSd and SigmaySd should be set to zero if in tension
	if self.Stress == 'COMPRESSION':
	    #
	    _SigmaxSd = self.SigmaxSd
	    _SigmaySd = self.SigmaySd
	#
	else:
	    #
	    _SigmaxSd = 0
	    _SigmaySd = 0
	#
	#
	# (7.40)
	_Lambdae = (math.sqrt((self.fy / self.SigmajSd) 
                              * math.pow(((_SigmaxSd / _fEpx)**_c 
                                          + (_SigmaySd / _fEpy)**_c 
                                          + (self.TauSd / _fEptau)**_c)
                                         , (1.0 / _c))))
	#
	self.msg = self.msg + '\n' + 'Lambdae =' + str(_Lambdae)
	#
	# (7.39)
	_fep = (self.fy / math.sqrt(1 + math.pow(_Lambdae, 4.0)))
	#
	self.msg = self.msg + '\n' + 'fep =' + str( _fep)
	#
	# (7.37)
	_Eta = min((self.SigmajSd / _fep), 1.0)
	#
	self.msg = self.msg + '\n' + 'Eta =' + str( _Eta)
	#
	# t = plate thickness
	# tf = thickness of flange
	# tW = thickness of web
	# (7.36)
	_C = ((self.hw / self.S) * (self.t / self.tw) 
              * math.sqrt( 1 - _Eta))
	#
	self.msg = self.msg + '\n' + 'C =' + str( _C)
	#
	# Beta = 1.0
	# 
	if self.S > self.L:
	    _Beta = 1.0
	# or may for stocky plates alternatively be calculated as
	# per eq. (7.35) for s <= l
	# (7.35)
	else:
	    _Beta = ((3 * _C + 0.20) / (_C + 0.20))
	#
	self.msg = self.msg + '\n' + 'Beta =' + str(_Beta)
	#
	#
	#
	#
	# 7.5.2 Torsional buckling of stiffeners
	# --------------------------------------
	#
	self.msg = self.msg + '\n' + ''
	self.msg = self.msg + '\n' + '7.5.2 Torsional buckling of stiffeners'
	#
	# The torsional buckling strength may be calculated as:
	#
	# Were :
	#
	# Beta = 1.0,
	#        or may for stocky plates alternatively be 
	#        calculated as per eq. (7.35) for s = l
	# Af = cross sectional area of flange
	# AW = cross sectional area of web
	# G = shear modulus
	# Ipo = polar moment of inertia = Integral (r**2 dA)
	# where r is measured from the connection between the
	#         stiffener and the plate
	# It = stiffener torsional moment of inertia (St. Venant
	#      torsion)
	# Iz = moment of inertia of the stiffeners neutral axis normal
	#      to the plane of the plate
	# b = flange width
	# ef = flange eccentricity, see Figure 7-3
	# hw = web height
	# hs = distance from stiffener toe (connection between
	#      stiffener and plate) to the shear centre of the stiffener
	# lT = distance between sideways supports of stiffener,
	#      distance between tripping brackets (torsional buckling
	#      length).
	# t = plate thickness
	# tf = thickness of flange
	# tW = thickness of web
	#
	# where
	#
	#
	# Generally fET may be calculated as:
	# (7.31)
	#_fET = (_Beta * (self.G * _It / _Ipo) 
	#        + math.pi**2 * ((self.E * _hs**2 * _Iz) 
	#                        / (_Ipo * _LT**2)))
	#
	_fET = {}
	_LambdaT = {}
	_MuT = {}
	_fT = {}
	#
	_fET[0], _LambdaT[0], _MuT[0], _fT[0]  = self.fT(self.Lt, _Beta, _Iz)
	#
	self.msg = self.msg + '\n' + 'fET(x) =' + str(_fET[0]) + ' LambdaT=' + str( _LambdaT[0]) + ' MuT=' +  str( _MuT[0]) + ' fT=' + str(_fT[0]) + ' x = Lt =' + str(self.Lt)
	#
	#
	_LT1 = min(0.40*self.L, self.Lt)
	_fET[1], _LambdaT[1], _MuT[1], _fT[1] = self.fT(_LT1, _Beta, _Iz)
	#
	self.msg = self.msg + '\n' + 'fET(x) =' + str(_fET[1]) + ' LambdaT=' + str(_LambdaT[1]) + ' MuT=' + str( _MuT[1])  + ' fT1=' + str( _fT[1])  + ' x = Lt1 ='+ str(_LT1)
	#
	#
	_LT2 = min(0.80*self.L, self.Lt)
	_fET[2], _LambdaT[2], _MuT[2], _fT[2] = self.fT(_LT2, _Beta, _Iz)
	#
	self.msg = self.msg + '\n' + 'fET(x) =' + str( _fET[2]) + ' LambdaT=' + str( _LambdaT[2]) + ' MuT=' + str( _MuT[2]) + ' fT2=' + str( _fT[2]) + ' x = Lt2 =' + str(_LT2)
	#
	#
	#
	#
	#
	self.msg = self.msg + '\n' + ''
	self.msg = self.msg + '\n' + '7.5.1 General'
	#
	# lk see eq. (7.74)
	#
	# (7.73)
	# Effective radius of gyration
	# Ie effective moment of inertia
	# Ae effective area
	_ie = math.sqrt(_Ie / _Ae)
	#
	self.msg = self.msg + '\n' + 'ie =' +  str(_ie)
	#
	# For a continuous stiffener the buckling length may be
	# calculated from the following equation where :
	#
	# pSd is design lateral pressure
	# Pf is the lateral pressure giving yield in outer-fibre 
	# at support.
	#
	# W = the smaller of Wep and Wes
	# l = span length
	#
	_W = min(_Wep, _Wes)
	#
	self.msg = self.msg + '\n' + 'W =' + str( _W)
	#
	# (7.75)
	_Pf = ((12.0 * _W / (self.L**2 * self.S)) 
               * (self.fy / self.GammaM))
	#
	self.msg = self.msg + '\n' + 'Pf =' + str( _Pf)
	#
	# In case of varying lateral pressure, pSd in eq. (7.74) should be
	# taken as the minimum of the value in the adjoining spans.
	#
	#
	# Lk see eq(7.74)
	if self.Stiffener == 'C':
	    #
	    _Lk = self.L * (1.0 - 0.50 * abs(self.PSd / _Pf))
	#
	elif self.Stiffener == 'S':
	    # For simple supported stiffener 
	    _Lk = 1.0 * self.L
	#
	else:
	    #
	    self.msg = self.msg + '\n' + '****error***'
	    exit
	#
	self.msg = self.msg + '\n' + 'Lk =' +  str(_Lk)
	#
	#
	# where:
	# for check at plate side :
	_frp = self.fy 
	#
	self.msg = self.msg + '\n' + 'frp =' + str( _frp)
	#
	# fT may be calculated according to sec. 7.5.2
	# LambdaT see eq. (7.30)
	#
	#
	# (7.24)
	_fE = math.pi**2 * self.E * (_ie / _Lk)**2
	#
	self.msg = self.msg + '\n' + 'fE =' + str( _fE)
	#
	#
	_side = 'P'
	_fr, _Lambda, _Mu, _fkp = self.fk(_side, _fT[0], _fE, _ie, _Zt, _Zp)
	#
	#
	self.msg = self.msg + '\n' + ''
	self.msg = self.msg + '\n' + '==> ' + 'Side=' + str(_side) + ' fr=' + str( _fr) + ' Lamda=' + str( _Lambda) + ' Mu=' + str( _Mu)
	self.msg = self.msg + '\n' + '_fkp' + str(  _fkp)
	#
	#
	#
	_side = 'S'
	_fr, _Lambda, _Mu, _fks = self.fk(_side, _fT[0], _fE, _ie, _Zt, _Zp)
	#
	#
	self.msg = self.msg + '\n' + ''
	self.msg = self.msg + '\n' + '==> ' + 'Side=' + str( _side) + ' fr=' + str( _fr) + ' Lamda=' + str(_Lambda) + ' Mu=' + str( _Mu)
	# print('_fk', _fks)
	#
	#
	#
	# 7.7.3 Resistance parameters for stiffeners.  
	# ------------------------------------------
	#
	self.msg = self.msg + '\n' + ' '
	self.msg = self.msg + '\n' + '7.7.3 Resistance parameters for stiffeners'
	#
	#
	# (7.65)
	_NRd = _Ae * self.fy / self.GammaM
	#
	self.msg = self.msg + '\n' + 'NRd =' + str( _NRd )
	#
	# where fks is calculated from sec 7.5 using eq(7.26)
	# (7.66)
	_NksRd = _Ae * _fks / self.GammaM
	#
	self.msg = self.msg + '\n' + '-----'
	self.msg = self.msg + '\n' + 'NksRd =' + str( _NksRd )
	# 
	# where fkp is calculated from sec 7.5 using eq(7.25)
	# (7.67)
	_NkpRd = _Ae * _fkp / self.GammaM
	#
	self.msg = self.msg + '\n' + 'NkpRd =' + str( _NkpRd )
	#
	# where fr is calculated from sec. 7.5 for stiffener side using
	#_LT = 0.4 * self.L 
	# or distance between lateral support if this is less.
	# (7.68)
	#
	#
	_side = 'S'
	_fr1, _Lambda, _Mu, _fks1 = self.fk( _side, _fT[1], _fE, _ie, _Zt, _Zp)
	#
	self.msg = self.msg + '\n' + '-----'
	self.msg = self.msg + '\n' + '_fk' + str( _fks1)
	#
	#
	_Ms1Rd = _Wes * _fr1 / self.GammaM
	#
	self.msg = self.msg + '\n' + 'Ms1Rd =' + str( _Ms1Rd)
	#
	# where fr is calculated from sec. 7.5 for stiffener side using
	#_LT = 0.8 * self.L 
	# or distance between lateral support if this is less.
	# (7.69)
	#
	#
	_side = 'S'
	_fr2, _Lambda, _Mu, _fks2 = self.fk( _side, _fT[2], _fE, _ie, _Zt, _Zp)
	#
	self.msg = self.msg + '\n' + '_fk' + str( _fks2)
	#
	_Ms2Rd = _Wes * _fr2 / self.GammaM
	#
	self.msg = self.msg + '\n' + 'Ms1Rd =' + str( _Ms2Rd)
	#
	# (7.70)
	_MstRd = _Wes * self.fy / self.GammaM
	#
	self.msg = self.msg + '\n' + '-----'
	self.msg = self.msg + '\n' + 'MstRd =' + str(_MstRd)
	#
	# (7.71)
	_MpRd = _Wep * self.fy / self.GammaM
	#
	self.msg = self.msg + '\n' + 'MpRd =' + str( _MpRd)
	#
	#
	# (7.72)
	_NE = (math.pi**2 * self.E * _Ae) / (_Lk / _ie)**2
	#
	self.msg = self.msg + '\n' + '-----'
	self.msg = self.msg + '\n' + 'NE =' + str( _NE)
	#
	#
	#
	# 7.6 Resistance of stiffened panels to shear stresses
	# ----------------------------------------------------
	#
	self.msg = self.msg + '\n' + ''
	self.msg = self.msg + '\n' + '7.6 Resistance of stiffened panels to shear stresses'
	#
	# The resistance towards shear stresses TauRd is found as the
	# minimum of TauRdy, TauRdl and TauRds according to the following:
	#
	# where Taucrl is obtained from eq. (7.6) and Taucrs is obtained
	# from:
	#
	# Is= moment of inertia of stiffener with full plate width.
	# (7.49)
	_Ip = (self.t**3 * self.S) / 10.90
	#
	self.msg = self.msg + '\n' + 'Ip =' + str( _Ip)
	#
	# (7.48)
	_Taucrs = ((36.0 * self.E / (self.S * self.t * self.L**2)) 
                   * math.pow((_Ip * self.Is**3), 1.0 / 4.0))
	#
	self.msg = self.msg + '\n' + 'Taucrs =' + str( _Taucrs)
	#
	# (7.47)
	_TauRds = _Taucrs / self.GammaM
	#
	self.msg = self.msg + '\n' + 'TauRds =' + str( _TauRds)
	#
	# (7.45)
	_TauRdy = self.fy / (math.sqrt(3) * self.GammaM)
	#
	self.msg = self.msg + '\n' + 'TauRdy =' + str( _TauRdy)
	#
	# (7.46)
	_TauRdl = self._Taucrl / self.GammaM
	#
	self.msg = self.msg + '\n' + 'TauRdl =' + str( _TauRdl)
	#
	# The resistance towards shear stress in then:
	self.TauRd = min(_TauRdy, _TauRdl, _TauRds)
	#
	self.msg = self.msg + '\n' + 'TauRd =' + str( self.TauRd)
	#
	#
	#
	# 7.7 Interaction formulas for axial compression
	#      and lateral pressure
	# ----------------------------------------------
	#
	# 7.7.1 Continuous stiffeners
	# ---------------------------
	#
	self.msg = self.msg + '\n' + ''
	self.msg = self.msg + '\n' + '7.7.1 Continuous stiffeners'
	#
	# For continuous stiffeners the following four interaction
	# equations need to be fulfilled in case of:
	#
	# When tension field action is assumed according to eq. (7.2)
	# then u = 0.
	# For resistance parameters see sec. 7.7.3 for stiffener and sec.
	# 8.3 for girders.
	#
	# (7.58)
	if self.Tautf == 0:
	    #
	    _u = (self.TauSd / self.TauRd)**2
	#
	else: _u = 0
	#
	self.msg = self.msg + '\n' + 'u ='+ str( _u)
	#
	#
	# qsd is given in eq. (7.8)
	# l  is the span length
	# z* is the distance from the neutral axis of the effective section
	#    to the working point of the axial force. z* may be varied in
	#    order to optimise the resistance. z* should then be selected so
	#    the maximum utilisation found from the equations (7.50) to
	#    (7.53) or (7.54) to (7.57) is at its minimum, see also
	#    Commentary Chapter 10. The value of z* is taken positive
	#    towards the plate. The simplification z* = 0 is always
	#    allowed.
	#
	#
	#
	self.msg = self.msg + '\n' + ''
	self.msg = self.msg + '\n' + '=====> ******'
	#
	self.CombinedUnityCheck(self.Stiffener, self.NSd, self.L, self.qsdp, 
                                self.qsds, _NksRd, _NkpRd, _MpRd, _MstRd, 
                                _Ms1Rd, _Ms2Rd, _NE, _NRd,  self.Z , _u)
	self.ur_combined_Plate = self.ur_combined
	#
	#
	# 7.8 Check for shear force
	# -------------------------
	#
	self.msg = self.msg + '\n' + ' '
	# print('7.8 Check for shear force')
	#
	# The stiffener should in all sections satisfy:
	#
	# where:
	# VSd = design shear force
	# VRd = design shear resistance
	# Anet = net shear area (shear area minus cut outs)
	#
	# If VSd > 0.5 VRd then the stiffener section modulus and
	# effective area need to be reduced to account for the
	# interaction of the shear with the moment and axial force in
	# the stiffener.
	#
	_Anet = _Ae
	#
	self.VRd = _Anet * (self.fy / (self.GammaM * math.sqrt(3.0)))
	#
	self.msg = self.msg + '\n' + 'VRd =' + str( self.VRd)
	#
	#_URv = 2 * self.VSd / self.VRd
	#self.msg = self.msg + '\n' + 'URv =' + str( _URv)
	#
	#
    #
    #
    # Section 8
    def BucklingOfGirders(self):
	#
	self.msg = self.msg + '\n' + ''
	self.msg = self.msg + '\n' + 'Section 8'
	#
	self.msg = self.msg + '\n' + ''
	self.msg = self.msg + '\n' + 'Girder Section property without associated plate'
	#
	_Asf, _AG,_Zpf, _Is, _IzG, _Wes, _Wep, _fyG = PlateStiffener(0.0, 0.0,
                                                                     self.hwG, self.twG, 
                                                                     self.bG, self.tfG, 
                                                                     self.LG,
                                                                     self.fyG, self.fyG, self.fyG,
                                                                     self.EG)
	#
	_AfG = self.bG * self.tfG
	_AwG = self.hwG * self.twG
	#
	self.msg = self.msg + '\n' + ''
	self.msg = self.msg + '\n' + 'fyG =' + str(_fyG)
	self.msg = self.msg + '\n' + 'Total Cross Sectional Area =' + str(  _Asf)
	self.msg = self.msg + '\n' + 'Cross Sectional Area of stiffener =' + str( _AG)
	self.msg = self.msg + '\n' + 'Neautral Axis Loacation =' + str( _Zpf)
	self.msg = self.msg + '\n' + 'Moment of Inertia IG =' + str( _Is) + ' ' + str( _IzG)
	self.msg = self.msg + '\n' + ''
	#
	#
	# 8.4 Effective widths of girders  
	# -------------------------------
	#
	self.msg = self.msg + '\n' + '8.4 Effective widths of girders '
	#
	# 8.4.1 General
	# -------------
	#
	# For the determination of the effective width the designer is
	# given two options denoted method 1 and method 2. These
	# methods are described in sec. 8.4.2 and 8.4.3 respectively:
	#
	#
	# 8.4.2 Method 1
	# --------------
	#
	self.msg = self.msg + '\n' + ''
	self.msg = self.msg + '\n' + '8.4.2 Method 1'
	self.msg = self.msg + '\n' + '--------------'
	#
	# Calculation of the girder by assuming that the stiffened plate
	# is effective against transverse compression (sy) stresses. See
	# also Commentary Chapter 10 and sec. 7.1.
	#
	# In this method the effective width may be calculated as:
	#
	# Cxs is found from eq. (7.14).
	# (8.20)
	_fkx = self.Cxs * self.fy
	#
	self.msg = self.msg + '\n' + 'fkx =' + str( _fkx)
	#
	# (8.19)
	_CxG_1 = math.sqrt(1 - (self.SigmaxSd / _fkx)**2)
	#
	self.msg = self.msg + '\n' + 'CxG 1 =' + str( _CxG_1)
	#
	# If the Sigmay stress in the plate is partly or complete in
	# compression CyG may be found from eq. (7.16).
	#
	# (8.22)
	if self.Stress == 'COMPRESSION':
	    #
	    #_CyG_1 = math.sqrt(1.0 - 3 * (self.TauSd / self.fy)**2)
	    _CyG_1 = self.Cys
	    #
	#
	# If the Sigmay stress in the girder is in tension due to the combined
	# girder axial force and bending moment over the total span of
	# the girder CyG may be taken as:
	#
	# (8.21)
	else:
	    #
	    _CyG_1 = 1.0
	    #
	#
	self.msg = self.msg + '\n' + 'CyG 1 =' + str( _CyG_1)
	#
	#
	# Le should not be taken larger than 0.3 LG for continuous
	# girders and 0.4 LG for simple supported girders when
	# calculating section modules Wep and WeG.
	#
	#
	# 8.4.3 Method 2
	# --------------
	#
	self.msg = self.msg + '\n' + ''
	self.msg = self.msg + '\n' + '8.4.3 Method 2'
	self.msg = self.msg + '\n' + '--------------'
	#
	#
	# Calculation of the girder by assuming that the stiffened plate
	# is not effective against transverse compression stresses (sy).
	# See also Commentary Chapter 10 and Sec. 7.1.
	#
	# In this case the plate and stiffener can be checked with sy
	# stresses equal to zero.
	#
	# In method 2 the effective width for the girder should be
	# calculated as if the stiffener was removed.
	# then:
	# (8.23)
	_CxG_2 = math.sqrt(1.0 - (self.SigmaxSd / self.fy)**2)
	#
	self.msg = self.msg + '\n' + 'CxG 2 ='  + str( _CxG_2 )
	#
	# where
	# SigmaxSd is based on total plate and stiffener area in x-direction.
	#
	# (8.25)
	_LambdaG = 0.525 * (self.L / self.t) * math.sqrt(self.fy / self.E)
	#
	self.msg = self.msg + '\n' + 'LambdaG =' + str( _LambdaG )
	#
	#
	# (8.24)
	if _LambdaG > 0.673 :
	    #
	    _CyG_2 = (_LambdaG - 0.22) / _LambdaG**2
	    #
	#
	# _LambdaG <= 0.673 
	else:
	    #
	    _CyG_2 = 1.0
	    #
	#
	#
	self.msg = self.msg + '\n' + 'CyG 2 =' + str( _CyG_2 )
	#
	#
	# Select correct CxG and CyG
	#
	self.msg = self.msg + '\n' + 'Effective PL Sigmay' + str( self.EffectivePLSigmay)
	#
	if self.EffectivePLSigmay == 'Y':
	    _CxG = _CxG_1
	    _CyG = _CyG_1
	#
	else:
	    _CxG = _CxG_2
	    _CyG = _CyG_2
	#
	self.msg = self.msg + '\n' + 'CxG =' + str(_CxG)
	self.msg = self.msg + '\n' + 'CyG =' + str( _CyG)
	#
	# (8.26)
	_CtauG = math.sqrt(1.0 - 3.0 * (self.TauSd / self.fy)**2)
	#
	self.msg = self.msg + '\n' + 'CtauG =' + str( _CtauG)
	#
	# The effective width for the plate of the girder is taken equal
	# to:
	#
	# (8.18)
	if self.GirderSupport == 'C':
	    # 
	    _Le = min((_CxG * _CyG * _CtauG * self.L), 0.30 * self.LG)
	#
	else:
	    #
	    _Le = min((_CxG * _CyG * _CtauG * self.L), 0.40 * self.LG)
	    #
	#
	self.msg = self.msg + '\n' + '_Le =' + str( _Le)
	#
	#
	#
	# Calculate Girder Section Property with Effective Plate
	# -----------------
	#
	self.msg = self.msg + '\n' + 'Section Property of Stiffener with Efective Plate'
	#
	_AGe, _ApG, _ZpG, _IGe, _IpG, _WeG, _WepG, _fy = PlateStiffener(_Le, self.t,
                                                                        self.hwG, self.twG, 
                                                                    self.bG, self.tfG, 
                                                                    self.LG,
                                                                    self.fyG, self.fyG, self.fyG,
                                                                    self.EG)
	#
	#
	_ZtG = (0.50*self.t + self.hwG + self.tfG - _ZpG)
	_WG = min(_WeG, _WepG)
	#
	self.msg = self.msg + '\n' + ''
	self.msg = self.msg + '\n' + 'fy ='  + str(_fy)
	self.msg = self.msg + '\n' + 'Total Cross Sectional Area AGe =' + str( _AGe)
	self.msg = self.msg + '\n' + 'Neautral Axis Loacation ZpG ='+ str(_ZpG)
	self.msg = self.msg + '\n' + 'Neautral Axis Loacation ZtG =' +  str( _ZtG)
	self.msg = self.msg + '\n' + 'Moment of Inertia IGe ='+ str( _IGe)
	self.msg = self.msg + '\n' + 'Torsional Moment of Inertia IpG ='+ str( _IpG)
	self.msg = self.msg + '\n' + 'Section Modulus = WeG' + str( _WeG)
	self.msg = self.msg + '\n' + 'Section Modulus = WepG' + str( _WepG)
	self.msg = self.msg + '\n' + 'Total Section Modulus = WG' + str( _WG)
	self.msg = self.msg + '\n' + ''
	#
	#
	#
	# 8.2 Girder forces
	# -----------------
	#
	# print('')
	# print('8.2 Girder forces')
	#
	#
	# The axial force should be taken as:
	# (8.1)
	_NySd = self.SigmaySd * (self.L * self.t + _AG)
	self.msg = self.msg + '\n' + 'NySd ='+ str( _NySd)
	#
	# The lateral line load should be taken as:
	# (8.2)
	#_qSd = (self.PSd + _Po) * self.L
	#
	# where
	# pSd = design lateral pressure
	# p0 = equivalent lateral pressure
	# AG = cross sectional area of girder
	#
	# The calculation of the additional equivalent lateral pressure
	# due to longitudinal compression stresses and shear shall be
	# calculated as follows:
	#
	# LP = length of panel
	# hwG = web height of girder
	# As = cross sectional area of stiffener
	# LG = girder span
	# s = stiffener spacing
	# Is = moment of inertia of stiffener with full plate width
	#
	# For linear variation of SigmaxSd, the maximum value within
	# 0.25LG to each side of the midpoint of the span may be used
	#
	# TauSd should correspond to the average shear flow over the
	# panel
	#
	#
	_Taucel = ((18.0 * self.E / (self.t * self.L**2)) 
                   * math.pow(self.t * self.Is / self.S, 0.75))
	#
	self.msg = self.msg + '\n' + 'Taucel =' + str( _Taucel)
	#
	_Tauceg = _Taucel * self.L**2 / self.Lp**2
	#
	self.msg = self.msg + '\n' + 'Tauceg =' + str( _Tauceg)
	#
	# Taucrg = critical shear stress of panel with girders removed,
	#          calculated from eq.(8.6) with Lambdatau calculated 
	#          using:
	#          Tauce =  Tauceg If the stiffener is not continuous
	#                   through the girder tcrg = 0.
	#
	#
	# Taucrl = critical shear stress of panel between girders 
	#          calculated from eq. (8.6) with Lambdatau calculated 
	#          using Tauce = Taucel
	#
	#
	_LambdaTau = lambda x: math.sqrt(0.60*self.fy / x)
	#
	#
	_LambdaTau_ce = {}
	#
	_LambdaTau_ce[0] = _LambdaTau(_Taucel)
	self.msg = self.msg + '\n' + 'LambdaTau cel =' + str( _LambdaTau_ce[0])
	#
	#
	_LambdaTau_ce[1] = _LambdaTau(_Tauceg)
	self.msg = self.msg + '\n' + 'LambdaTau ceg =' + str( _LambdaTau_ce[1])
	#
	#
	_Taucr = {}
	#
	_Taucr[0] = self.TauCR( _LambdaTau_ce[0])
	self.msg = self.msg + '\n' + 'Taucr cel =' + str( _Taucr[0])
	#
	_Taucr[1] = self.TauCR(_LambdaTau_ce[1])
	self.msg = self.msg + '\n' + 'Taucr ceg =' + str( _Taucr[1])
	#
	#
	# (8.17)
	_iGe = math.sqrt(_IGe / _AGe)
	#
	self.msg = self.msg + '\n' + 'iGe =' + str( _iGe)
	#
	# (8.11)
	_fEG = math.pi**2 * self.E * (_iGe / self.LGk)**2
	#
	self.msg = self.msg + '\n' + 'fEG =' + str( _fEG)
	#
	# fEG is given in eq (8.11)
	#
	_LambdaG = math.sqrt(self.fy / _fEG)
	#
	self.msg = self.msg + '\n' + 'Lambda G =' + str( _LambdaG)
	# 
	#
	_Q1 = max(_LambdaG - 0.20, 0)
	_Q = min(_Q1 , 1.0)
	#
	self.msg = self.msg + '\n' + 'Q =' + str(_Q)
	#
	#
	if self.Stiffener == 'S':
	    _TaucrG = 0
	#
	else:
	    _TaucrG = _Taucr[1]
	#
	self.msg = self.msg + '\n' + '_TaucrG =' + str( _TaucrG)
	#
	#
	print ('1')
	if self.TauSd > _TaucrG :
	    # (8.4)
	    _CG = (_Q * (7.0 - 5.0 * (self.S / self.L)**2) 
                   * ((self.TauSd - _TaucrG) / self._Taucrl)**2)
	    #
	#
	# self.TauSd <= _Taucrg
	else:
	    # (8.5)
	    _CG = 0
	    #
	#
	self.msg = self.msg + '\n' + 'C =' + str( _CG)
	#
	#
	# For tension in the x-direction:
	if self.SigmaxSd < 0:
	    # (8.7)
			#_Po = (((0.40 * (self.t + _As / self.S)) 
	    _Po = (((0.40 * (self.t + self.As / self.S)) 
                    / (self.hwG * (1.0 - self.S / self.LG))) 
                   * (self.fy / self.E) * (self.LG / self.L)**2 
                   * (_CG * self.TauSd))
	#
	# For compression in the x-direction:
	else:
	    #
	    # (8.3)
	    _Po1 = (((0.40 * (self.t + self.As / self.S)) 
                     / (self.hwG * (1.0 - self.S / self.LG))) 
                    * (self.fy / self.E) * (self.LG / self.L)**2 
                    * (self.SigmaxSd + _CG * self.TauSd))
	    #
	    self.msg = self.msg + '\n' + 'Po1 =' + str( _Po1) + ' C*Tsd = ' +str(( _CG * self.TauSd))
	    # But not less than
	    _Po2 = (0.020 * ((self.t + (self.As / self.S)) / self.L) 
                    * (self.SigmaxSd + _CG * self.TauSd))
	    #
	    self.msg = self.msg + '\n' + 'Po2 =' + str( _Po2)
	    #
	    _Po = max(_Po1, _Po2)
	    #
	#
	# Po shall be applied in the direction of external pressure
	# Psd. For situtations where Psd is less than Po. the girder
	# need to be checked for Po applied in both directions
	# (i.e. at plate side and stiffener side)
	#
	self.msg = self.msg + '\n' + 'SigmaxSd =' + str( self.SigmaxSd )
	self.msg = self.msg + '\n' + 'PSd =' + str(self.PSd)
	self.msg = self.msg + '\n' + 'p0 =' + str(_Po)
	#
	#
	# Equivalent lateral line load on plate side
	#
	# (8.2)
	#
	if self.PSd > 0:
	    _qSdp = (self.PSd + _Po) * self.L
	#
	else:
	    #
	    if _Po < abs(self.PSd):
		_qSdp = 0.0
	    #
	    else:
		_qSdp = (self.PSd + _Po) * self.L
	    #
	#
	self.msg = self.msg + '\n' + '_qSdp =' + str( _qSdp)
	#
	# Equivalent lateral line load on stiffener side
	#
	#
	if self.PSd > 0:
	    #
	    _qSds = ( _Po - self.PSd ) * self.L
	#
	else:
	    #
	    if _Po < self.PSd:
		_qSds = 0.0
	    #
	    else:
		_qSds = ( _Po - self.PSd ) * self.L
	    #
	#
	self.msg = self.msg + '\n' + '_qSds =' + str( _qSds)
	#
	#
	#
	#
	# 8.5 Torsional buckling of girders
	# ---------------------------------
	#
	self.msg = self.msg + '\n' + ''
	# print('8.5 Torsional buckling of girders')
	#
	#
	# SigmaySd = compressive stress in the free flange
	# (8.32)
	self.PSd = 0.020 * self.SigmaySd * (_AfG + _AwG / 3.0)
	#
	self.msg = self.msg + '\n' + 'PSd =' + str( self.PSd )
	#
	#
	# Torsional buckling need not to be considered if tripping
	# brackets are provided so that the laterally unsupported 
	# length LGT, does not exceed the value LGT0 defined by:
	# where
	# b = flange width
	# C = 0.55 for symmetric flanges
	#     1.10 for one sided flanges
	#
	if self.efG == 0:
	    _CGTO = 0.55
	#
	else:
	    _CGTO = 1.10
	#
	self.msg = self.msg + '\n' + 'C GTO =' + str( _CGTO)
	#
	# (8.31)
	_LGTO = (self.bG * _CGTO *
                 ( math.sqrt(self.EG * _AfG 
                             / (self.fyG * (_AfG + _AwG / 3.0)))))
	#
	self.msg = self.msg + '\n' + 'L GTO =' + str(_LGTO)
	#
	self.msg = self.msg + '\n' + ''
	if _LGTO < self.LGt:
	    self.msg = self.msg + '\n' + 'Girder torsional buckling check need to be considered'
	    pass
	#
	else:
	    self.msg = self.msg + '\n' + 'Girder torsional buckling not required'
	    pass
	#
	# The torsional buckling strength of girders may be determined
	# as:
	#
	# LGT = distance between lateral supports
	# Af, Aw = cross sectional area of flange and web of girder
	# Iz = moment of inertia of girder (exclusive of plate flange)
	#      about the neutral axis perpendicular to the plate
	#
	_LGt1 = min(0.40*self.LG, self.LGt)
	#
	self.msg = self.msg + '\n' + 'LGt 1 ' + str( _LGt1)
	#
	#
	_LGt2 = min(0.80*self.LG, self.LGt)
	#
	self.msg = self.msg + '\n' + 'LGt 2 ' +str(_LGt2)
	#
	#
	self.msg = self.msg + '\n' + ''
	#
	# (8.30)
	fETG = lambda x: ((math.pi**2 * self.EG * _IzG) 
                          / ((_AfG + _AwG / 3.0) * x**2))
	#
	_fETG = {}
	#
	_fETG[0] = fETG(self.LGt)
	self.msg = self.msg + '\n' + 'fETG-LGt ='+str( _fETG[0]) + ' ' + str(self.LGt)
	#
	_fETG[1] = fETG(_LGt1)
	self.msg = self.msg + '\n' + 'fETG-_LGt1 =' + str( _fETG[1]) + ' ' + str(_LGt1)
	#
	_fETG[2] = fETG(_LGt2)
	self.msg = self.msg + '\n' + 'fETG-_LGt2 =' + str( _fETG[2]) + ' ' + str(_LGt1)
	#
	#
	#
	# print('')
	#
	# (8.28)
	LambdaTG = lambda x: math.sqrt(self.fy / x)
	#
	_LambdaTG = {}
	#
	_LambdaTG[0] = LambdaTG(_fETG[0])
	self.msg = self.msg + '\n' + 'LambdaTG - LGt ='+ str( _LambdaTG[0])
	#
	_LambdaTG[1] = LambdaTG(_fETG[1])
	self.msg = self.msg + '\n' + 'LambdaTG - LGt1 =' + str( _LambdaTG[1])
	#
	_LambdaTG[2] = LambdaTG(_fETG[2])
	self.msg = self.msg + '\n' + 'LambdaTG - LGt2 =' + str( _LambdaTG[2])
	#
	#
	# (8.29)
	MuTG = lambda x: 0.35 * (x - 0.60)
	#
	# print('')
	#
	_MuTG = {}
	#
	_MuTG[0] = MuTG(_LambdaTG[0])
	self.msg = self.msg + '\n' + 'MuTG - LGt =' + str( _MuTG[0])
	#
	_MuTG[1] = MuTG(_LambdaTG[1])
	self.msg = self.msg + '\n' + 'MuTG - LGt1 ='+ str( _MuTG[1])
	#
	_MuTG[2] = MuTG(_LambdaTG[2])
	self.msg = self.msg + '\n' + 'MuTG - LGt2 =' + str( _MuTG[2])
	#
	#
	self.msg = self.msg + '\n' + ''
	#
	_fTG = {}
	#
	_fTG[0] = self.fTG( _LambdaTG[0], _MuTG[0])
	self.msg = self.msg + '\n' + 'fTG - LGt =' + str( _fTG[0])
	#
	_fTG[1] = self.fTG( _LambdaTG[1], _MuTG[1])
	self.msg = self.msg + '\n' + 'fTG - LGt 1 =' + str( _fTG[1])
	#
	_fTG[2] = self.fTG( _LambdaTG[2], _MuTG[2])
	self.msg = self.msg + '\n' + 'fTG - LGt 2 =' + str( _fTG[2])
	#
	#
	#
	# Tripping brackets are to be designed for a lateral force PSd,
	# which may be taken equal to (see Figure 8-2 ):
	#
	#
	#
	# 8.3 Resistance parameters for girders
	# -------------------------------------
	#
	self.msg = self.msg + '\n' + ''
	self.msg = self.msg + '\n' + '8.3 Resistance parameters for girders'
	self.msg = self.msg + '\n' + ''
	#
	#
	#
	_fkG = {}
	_frG = {}
	#
	_frG[0], _LambdaG0, _MuG0, _fkG[0] = self.fk( 'S', _fTG[0], _fEG, _iGe, _ZtG, _ZpG)
	self.msg = self.msg + '\n' + 'fk - S =' + str( _fkG[0])
	#
	_frG[1], _LambdaG1, _MuG1, _fkG[1] = self.fk( 'P', _fTG[0], _fEG, _iGe, _ZtG, _ZpG)
	self.msg = self.msg + '\n' + 'fk - P =' + str( _fkG[1])
	#
	_frG[2], _LambdaG1, _MuG1, _fkG[2] = self.fk( 'S', _fTG[1], _fEG, _iGe, _ZtG, _ZpG)
	self.msg = self.msg + '\n' + 'fk 1 - S =' + str( _fkG[2])
	#
	_frG[3], _LambdaG1, _MuG1, _fkG[3] = self.fk( 'S', _fTG[2], _fEG, _iGe, _ZtG, _ZpG)
	self.msg = self.msg + '\n' + 'fk 2 - S =' + str( _fkG[3])
	self.msg = self.msg + '\n' + ''
	#
	#
	# The resistance of girders may be determined by the
	# interaction formulas in sec. 7.7 using the following
	# resistance
	#
	# (8.8)
	_NRd = (_AG + _Le * self.t) * (self.fy / self.GammaM)
	#
	self.msg = self.msg + '\n' + 'NRd =' + str( _NRd)
	#
	# fk is calculated from sec 7.5 using Mu according to eq (7.26)
	# (8.9)
	_NksRd = (_AG + _Le * self.t) * (_fkG[0] / self.GammaM)
	self.msg = self.msg + '\n' + 'NksRd =' + str( _NksRd)
	#
	# 
	# fk is calculated from sec. 7.5 using Mu according to eq.
	# (7.25) using:
	# fr = fy for check at plate side
	# fr = fTG for check at girder flange side
	#
	# (8.10)
	_NkpRd = (_AG + _Le * self.t) * (_fkG[1] / self.GammaM)
	self.msg = self.msg + '\n' + 'NkpRd =' + str( _NkpRd)
	#
	#
	# (8.12)
	_Ms1Rd = _WeG * (_frG[2] / self.GammaM)
	self.msg = self.msg + '\n' + 'Ms1Rd =' + str( _Ms1Rd)
	#
	# (8.13)
	_Ms2Rd = _WeG * (_frG[3] / self.GammaM)
	self.msg = self.msg + '\n' + 'Ms2Rd =' + str( _Ms2Rd)
	#
	# (8.14)
	_MstRd = _WeG * (self.fyG / self.GammaM)
	self.msg = self.msg + '\n' + 'MstRd =' + str( _MstRd)
	#
	# (8.15)
	_MpRd = _WepG * (self.fyG / self.GammaM)
	self.msg = self.msg + '\n' + 'MpRd =' + str( _MpRd)
	#
	# (8.16)
	_NE = (math.pi**2 * self.E * _AGe / (self.LGk / _iGe)**2)
	self.msg = self.msg + '\n' + 'NE =' + str( _NE)
	#
	#
	#
	# 8.1 General
	# -----------
	#
	self.msg = self.msg + '\n' + ''
	self.msg = self.msg + '\n' + '8.1 General'
	self.msg = self.msg + '\n' + ''
	#
	# The check for girders is similar to the check for stiffeners of
	# stiffened plates in equations (7.50) to (7.57) or (7.59) to
	# (7.64) for continuous or sniped girders, respectively. Forces
	# shall be calculated according to sec. 8.2 and cross section
	# properties according to 8.4. Girder resistance should be
	# found from sec. 8.3. Torsional buckling of girders may be
	# assessed according to sec. 8.5.
	#
	#
	self.CombinedUnityCheck( self.GirderSupport, _NySd, self.LG, _qSdp, _qSds, _NksRd, _NkpRd, _MpRd, _MstRd, _Ms1Rd, _Ms2Rd, _NE, _NRd,  _z = 0.0 , _u = 0)
	self.ur_combined_Girder = self.ur_combined
	#
	#
	#
    #
    #
    # Section 9
    def LocalBuckling(self):
	#
	# 9 Local buckling of stiffeners, girders and
	# brackets
	# -------------------------------------------
	#
	# 9.1 Local buckling of stiffeners and girders
	# --------------------------------------------
	#
	# 9.1.1 General
	# -------------
	#
	# The methodology given in Chapter 7 and Chapter 8 is only
	# valid for webs and flanges that satisfy the the following
	# requirements or fulfils requirements to cross section type
	# III defined in Appendix A of DNV-OS-C101
	#
	# In lieu of more refined analysis such as in Chapter 7, web
	# stiffeners should satisfy the requirements given in sec. 9.1.2
	# and sec. 9.1.3. 
	#
	# Flange outstand for T or L stiffeners or girders should
	# satisfy:
	#
	_Epsilon = math.sqrt(235.0 / self.fy)
	#
	# For definition of c see Figure 7-3
	#
	# for welded sections
	# (9.1)
	if self.SectionType == 'WELDED':
	    #
	    if self.c > (14.0 * self.tf * _Epsilon):
		self.msg = self.msg + '\n' + '===> Fail'

	#
	# for rolled sections
	else:
	    #
	    if self.c > (15.0 * self.tf * _Epsilon):
		self.msg = self.msg + '\n' + '===> Fail'
	#
	#
	# Web of stiffeners and girders should satisfy:
	# (9.2)
	if self.hw > (42.0 * self.tw * _Epsilon):
	    #
	    self.msg = self.msg + '\n' + '===> Fail'
	#
	#
	# 9.1.2 Transverse web stiffeners:
	# -------------------------------
	#
	# Is = moment of inertia of web stiffener with full web plate
	#      flange s
	# Lt = length of transverse web stiffener
	# S = distance between transverse web stiffeners
	#
	_Is93 = ((0.30 * self.Lt * self.S**2 * self.tw) 
                 * ((2.50 * self.Lt / self.S) - (2.0 * self.S / self.Lt))
                 * (self.fy / self.E))
	#
	# (9.3)
	if self.Is < _Is93 :
	    #
	    self.msg = self.msg + '\n' + '===> Fail'
	    #
	#
	#
	# 9.1.3 Longitudinal web stiffener:
	# --------------------------------
	#
	# Is = moment of inertia of web stiffener with full web plate
	#      flange s.
	# As = cross sectional area of web stiffener exclusive web
	#      plating.
	# Ll = length of longitudinal web stiffener
	#  S = distance between longitudinal web stiffeners
	#
	_Is94 = ((0.25 * self.Ll**2) * (_As + self.S * self.tw) 
                 * (self.fy / self.E))
	#
	# (9.4)
	if self.Is < _Is94 :
	    #
	    self.msg = self.msg + '\n' + '===> Fail'
	    #
	#
	#
	# 9.2 Buckling of brackets
	# Brackets should be stiffened in such a way that:
	#
	# tb = plate thickness of bracket.
	# Stiffeners as required in eq. (9.6) or eq. (9.7) may be
	# designed in accordance with Chapter 7. See Figure 9-3
	#
	# (9.5)
	if self.BracketType == 'FREE':
	    #
	    if self.d0 > (0.70 * self.tb * math.sqrt(self.E / self.fy)):
		#
		self.msg = self.msg + '\n' + '===> Fail'
	#
	# (9.6)
	elif self.BracketType == 'SINGLE':
	    #
	    if self.d1 > (1.650 * self.tb * math.sqrt(self.E / self.fy)):
		#
		self.msg = self.msg + '\n' + '===> Fail'
	#
	# (9.7)
	elif self.BracketType == 'DOUBLE':
	    #
	    if self.d2 > (1.350 * self.tb * math.sqrt(self.E / self.fy)):
		#
		self.msg = self.msg + '\n' + '===> Fail'
		#
	#
	#
	else:
	    #
	    self.msg = self.msg + '\n' + 'Bracket type not supported'
	    #
	#
	#
    #
    #
    #
    #
#
#
#     Stiffener is at top/bottom of panel as defined by ILRFLG below
#
#              6                5                  4
#           x  |=================|=================|  x
#           :  |                 |                 |  :
#           :  |                 |                 |  :
#    Plate  :  |=================|=================|  x
#    Length :  |                 |                 |  : Stiffener Spacing
#   Vertical:  |                 |                 |  :       (Ls)
#     (LG)  :  |=================|=================|  x
#           : Y|                 |                 |  :
#           : ^|                 |                 |  :
#           x  |=================|=================|  x
#              1  > X            2                  3
#              x-----------------x-----------------x 
#              :      (L1)       :      (L2)       : Girder Spacing
#              :                 :                 :
#              x-----------------x                 :
#              :       (L)       : Stiffener Length
#              :                 :                 :
#              x-----------------------------------x Panel Length Htal
#                              (Lp)
#
#     Note:
#           b1,b2   : Sub-panel widths for panels adjoining stiffener
#           tf1, tf2: Panel thickness for panels adjoining stiffener
#           PANEL1  : Flags if panels present on each side of stiffener
#           PANEL2
#
#
class CodeCheckPanel(DNVRPC201):
    #
    """
    This programm checks stiffener buckling strength of stiffened panels 
    according to DNV recommended practice "Buckling Strength of Plated 
    Structures", DNV-RP-C201, October 2002.with amendments in October 2008
    Note: compressive stresses are positive,  
    tensile stresses are negtive.
    """
    #
    def __init__(self):
	#
	self.DesignCode = 'DNV'
	#
	#
	self.Stress = 'TENSION'
	self.PlateDescription = 'N/A'
	self.LoadDescription = 'N/A'
	self.BucklingCheck = 'YES'
	#
	# Material Factor
	self.GammaM = 1.150
	#
	# Longitudinal stress
	self.Sigmax1Sd = 0
	self.Sigmax2Sd = 0
	#
	# Tranverse stress 
	self.Sigmay1Sd = 0
	self.Sigmay2Sd = 0
	#
	# Shear stress 
	self.TauSd = 0
	#
	# Lateral Pressure 
	self.PSd = 0
	#
	# lenght to reference point
	#
	self.L = 0
	self.L1 = 0
	self.S = 0
	self.S1 = 0
	self.LG = 0
	# 
	# for simplification Z = 0
	self.Z = 0
	#
	self.ur = 0.0
	self.msg = ''
	#
    #
    #
    def GeneralData(self, Stiffener = 'C', Girder = 'C', EPLSy = 'N', GM = 1.15):
	#
	#
	# Stiffener can be :
	# C for continuos stiffener
	# S for simply supported stiffener (with sniped ends)
	# 
	self.Stiffener = Stiffener.upper()
	#
	# Girder can be :
	# C for continous girder 
	# S for simply supported girder (with sniped ends)
	self.Girder = Girder.upper()
	#
	# Y if assuming that the stiffened plate is effective
	#   against transverse compression stress Sigmay
	# N if otherwise
	# (refer section 8.4.2 & 8.4.3)
	self.EffectivePLSigmay = EPLSy.upper()
	#
	# Material Factor
	self.GammaM = GM
	#
	#
    #
    #
    # Plate Pannel Section
    #
    def PlateGeometry(self, thk, LP, SP):
	#
	# Plate thickness
	self.t = thk
	# Plate Length horizontal
	self.Lp = LP
	# Plate Length Vertical
	self.LG = SP
	#
	self.PlateDescription = 'PLATE + ' + str(self.PlateDescription)
	#
    #
    #
    def PlateMaterial(self, Fyp = 265.0, E = 210000, Nu = 0.3):
	#
	# Young's Modulus
	self.E = E
	# Yield Strength
	self.fyp = Fyp
	# Poisson's ratio
	self.Poisson = Nu
	#
    #
    #
    # Stiffener Section
    #
    def StiffenerCrossSection(self, StifType, hw, b, tw, tf, S1, S2 = 0):
	#
	#
	# Type for stiffener can be :
	# I for I-Section
	# L for angle
	# T for T-Section
	# F for flat bar
	# B for Box-Section
	self.StiffenerType = StifType.upper()
	#
	# Stiffener Section Detail
	#
	# height of stiffener web
	self.hw = hw
	# web thickness
	self.tw = tw
	#
	# width of flange (top)
	self.bf = b
	# flange thickness
	self.tf = tf
	# self.ef = self.tw /2
	#
	# Stiffener Spacing (mm)
	# width of flange (bottom)
	self.s1 = S1
	self.s2 = S2
	#
	#
	self.PlateDescription = 'STIFFENER + ' + str(self.PlateDescription)
	#
    #
    #
    def StiffenerBucklingLength(self, LS, LT):
	#
	# Stiffener Length
	self.L = LS
	# Stiffener torsional buckling length
	self.Lt = LT
	#
	#
    #
    #
    def StiffenerMaterial(self, FyS = 265.0, ES = 210000, Nu = 0.3):
	#
	# Yield Strength (N/mm2)
	self.fyS = FyS
	# Young's Modulus
	self.ES = ES
	# Poisson's ratio
	self.PoissonS = Nu
	#
	#
    #
    #
    # Girder Section
    #
    def GirderCrossSection(self, GType , hwG, twG, bG, tfG, L1 = 0, L2 = 0):
	#
	#
	# Type for stiffener can be :
	# I for I-Section
	# L for angle
	# T for T-Section
	# F for flat bar
	# B for Box-Section
	self.GirderType = GType.upper()
	#
	# height of stiffener web
	self. hwG = hwG
	# web thickness
	self.twG = twG
	# width of flange
	self.bG = bG
	# flange thickness
	self.tfG = tfG
	#
	# Girder Spacing (mm)
	self.L1 = L1
	self.L2 = L2
	#
	self.PlateDescription = 'GIRDER + ' + str(self.PlateDescription)
	#
    #
    #
    def GirderBucklingLength(self, LG, LGK, LGT):
	#
	# Girder Length
	self.LG = LG
	#
	# Girder Overal Buckling Length
	self.LGk = LGK
	# Girder Torsional Buckling Length
	# (Length between tripping brackets)
	self.LGt = LGT
	#
    #
    #
    def GirderMaterial(self, FyG = 265.0, EG = 210000, Nu = 0.3):
	#
	# Yield Strength (N/mm2)
	self.fyG = FyG
	# Young's Modulus
	self.EG = EG
	# Poisson's ratio
	self.PoissonG = Nu
	#
    #
    #
    # Web Stiffener Section
    #
    def WebStiffener(self, WSType, Tws, Lws, SWs):
	# 
	# Web Stiffener Type:
	# TRANVERSE
	# LONGITUDINAL
	#
	self.WebStiffenerType = WSType
	#
	# thickness of web stiffener
	self.tws = Tws
	#
	# Length of tranverse web stiffener
	self.Lws = Lws
	#
	# Distance between transverse web stiffener
	self.Sws = SWs
	#
	self.PlateDescription = 'WEB_STIFFENER + ' + str(self.PlateDescription)
	#
    #
    #
    def WebStiffenerMaterial(self, FyWS = 265.0, EWS = 210000, Nu = 0.3):
	#
	# Yield Strength (N/mm2)
	self.fyWs = FyWS
	# Young's Modulus
	self.EWs = EWS
	# Poisson's ratio
	self.PoissonWs = Nu
	#
    #
    #
    # Bracket Section
    #
    def BracketGeometry(self, BType, Tb, D1, D2):
	#
	# Bracket type:
	# FREE-EDGE
	# SINGLE
	# DOUBLE
	#
	self.BracketType = BType
	#
	# plate thickness of the bracket
	#
	self.tb = Tb
	#
	# see Figure 9-3
	# distance d0 or d1
	self.d1 = D1
	# distance d2 for double brackets
	self.d2 =D2
	#
	self.PlateDescription = 'BRACKET + ' + str(self.PlateDescription)
	#
    #
    #
    def BracketMaterial(self, FyB = 265.0, EB = 210000, Nu = 0.3):
	#
	# Yield Strength (N/mm2)
	self.fyB = FyB
	# Young's Modulus
	self.EB = EB
	# Poisson's ratio
	self.PoissonB = Nu
	#
    #
    #
    #
    # Stress / Pressure Section
    #
    #
    def LongitudinalStress(self, Sigmax1Sd, Sigmax2Sd = 0, S1 = 0):
	#
	# Longitudinal stress larger
	# (compression positive)
	self.Sigmax1Sd = Sigmax1Sd
	#
	# Longitudinal stress smaller
	# (compression positive)
	self.Sigmax2Sd = Sigmax2Sd
	#
	# lenght to reference point
	self.S1 = S1
	#
	self.LoadDescription = 'LONGITUDINAL + ' + str(self.LoadDescription)
	#
	#
    #
    #
    def TransverseStress(self, Sigmay1Sd, Sigmay2Sd = 0, L1 = 0):
	#
	# Tranverse stress larger 
	# (compression positive)
	self.Sigmay1Sd = Sigmay1Sd
	#
	# Tranverse stress smaller 
	# (compression positive)
	self.Sigmay2Sd = Sigmay2Sd
	#
	# lenght to reference point
	self.L1 = L1
	#
	self.LoadDescription = 'TRANSVERSE + ' + str(self.LoadDescription)
	#
    #
    #
    def ShearStress(self, TauSd):
	#
	# Shear stress 
	self.TauSd = TauSd
	#
	self.LoadDescription = 'SHEAR + ' + str(self.LoadDescription)
	#
    #
    #
    def LateralPressure(self, PSd):
	# 
	# Lateral Pressure 
	# *Note - Lateral pressure should be
	#  input as:
	# Positive for pressure on Plate Side
	# Negative for pressure on stiffener side
	#
	#self.PressureSide = PSide
	#
	self.PSd = PSd
	#
	#
	self.LoadDescription = 'LATERAL-PRESSURE + ' + str(self.LoadDescription)
	#
    #
    #
    def PrintResults (self):
	#
	#
	# Plate Description Selection
	#
	self.msg = self.msg + '\n' + 'Decription: ' + str(self.PlateDescription)
	#
	self.PlateDescription = PlateCase(self.PlateDescription)
	#
	self.msg = self.msg + '\n' + 'Plate Description: ' + str( self.PlateDescription)
	#
	#
	# Load Description Selection
	#
	self.msg = self.msg + '\n' + 'Load Description =' + str( self.LoadDescription)
	#
	self.LoadDescription = LoadCase(self.LoadDescription)
	#
	self.msg = self.msg + '\n' + 'Load Description =' + str( self.LoadDescription) + str(len(self.LoadDescription.split()))
	#
	#
	#
	#
	if 'UNSTIFFENED' in self.PlateDescription:
	    #
	    # print('Unstiffened plate')
	    #
	    self.S = self.LG
	    #
	    self.L = self.Lp
	    #
	    self.fy = self.fyp
	    #
	#
	else:
	    #  stiffener
	    self.msg = self.msg + '\n' + ''
	    self.msg = self.msg + '\n' + 'Stiffened plate'
	    #
	    # flange excentricity
	    self.ef = self.tw / 2.0
	    #
	    # plate width, stiffener spacing
	    self.S = (self.s1 + self.s2) / 2.0
	    #
	    # define mc factor
	    #
	    if self.Stiffener.upper() == 'S':
		self.mc = 8.90
	    #
	    elif self.Stiffener.upper() == 'C':
		self.mc = 13.30
	    #
	    else:
		'Stiffener type not recognized'
		exit
	    #
	    # Girder
	    if 'GIRDER' in self.PlateDescription:
		#
		# Girder excentricity
		self.efG = self.twG / 2.0
	    #
	    #
	    _b = max(self.s1 , self.s2)
	    #
	    self.Asf, self.As,self.Zpf, self.Is, _Iz, _Wes, _Wep, self.fy = PlateStiffener(_b, self.t,
                                                                                           self.hw, self.tw, 
                                                                                           self.bf, self.tf, 
                                                                                           self.L,
                                                                                           self.fyp, self.fyS, self.fyS,
                                                                                           self.ES)
	    #
	    #
	    self.msg = self.msg + '\n' + 'fy ='+str(self.fy)
	    self.msg = self.msg + '\n' + 'Total Cross Sectional Area =' + str( self.Asf)
	    self.msg = self.msg + '\n' + 'Cross Sectional Area of stiffener =' + str( self.As)
	    self.msg = self.msg + '\n' + 'Neautral Axis Loacation =' + str( self.Zpf)
	    self.msg = self.msg + '\n' + 'Moment of Inertia I =' + str( self.Is)
	    self.msg = self.msg + '\n' + ''
	    #
	#
	#
	# Calculate equivalent design transverse stress
	#
	# Longitudinal stress
	#
	#
	self.Sigmax2Sd, self.SigmaxSd, self.SigmaxType = StressDefinition(self.Sigmax1Sd, 
                                                                          self.Sigmax2Sd, 
                                                                         self.S, self.S1)
	#
	self.msg = self.msg + '\n' + ''
	self.msg = self.msg + '\n' + 'SigmaxSd =' + str( self.SigmaxSd) + ' Load Type: '  + str(self.SigmaxType)
	self.msg = self.msg + '\n' + ' '
	#
	#
	# this need to pass to 6.8
	# Transverse stress
	#
	_L1 = min(0.25*self.L, 0.50*self.S)
	self.msg = self.msg + '\n' + 'Equiv L1 =' + str( _L1)
	#
	self.Sigmay2Sd, self.SigmaySd, self.SigmayType = StressDefinition(self.Sigmay1Sd, 
                                                                          self.Sigmay2Sd, 
                                                                         self.L, _L1)
	#
	self.msg = self.msg + '\n' + 'SigmaySd =' + str( self.SigmaySd) + ' Load Type: ' + str(self.SigmayType)
	self.msg = self.msg + '\n' + ' '
	#
	#
	self.Stress = StressSign(self.SigmaxSd, self.SigmaySd)
	#
	self.msg = self.msg + '\n' + 'Stress Sign: ' + str( self.Stress)
	self.msg = self.msg + '\n' + ''
	#
	self.Epsilon = math.sqrt(235.0 / self.fy)
	#
	self.msg = self.msg + '\n' + 'Epsilon =' + str( self.Epsilon)
	self.msg = self.msg + '\n' + ''
	#
	# Shear Modulus
	self.G = self.E / (2.0 * (1 + self.Poisson))
	self.msg = self.msg + '\n' + 'Shear Modulus, G =' + str(self.G)
	#
	#
	if self.DesignCode == 'DNV':
	    #
	    #DNVRPC201.PlateBucklingCheck(self)
	    #
	    #
	    if 'UNSTIFFENED' in self.PlateDescription:
		#
		if self.S > self.L:
		    #
		    self.msg = self.msg + '\n' + 'Fail ==> s > l'
		    exit()
		    #
		#
		else:
		    #
		    #
		    if len(self.LoadDescription.split()) == 1:
			#
			self.msg = self.msg + '\n' + 'Single'
			#
			if self.SigmaxSd != 0:
			    self.msg = self.msg + '\n' + 'SigmaxSd'
			    #
			    if self.SigmaxType == 'UNIFORM':
				self.msg = self.msg + '\n' + 'UNIFORM'
				DNVRPC201.BucklingOfUnstiffenedPlatesLongitudinalUniformCompression(self)
				#
			    #
			    else:
				self.msg = self.msg + '\n' + 'Varying'
				DNVRPC201.BucklingOfUnstiffenedPlatesVaryingLongStress(self)
				#
			    #
			#
			#
			elif self.SigmaySd != 0:
			    self.msg = self.msg + '\n' + 'SigmaySd '
			    #
			    if self.SigmayType == 'UNIFORM':
				self.msg = self.msg + '\n' + 'UNIFORM'
				DNVRPC201.BucklingOfUnstiffenedPlatesTransverseCompression(self)
				#
			    #
			    else:
				self.msg = self.msg + '\n' + 'Varying'
			    #
			#
			#
			else:
			    self.msg = self.msg + '\n' + 'Shear'
			    DNVRPC201.BucklingOfUnstiffenedPlatesShear(self)
			    #
			#
			#
		    #
		    else:
			self.msg = self.msg + '\n' + 'Multiple'
			DNVRPC201.BucklingOfUnstiffenedPlatesLongitudinalUniformCompression(self)
			DNVRPC201.BucklingOfUnstiffenedPlatesTransverseCompression(self)
			DNVRPC201.BucklingOfUnstiffenedPlatesBiaxiallyLoadedShear(self)
			#
		    #
		#
	    #
	    #
	    else:
		#
		self.msg = self.msg + '\n' + 'Buckling check not necessary'
		self.BucklingCheck = 'NO'
		#
		self.msg = self.msg + '\n' + 'stiffened'
		#
		#self.LG = self.S
		# Chp 5
		DNVRPC201.LateralLoadedPlates(self)
		# Chp 6.3
		DNVRPC201.BucklingOfUnstiffenedPlatesTransverseCompression(self)
		# Chp 7
		DNVRPC201.BucklingOfStiffenedPlates(self)
		#
		#
		#
		if 'GIRDER' in self.PlateDescription:
		    #
		    self.msg = self.msg + '\n' + 'stiffened + girder'
		    #
		    # need to fix this
		    self.GirderSupport = 'C'
		    self.msg = self.msg + '\n' + 'Girder Support Type =' + str(self.GirderSupport)
		    #
		    #
		    DNVRPC201.BucklingOfGirders(self)
		    #
		    #
		    if 'WEB_STIFFENER' in self.PlateDescription:
			#
			self.msg = self.msg + '\n' + 'WEB_STIFFENER'
			#
		    #
		    #
		    if 'BRACKET' in self.PlateDescription:
			#
			self.msg = self.msg + '\n' + 'BRACKET'
		    #
		    #
		#
		#
	#self.ur = DNVRPC201.ur
	#print('DNV ' + str(DNVRPC201.ur))
	#print('self ' + str(self.ur))
	    #
	    #
	#
	#
    #
    #
    #
