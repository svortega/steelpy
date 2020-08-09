# 
# Copyright (c) 2020 fem2ufo
# 

# Python stdlib imports
import re


# package imports
from steelpy.process.units.buckingham import Number



#
class Units:
    
    def __init__(self):
        """
        Units [length, mass, time, temperature, force, pressure/stress]
        """
        #self._units = ["", "", "second", "", "", ""]
        pass
    #
    def __getattr__(self, key):
        """
        """
        #_item = unit_control.find_unit_case(key)
        #
        # mass
        if re.match(r"\b(g(ram(me|o)?(s)?)?)\b", str(key), re.IGNORECASE):
            return Number(1, dims = 'gram')
        elif re.match(r"\b(k(ilo)?g(ram(me|o)?(s)?)?)\b", str(key), re.IGNORECASE):
            return Number(1, dims = 'kilogram')
        elif re.match(r"\b((metric(\_)?)?t(on(ne|elada)?(s)?)?)\b", str(key), re.IGNORECASE):
            return Number(1, dims = 'megagram')
        elif re.match(r"\b((lb|pound)(s)?)\b", str(key), re.IGNORECASE):
            return Number(1, dims = 'pound')
        elif re.match(r"\b(slug(s)?)\b", str(key), re.IGNORECASE):
            return Number(1, dims = 'lbf * second^2 / foot')
        elif re.match(r"\b(st(one)?(s)?((\_)?weight)?)\b", str(key), re.IGNORECASE):
            return Number(14, dims = 'pound')
        elif re.match(r"\b((long|imperial|displacement)(\_)?ton(s)?)\b", str(key), re.IGNORECASE):
            return  Number(2_440, dims = 'pound')
        elif re.match(r"\b((short)?(\_)?ton(s)?)\b", str(key), re.IGNORECASE):
            return Number(2_000, dims = 'pound')
        # grav
        elif re.match(r"\b(g(0|n|ravity))\b", str(key), re.IGNORECASE):
            return Number(9.80665, dims = 'metre/second^2')
        # time
        elif re.match(r"\b(s(ec(ond)?|egundo)(s)?)\b", str(key), re.IGNORECASE):
            return Number(1, dims = 'second')
        elif re.match(r"\b(min(ut(e|o))?(s)?)\b", str(key), re.IGNORECASE):
            return Number(1, dims = 'minute')
        elif re.match(r"\b(h((ou|o)?r(a)?)(s)?)\b", str(key), re.IGNORECASE):
            return Number(1, dims = 'hour')
        # length
        elif re.match(r"\b(m(et(re|er|ro))?(s)?)\b", str(key), re.IGNORECASE):
            return Number(1, dims = 'metre')
        elif re.match(r"\b(k(ilo)?m(et(re|er|ro))?(s)?)\b", str(key), re.IGNORECASE):
            return Number(1, dims = 'kilometre')
        elif re.match(r"\b(d(eci)?m(et(re|er|ro))?(s)?)\b", str(key), re.IGNORECASE):
            return Number(1, dims = 'decimetre')
        elif re.match(r"\b(c(enti)?m(et(re|er|ro))?(s)?)\b", str(key), re.IGNORECASE):
            return Number(1, dims = 'centimetre')
        elif re.match(r"\b(m(il(l)?i)?m(et(re|er|ro))?(s)?)\b", str(key), re.IGNORECASE):
            return Number(1, dims = 'millimetre')
        elif re.match(r"\b(inch(es)?|pulgada(s)?)\b", str(key), re.IGNORECASE):
            return Number(1, dims = 'inch')
        elif re.match(r"\b(f(o|ee)?t|pie(s)?)\b", str(key), re.IGNORECASE):
            return Number(1, dims = 'foot')
        elif re.match(r"\b(y(ar)?d(a)?(s)?)\b", str(key), re.IGNORECASE):
            return Number(1, dims = 'yard')
        elif re.match(r"\b(mi(le|lla)?(s)?)\b", str(key), re.IGNORECASE):
            return Number(1, dims = 'mile')
        elif re.match(r"\b(n(autical(\_)?)?(mi(le|lla)?(s)?|q))\b", str(key), re.IGNORECASE):
            return Number(6080, dims = 'foot')
        # speed
        elif re.match(r"\b(mi(le|lla)?(s)?p(er)?h((ou|o)?r(a)?)(s)?)\b", str(key), re.IGNORECASE):
            #elif key.lower() in ['mph']:
            return Number(1, dims = 'mile/hour')
        elif re.match(r"\b(knot(s)?)\b", str(key), re.IGNORECASE):
            #elif key.lower() in ['knot']:
            return Number(6080, dims = 'foot/hour')
        # temperature
        elif re.match(r"\b(k(elvin)?(s)?)\b", str(key), re.IGNORECASE):
            #elif key.lower() in ['k', 'kelvin']:
            return Number(1, dims = 'kelvin')
        elif re.match(r"\b((o)?c(elsius)?((\_)?deg(rees)?)?)\b", str(key), re.IGNORECASE):
            #elif key.lower() in ['c','oc', 'c_degrees', 'cdeg', 'celsius']:
            return Number(1, dims = 'celsius')
        elif re.match(r"\b((o)?f(ahrenheit(s)?)?((\_)?deg(rees)?)?)\b", str(key), re.IGNORECASE):
            #elif key.lower() in ['of', 'f_degrees', 'fdeg', 'fahrenheit']:
            return Number(1, dims = 'fahrenheit')
        elif re.match(r"\b(r(ankine(s)?)?((\_)?deg(rees)?)?)\b", str(key), re.IGNORECASE):
            #elif key.lower() in ['or', 'r_degrees', 'rdeg', 'rankine']:
            return Number(1, dims = 'rankine')
        # force
        elif re.match(r"\b(n(ewton)?(s)?)\b", str(key), re.IGNORECASE):
            #elif key.lower() in ['n', 'newton']:
            return Number(1, dims = 'newton')
        elif re.match(r"\b(k(ilo)?n(ewton)?(s)?)\b", str(key), re.IGNORECASE):
            #elif key.lower() in ['kn', 'kilonewton']:
            return Number(1, dims = 'kilonewton')
        elif re.match(r"\b(m(ega)?n(ewton)?(s)?)\b", str(key), re.IGNORECASE):
            #elif key.lower() in ['mn', 'meganewton']:
            return Number(1, dims = 'meganewton')
        elif re.match(r"\b(g(iga)?n(ewton)?(s)?)\b", str(key), re.IGNORECASE):
            #elif key.lower() in ['gn', 'giganewton']:
            return Number(1, dims = 'giganewton')
        elif re.match(r"\b((lbf|pound)(s)?(\_)?f(orce)?)\b", str(key), re.IGNORECASE):
            #elif key.lower() in ['lbf', 'pound_force']:
            return Number(1, dims = 'lbf')
        elif re.match(r"\b(k(ilo)?(ip|lbf|pound)(s)?(\_)?(f(orce)?)?)\b", str(key), re.IGNORECASE):
            #elif key.lower() in ['kip', 'kips', 'kipf', 'klbf']:
            return  Number(1, dims = 'kilolbf')
        elif re.match(r"\b(p(oun)?d(a)?l(s)?)\b", str(key), re.IGNORECASE):
            #elif key.lower() in ['pdl', 'poundal']:
            return Number(1, dims = 'pound * foot / second^2')
        # pressure
        elif re.match(r"\b(pa(scal(s)?)?)\b", str(key), re.IGNORECASE):
            #elif key.lower() in ['pa', 'pascal']:
            return Number(1, dims = 'pascal')
        elif re.match(r"\b(k(ilo)?pa(scal(s)?)?)\b", str(key), re.IGNORECASE):
            #elif key.lower() in ['kpa', 'kilopascal']:
            return Number(1, dims = 'kilopascal')
        elif re.match(r"\b(m(ega)?pa(scal(s)?)?|bar)\b", str(key), re.IGNORECASE):
            #elif key.lower() in ['mpa', 'bar', 'megapascal']:
            return  Number(1, dims = 'megapascal')
        elif re.match(r"\b(g(iga)?pa(scal(s)?)?)\b", str(key), re.IGNORECASE):
            #elif key.lower() in ['gpa', 'gigapascal']:
            return Number(1, dims = 'gigapascal')
        elif re.match(r"\b(mb(ar)|hectopascal(s)?)\b", str(key), re.IGNORECASE):
            #elif key.lower() in ['mbar', 'mb', 'hectopascal']:
            return Number(1, dims = 'hectopascal')
        elif re.match(r"\b(psf)\b", str(key), re.IGNORECASE):
            #elif key.lower() in ['psf']:
            return Number(1, dims = 'lbf / foot^2')
        elif re.match(r"\b(psi)\b", str(key), re.IGNORECASE):
            #elif key.lower() in ['psi']:
            return Number(1, dims = 'psi')
        elif re.match(r"\b(ksi)\b", str(key), re.IGNORECASE):
            #elif key.lower() in ['ksi']:
            return Number(1, dims = 'kilopsi')
        # energy 
        elif re.match(r"\b(w(att)?)\b", str(key), re.IGNORECASE):
            #elif key.lower() in ['w', 'watt']:
            return Number(1, dims = 'watt')
        elif re.match(r"\b(j(oule)?)\b", str(key), re.IGNORECASE):
            #elif key.lower() in ['j', 'joule']:
            return Number(1, dims = 'joule')
        elif re.match(r"\b(hp)\b", str(key), re.IGNORECASE):
            #elif key.lower() in ['hp']:
            return Number(33000, dims = 'foot * lbf / minute')        
        # angle
        elif re.match(r"\b(rad(ian)?(s)?)\b", str(key), re.IGNORECASE):
            #elif key.lower() in ['rad', 'radian', 'radians']:
            return Number(1, dims = 'radian')
        elif re.match(r"\b(deg(ree)?(s)?)\b", str(key), re.IGNORECASE):
            #elif key.lower() in ['deg', 'degree', 'degrees']:
            return Number(1, dims = 'degree')
        elif re.match(r"\b(arcmin)\b", str(key), re.IGNORECASE):
            #elif key.lower() in ['arcmin']:
            return Number(1/60., dims='degree')
        elif re.match(r"\b(arcsec)\b", str(key), re.IGNORECASE):
            #elif key.lower() in ['arcsec']:
            return Number(1/120., dims='degree')
        else:
            raise Exception(" unit item {:} not recognized".format(key))
    #



    
