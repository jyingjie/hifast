import ctypes as ct
from ctypes import c_double,c_int
from glob import glob
import os

_sofa_lib= glob(os.path.dirname(__file__)+'/_sofa_c.*.so')
if len(_sofa_lib) == 0:
    raise ImportError('can not find sofa lib file.')
_sofa= ct.CDLL(_sofa_lib[0])

def Atoc13(aob, zob, utc1, utc2, dUT1, elong, phi, hm, phpa, temperature, humidity, wl):
    rc = c_double()
    dc = c_double()

    tmp= _sofa.iauAtoc13(b"A", c_double(aob), c_double(zob), c_double(utc1), c_double(utc2), c_double(dUT1), c_double(elong), 
                c_double(phi), c_double(hm), c_double(0.), c_double(0.), c_double(phpa), c_double(temperature), c_double(humidity),
                c_double(wl), ct.byref(rc), ct.byref(dc))
    
    if tmp != 0:
        return [None, None]
    return [rc.value, dc.value]


def Dtf2d(year, month, day, hour, minute, second):
    utc1 = c_double()
    utc2 = c_double()

    tmp= _sofa.iauDtf2d(b"UTC", c_int(year), c_int(month), c_int(day), c_int(hour),
                            c_int(minute), c_double(second),
                   ct.byref(utc1), ct.byref(utc2))
    if tmp != 0:
        return [None, None]
    return utc1.value,utc2.value


_sofa.iauA2tf.argtypes = [
    ct.c_int, # tt1
    ct.c_double, # tt2
    ct.POINTER(ct.c_char), #ut11
    ct.POINTER(ct.c_int)
]
def A2tf(ndp=2,rc=None,plusChar = b'+'):
    """
    rc: angle in radians
    Return: hours, minutes, seconds, fraction
    """
    plusChar = b'+'
    #int hour, minute, second, milisecond;

    ihmsf= (ct.c_int*4)()

    tmp= _sofa.iauA2tf(2, rc, plusChar, ihmsf)     
    
    if tmp != 0:
        return [None]*4
    
    return list(ihmsf)