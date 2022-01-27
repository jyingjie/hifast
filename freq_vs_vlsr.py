#!/usr/bin/env python
# coding: utf-8
import numpy as np
from fast_python.util import freq2vlsr
from astropy.time import Time
from astropy import units as u

obsdate= Time('2020-01-01T20:00:00',format='isot')- 8*u.hour # CST to UTC
ra,dec= 10.68470833, 41.26875 #deg, M31
freq= np.arange(1410,1440,1)

vlsr= freq2vlsr(freq, ra, dec, obsdate.mjd)

print("M31")
print("obs date: ", obsdate.to_datetime().strftime('%Y %m %d'),'UTC')
print('freq    V_LSR(km/s)')
print(np.vstack([freq,vlsr]).T)
