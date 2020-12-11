#!/usr/bin/env python
# coding: utf-8

import scipy.interpolate as interp
import numpy as np
import os
from .ugdopplerfast import ugdopplerfast

def freq2vel(freq, vtype='radio'):
    """
    rest freq to velocity
    
    Parameters
      ----------
    freq: array; shape: (n,)
    vtype: str
          radio or optical, default is radio
    """
    #velocity
    restfreq = 1420.405751
    if vtype == 'radio':
        vel = 299792.458*(restfreq - freq)/ restfreq
    elif vtype == 'optical':
        vel = 299792.458*(restfreq - freq)/ freq
    return vel

def frame_correct(Ta, freq, mjd, ra, dec, frame='LSRK', interp_kind='linear'):
    """
    rest frame
    
    Parameters
      ----------
    Ta: array; shape: (m,n,2) or (m,n)
    freq: array; shape: (n,)
    ra, dec: array; shape: (m,)
    frame: str
          rest frame, HELIOCENT or LSRK
    """
    if Ta.ndim ==2:
        Ta = Ta[:,:,None]
        
    # obs vel relate to rest frame
    jd = mjd + 2400000.5
    if frame == 'HELIOCENT':
        vobs =  ugdopplerfast(ra, dec, jd, frame=frame) # velocity of observer (telescope) with respect to frame
    elif frame == 'LSRK':
        vobs =  - ugdopplerfast(ra, dec, jd, frame=frame) # velocity of observer (telescope) with respect to frame
    
    # new freq sample points
    f_d = freq[1] - freq[0]
    wvlo = freq_correct(freq, np.median(vobs))
    tmp = np.hstack([np.arange(freq[0]-f_d,wvlo.min(),-f_d)[::-1], freq, np.arange(freq[-1]+f_d,wvlo.max(),f_d)])
    wvlo = tmp[(tmp> wvlo.min()) & (tmp< wvlo.max())]
    
    
    Ta_new = np.zeros((Ta.shape[0],len(wvlo),Ta.shape[2]))
    for i in range(Ta.shape[0]):
        for j in range(Ta.shape[2]):
            Ta_new[i,:,j], _ = doppler_correct(freq, Ta[i,:,j], vobs[i], wvlo, method='interp', interp_kind=interp_kind)
    freq_new = wvlo
    if Ta_new.shape[2] ==1:
        Ta_new = Ta_new[...,0]
    return Ta_new, freq_new

def doppler_correct(wvl, flux, v, wvlo=None, method='interp', **kwargs):
    if method=='interp':
        return doppler_correct_inter(wvl, flux, v, wvlo=wvlo, **kwargs)

# def freq_correct(wvl,v):
#     return wvl * (1.0 - v/299792.458)
def freq_correct(wvl,v):
    return wvl * np.sqrt((299792.458-v)/(299792.458+v))

def doppler_correct_inter(wvl, flux, v, wvlo=None, edgeHandling=None, fillValue=None, interp_kind='linear'):
    """
      Doppler Correct a given spectrum.
      Modified from PyAstronomy.pyasl.dopplerShift, https://pyastronomy.readthedocs.io/en/latest/pyaslDoc/aslDoc/dopplerShift.html
      
      .. warning:: Shifting a spectrum using linear
                  interpolation has an effect on the
                  noise of the spectrum. No treatment
                  of such effects is implemented in this
                  function.

      Parameters
      ----------
      wvl : array
          Input wavelengths.
      flux : array
          Input flux.
      v : float
          Doppler shift in km/s
      wvlo : array
          Output wavelengths in A.
      edgeHandling : string, {"fillValue", "firstlast"}, optional
          The method used to handle the edges of the
          output spectrum.
      fillValue : float, optional
          If the "fillValue" is specified as edge handling method,
          the value used to fill the edges of the output spectrum.

      Returns
      -------
      nflux : array
          The shifted flux array at the *old* input locations if wvlo is None, else at locations of wvlo.
      wlprime : array
          The shifted wavelength axis.
    """
    # Shifted wavelength axis
    wlprime = freq_correct(wvl,v)
    fv = np.nan
    if edgeHandling == "fillValue":
        if fillValue is None:
            raise(ValueError("Fill value not specified, If you request 'fillValue' as edge handling method, you need to specify the 'fillValue' keyword."))
        fv = fillValue

    f = interp.interp1d(wlprime, flux, kind=interp_kind, bounds_error=False, fill_value=fv)
    if wvlo is None:
        nflux = f(wvl)
    else:
        nflux = f(wvlo)
    if edgeHandling == "firstlast":
        firsts = []
        # Search for first non-NaN value save indices of
        # leading NaN values
        for i in range(len(nflux)):
            if np.isnan(nflux[i]):
                firsts.append(i)
            else:
                firstval = nflux[i]
                break
        # Do the same for trailing NaNs
        lasts = []
        for i in range(len(nflux)-1, 0, -1):
            if np.isnan(nflux[i]):
                lasts.append(i)
            else:
                lastval = nflux[i]
                break
        # Use first and last non-NaN value to
        # fill the nflux array
        nflux[firsts] = firstval
        nflux[lasts] = lastval
    return nflux, wlprime


