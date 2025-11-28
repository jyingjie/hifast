

__all__ = ['get_vcorr', 'freq2vel', 'vel2freq', 'frame_correct_freq', 'spectrum_resample', 'correct_spec']


import numpy as np
import scipy.interpolate as interp

import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.coordinates import ICRS, LSR, LSRK, LSRD

from . import conf


def get_vcorr(ra, dec, mjd, frame='BARYCENT'):
    """

    Parameters
    ----------
    ra, dec: scalar or ndarray
             unit in deg
    mjd: scalar or ndarray

    frame: str
        'BARYCENT', 'HELIOCEN', 'LSR', 'LSRK', 'LSRD'

    Returns
    -------
    vcorr: scalar or ndarray
           unit in km/s
    """

    if frame == 'BARYCENT' or 'LSR' in frame:
        kind = 'barycentric'
    elif frame == 'HELIOCEN':
        kind = 'heliocentric'
    else:
        raise(ValueError(f'frame {frame} not support'))
    # telescope location from conf
    geog = EarthLocation.from_geodetic(lat=conf.lat*u.rad, lon=conf.long*u.rad, height=conf.height*u.m)
    c = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame="icrs")
    vcorr = c.radial_velocity_correction(kind=kind, obstime=Time(mjd, format='mjd'), location=geog)
    if frame in ['BARYCENT', 'HELIOCEN']:
        return vcorr.to(u.km/u.s).value

    # to Local standard of rest
    frame_dict = {'LSR': LSR,
              'LSRK': LSRK,
              'LSRD': LSRD,
             }

    if not np.isscalar(ra):
        pm_ra_cosdec = np.full_like(ra, 0)*u.mas/u.yr # must be zero
        pm_dec = np.full_like(ra, 0)*u.mas/u.yr # must be zero
        distance = np.full_like(ra, 1)*u.pc # arbitrary value
    else:
        pm_ra_cosdec = 0*u.mas/u.yr
        pm_dec = 0*u.mas/u.yr
        distance = 1*u.pc
    vel = ICRS(ra=ra*u.deg, dec=dec*u.deg, pm_ra_cosdec=pm_ra_cosdec, pm_dec=pm_dec, radial_velocity=vcorr, distance=distance)
    vcorr = vel.transform_to(frame_dict[frame]()).radial_velocity

    return vcorr.to(u.km/u.s).value



def freq2vel(freq, vtype='radio'):
    """
    freq to velocity

    Parameters
    ----------
    freq: array
         shape: (n,); unit in Mhz
    vtype: str
          'radio' or 'optical', default is 'radio'
    """
    restfreq = conf.restfreq
    if vtype == 'radio':
        vel = (restfreq - freq)/restfreq*conf.vlight
    elif vtype == 'optical':
        vel = (restfreq - freq)/freq*conf.vlight
    return vel

def vel2freq(vel, vtype='radio'):
    """
    velocity to freq

    Parameters
    ----------
    vel: array
         shape: (n,); unit in km/s
    vtype: str
          'radio' or 'optical', default is 'radio'
    """
    restfreq = conf.restfreq
    if vtype == 'radio':
        freq = restfreq * (1 - vel / conf.vlight)
    elif vtype == 'optical':
        freq = restfreq / (vel / conf.vlight + 1)
    return freq

def _freq_vcorr(freq, vcorr):
    """
    doppler correct ``freq`` (Mhz) with ``vcorr`` (km/s)
    """
    return freq * np.sqrt((conf.vlight - vcorr)/(conf.vlight + vcorr))

def frame_correct_freq(freq, ra, dec, mjd, frame='BARYCENT'):
    """
    dopper correct the ``freq`` of the spectra at (ra, dec) observed at time ``mjd``

    freq : array like (m,)
          unit in Mhz
    mjd, ra, dec: scalar or array (n,)
    frame: str; 'BARYCENT', 'HELIOCEN', 'LSR', 'LSRK', 'LSRD'

    return: shape (m,) if mjd, ra, dec are scalar; shape (n, m) if array
    """
    isscalar = False
    if np.isscalar(mjd):
        mjd = np.array([mjd])
        isscalar = True
    if np.isscalar(ra):
        ra = np.array([ra])
        isscalar = True
    if np.isscalar(dec):
        dec = np.array([dec])
        isscalar = True

    vcorr = get_vcorr(ra, dec, mjd, frame)
    if not isscalar:
        vcorr = vcorr[:, None]
    return _freq_vcorr(freq, vcorr)



def spectrum_resample(freq, flux, freq_new, method='interp', interp_kind='linear'):
    """
    resample spectrum with ``flux`` at ``freq`` to `freq_new`
    """
    if method == 'interp':
        flux_new = interp.interp1d(freq, flux, kind=interp_kind, bounds_error=False, fill_value=np.nan)(freq_new)
        return flux_new

def correct_spec(Ta, freq, ra, dec, mjd, frame='LSRK', method='interp', interp_kind='linear'):
    """
    doppler correct spectra with ``Ta`` at ``freq`` and resample the spectra at same frequency sample points

    Parameters
    ----------
    Ta: array; shape: (m,n,2) or (m,n)
    freq: array; shape: (n,) in ascending order
    ra, dec: array; shape: (m,)
    frame: str
          rest frame, 'BARYCENT', 'HELIOCEN', 'LSR', 'LSRK', 'LSRD'
    """
    if Ta.ndim == 2:
        Ta = Ta[:,:,None]
        ndim_ori = 2
    elif Ta.ndim == 3:
        ndim_ori = 3
    else:
        raise(ValueError('input ndim 2 or 3'))

    # velocity of observer (telescope) with respect to frame
    vcorr = get_vcorr(ra, dec, mjd, frame)

    # new freq sample points
    f_d = freq[1] - freq[0]
    fcorr = _freq_vcorr(freq, np.median(vcorr))
    tmp = np.hstack([np.arange(freq[0]-f_d, fcorr.min(),-f_d)[::-1], freq, np.arange(freq[-1]+f_d, fcorr.max(),f_d)])
    freq_new = tmp[(tmp> fcorr.min()) & (tmp< fcorr.max())]

    Ta_new = np.zeros((Ta.shape[0], len(freq_new), Ta.shape[2]))
    for i in range(Ta.shape[0]):
        for j in range(Ta.shape[2]):
            fcorr_i = _freq_vcorr(freq, vcorr[i])
            Ta_new[i,:,j] = spectrum_resample(fcorr_i, Ta[i,:,j], freq_new, method=method, interp_kind=interp_kind)
    if ndim_ori == 2:
        Ta_new = Ta_new[..., 0]
    return Ta_new, freq_new
