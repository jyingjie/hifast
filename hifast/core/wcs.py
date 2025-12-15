

__all__ = ['gen_header', 'header_3to2', 'gen_grid_radec']


import numpy as np
from astropy.wcs import WCS
from astropy import units as u
from . import conf
from .corr_vel import freq2vel, vel2freq


def _adjust_header(header, ra_range, dec_range, n1=3, n2=3):
    """
    adjust header to fit ra dec range, inplace change
    """
    # points used to test
    radec_t = map(lambda x:[x[0], x[-1], (x[0]+x[-1])/2], [ra_range, dec_range])
    world_t = np.vstack(tuple(map(np.ravel, np.meshgrid(*radec_t) + [np.full(9,0)]))).T
    # adjust
    for i in range(n1+n2):
        w = WCS(header)
        pix_t = w.wcs_world2pix(world_t,0) # pixel of the test points
        pix_t_min = np.min(pix_t,axis=0)
        pix_t_max = np.max(pix_t,axis=0)
#         print('init',end=':')
#         print(' right ', [header['NAXIS1'],header['NAXIS2']]- np.max(pix_t,axis=0)[:2])
#         print('    left  ',np.min(pix_t,axis=0)[:2])

        if i < n1:
            naxis_new = np.ceil(pix_t_max - pix_t_min).astype('int')
            header['CRPIX1'] = naxis_new[0]/2
            header['CRPIX2'] = naxis_new[1]/2
            header['NAXIS1'] = naxis_new[0]
            header['NAXIS2'] = naxis_new[1]
            #print('v1',end=':')
        else:
            pix_t_min = np.floor(pix_t_min).astype('int') - 1 # number 1 is for redundancy
            pix_t_max = np.ceil(pix_t_max).astype('int') + 1
            header['CRPIX1'] -= pix_t_min[0]
            header['CRPIX2'] -= pix_t_min[1]
            header['NAXIS1'] += (pix_t_max[0] - header['NAXIS1'] - pix_t_min[0])
            header['NAXIS2'] += (pix_t_max[1] - header['NAXIS2'] - pix_t_min[1])
            #print('v2',end=':')
        #print(header)
    return header


def gen_header(ra_range, dec_range, ra_delta, dec_delta, z, proj='SIN', type3='vrad', frame='',
               histories=None, beam_fwhw=2.9,
               ra_center=None, dec_center=None):

    ra = np.arange(ra_range[0], ra_range[1] + ra_delta, ra_delta)[::-1]
    dec = np.arange(dec_range[0], dec_range[1] + dec_delta, dec_delta)
    # Create a new WCS object.
    w = WCS(naxis=3)
    if type3.upper() == 'FREQ':
        ctype3 = 'FREQ'
    elif type3.upper() == 'VOPT':
        ctype3 = 'VOPT-F2W'
    elif type3.upper() == 'VRAD':
        ctype3 = 'VRAD'
    else:
        raise(ValueError('type3'))

    if proj.lower() == 'rcar':
        w.wcs.ctype = [f"RA---CAR", f"DEC--CAR", ctype3]
        w.wcs.lonpole = 0.
        ra_center = 120.
        dec_center = 0.
    else:
        w.wcs.ctype = [f"RA---{proj}", f"DEC--{proj}", ctype3]
    # center pixel
    w.wcs.crpix = [len(ra)/2,
                   len(dec)/2,
                   1,]
    # coordinate and z value of that pixel.
    w.wcs.crval = [(ra[0] + ra[-1])/2,
                   (dec[0] + dec[-1])/2,
                   z[0]]
    # the pixel scale in (ra,dec, z)
    w.wcs.cdelt = list(map(lambda x:x[1]-x[0], [ra, dec, z]))
    if type3.upper() == 'VOPT':
        freq_r = vel2freq(z, 'optical')
        cdelt3 = -conf.restfreq/(freq_r[0]**2)*(freq_r[1]-freq_r[0])*conf.vlight
        w.wcs.cdelt[2] = cdelt3

    if ra_center is not None:
        w.wcs.crpix[0] = (ra_center - w.wcs.crval[0])/w.wcs.cdelt[0] + w.wcs.crpix[0]
        w.wcs.crval[0] = ra_center
    if dec_center is not None:
        w.wcs.crpix[1] = (dec_center - w.wcs.crval[1])/w.wcs.cdelt[1] + w.wcs.crpix[1]
        w.wcs.crval[1] = dec_center

    w.wcs.cunit = ['deg', 'deg', 'km/s']
    if type3.upper() == 'FREQ':
        w.wcs.cunit[2] = 'MHz'
    # projection
    w.wcs.specsys = frame # or HELIOCENT
    w.wcs.restfrq = conf.restfreq*1e6 # Mhz to Hz # 1.420405752E+9 hz for HI

    header = w.to_header()
    if proj.lower() == 'rcar':
        header["LONPOLE"] = 0.0

    header["NAXIS"] = 3
    header["NAXIS1"] = len(ra)
    header["NAXIS2"] = len(dec)
    header["NAXIS3"] = len(z)
    header = _adjust_header(header, ra_range, dec_range)


    # additional
    header["EQUINOX"] = 2000.0
    header["LINE"] = 'HI'
    header["BMAJ"] = beam_fwhw/60
    header["BMIN"] = beam_fwhw/60
    header["BPA"] = 0.0
    header["BUNIT"] = ''

    #add history
    if histories is not None:
        for his in histories:
            header['HISTORY']= his
    #print(header)
    return header


def header_3to2(header_3d):
    import copy
    header_2d = copy.deepcopy(header_3d)

    keys_rm = ['CRPIX3','CDELT3','CUNIT3','CTYPE3','CRVAL3','NAXIS3',]
    for key in keys_rm:
        header_2d.pop(key)
    header_2d['WCSAXES'] = 2
    header_2d['NAXIS'] = 2

    return header_2d

def gen_grid_radec(header):
    """
    https://github.com/radio-astro-tools/spectral-cube/blob/master/spectral_cube/base_class.py; world
    """
    wcs = WCS(header)
    inds = np.ogrid[[slice(0, s) for s in (1,header['NAXIS2'],header['NAXIS1'])]]
    inds = np.broadcast_arrays(*inds)
    # view=tuple([0 for ii in range(3 - 2)] + [slice(None)] * 2)
    # inds = [i[view] for i in inds[::-1]]
    inds= [i[0,:,:] for i in inds[::-1]]

    shape = inds[0].shape
    inds = np.column_stack([i.ravel() for i in inds])
    world = wcs.all_pix2world(inds,0).T
    world = [w.reshape(shape) for w in world]
    world = [w * u.Unit(wcs.wcs.cunit[i]) for i, w in enumerate(world)]
    world = world[::-1]
    _, dec_grid, ra_grid= world
    return ra_grid, dec_grid
