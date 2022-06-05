#!/usr/bin/env python
# coding: utf-8

import numpy as np
from glob import glob
import os
import sys
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
import h5py
from .core import grid

def _adjust_header(header, ra_range, dec_range):
    """
    adjust header to fit ra dec range, inplace change
    """
    # points used to test
    radec_t = map(lambda x:[x[0], x[-1], (x[0]+x[-1])/2], [ra_range, dec_range])
    world_t = np.vstack(tuple(map(np.ravel, np.meshgrid(*radec_t) + [np.full(9,0)]))).T
    # adjust
    n1 = 2
    n2 = 2
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
    return header

def gen_header(ra_range, dec_range, x_delta, y_delta, z, proj, vel_type, frame, histories=None):
    
    x= np.arange(ra_range[0],ra_range[1]+x_delta,x_delta)[::-1] # reverse ra
    y= np.arange(dec_range[0],dec_range[1]+y_delta,y_delta)
    # Create a new WCS object.
    w = WCS(naxis=3)
    # center pixel
    w.wcs.crpix = list(map(lambda x:len(x)/2, [x,y,z]))
    # coordinate and z value of that pixel.
    w.wcs.crval = list(map(lambda x:(x[-1]+x[0])/2, [x,y,z]))
    # the pixel scale in (ra,dec, z)
    w.wcs.cdelt = list(map(lambda x:(x[-1]-x[0])/(len(x)-1), [x,y,z]))
    # projection
    w.wcs.ctype = [f"RA---{proj}", f"DEC--{proj}", vel_type]
    w.wcs.specsys= frame # or HELIOCENT
    
    header = w.to_header()
    
    header["NAXIS"] = 3
    header["NAXIS1"] = len(x)
    header["NAXIS2"] = len(y)
    header["NAXIS3"] = len(z)
    header = _adjust_header(header, ra_range, dec_range)
    header["CUNIT1"] = 'deg'
    header["CUNIT2"] = 'deg'
    header["CUNIT3"] = 'km/s'
    
    # additional
    header["EQUINOX"] = 2000.0                                                  
    header["LINE"]    = 'HI'
    header["RESTFRQ"]  =   1.420405751E+9
    header["BMAJ"] = 2.9/60
    header["BMIN"] = 2.9/60
    header["BPA"] = 0.0
    
    #add history
    if histories is not None:
        for his in histories:
            header['HISTORY']= his
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

def _stack_spec(vel, Ta, vrange=None):
    """
    vlrs, Ta: list; vel in descending order
    """
    if vrange is None:
        ranges= np.array([np.min([i[0]for i in vel]), np.max([i[-1] for i in vel])])
        delta = (vel[0][0] - vel[0][1])/2.
        vmax,vmin= ranges.max() + delta, ranges.min()-delta
    else:
        vmax,vmin=vrange[1],vrange[0]
    _vel = []
    _Ta = []
    for i,j in zip(vel,Ta):
        is_use= (i>=vmin)&(i<=vmax)
        _vel+= [i[:][is_use]]
        _Ta+= [j[:][:,is_use]]

    vel_len= np.array([len(i) for i in _vel])
    len_use= vel_len.min()

    vel= np.vstack([i[:len_use] for i in _vel])
    Ta= np.vstack([i[:,:len_use] for i in _Ta])
    return vel,Ta

def _preprocess(vals, threshold = None):
    """
    process vals in place
    Parameters
    -----------
    vals: Ta or flux
    threshold: vals less than threshold will be masked as nan
    
    """
    if threshold is not None:
        vals[np.abs(vals) > threshold] = np.nan
    return vals

def remove_nan(vel, Ta):
    is_use= np.sum(np.isnan(Ta),axis=0)==0
    return vel[:,is_use],Ta[:,is_use]

def get_ra_range(ra):
    ra_s= np.sort(ra)
    diff_s= np.diff(ra_s)
    ind_max= np.argmax(diff_s)
    if diff_s[ind_max] >= ra_s[0]+360-ra_s[-1]:
        ra= np.copy(ra)
        is_c= ra>=ra_s[ind_max+1]
        ra[is_c]= ra[is_c]-360
    return [np.nanmin(ra), np.nanmax(ra)]
    
if __name__ == '__main__':
    import os
    import argparse
    from hifast.utils.io import PolarMjdChan_to_MjdChanPolar
    from hifast.utils.io import formatter_class
    parser = argparse.ArgumentParser(allow_abbrev=False, formatter_class=formatter_class)
    parser.add_argument('fname',nargs='+',
                        help='file name')
    parser.add_argument('--ra_range', type=float,nargs=2,
                       help='ra range; unit: deg')
    parser.add_argument('--dec_range', type=float,nargs=2,
                       help='dec_range; unit: deg')
    parser.add_argument('--range3', type=float,nargs=2,
                       help='freq or vel range')
    parser.add_argument('--type3', choices=['vel', 'freq'], default='vel',
                       help='third axis in cube,  freq or vel')
    parser.add_argument('--bwidth', type=float, default=[60.],nargs='+',
                       help='unit: arc second')
    parser.add_argument('--outname', required=True,
                       help='output file name, full path')
    parser.add_argument('-f', action='store_true', dest='force',
                        help='if set, overwriting file if output file exists')
    parser.add_argument('-k', '--key', choices=['flux', 'Ta'],
                       help='properties to make fits cube. flux or Ta')
    parser.add_argument('-p','--proj', default='AIT', #choices=['AIT', 'SIN', 'TAN'],
                       help='fits wcs projection')
    parser.add_argument('--r_cut', type=float, default=90,
                       help='spectra inside r_cut from the grid point will be considered; unit: arc second')
    parser.add_argument('-m','--method', default='bessel_gaussian', choices=['reweight', 'mean', 'median', 'gaussian', 'bessel_gaussian', 'sinc_gaussian'],
                       help='method to process the spec in r_cut')
    parser.add_argument('--frac_finite_min', '--frac-finite-min', type=float, default=1,
                       help='For a grid point having ``n``` spectra in ``r_cut``, if the number of ``finite value`` in a channel is zero or smaller than ``frac_finite_min*n``, the output value in the channel will be set as ``nan``')
    parser.add_argument('-t','--threshold', type=float,
                       help='vals less than threshold will be masked as nan')
    parser.add_argument('--w_on_t', action='store_true',
                       help='weight on sample time')

    args = parser.parse_args()
    fname= args.fname
    ra_range= args.ra_range
    dec_range= args.dec_range
    bwidth= args.bwidth
    outname= args.outname
    key= args.key
    type3 = args.type3
    proj= args.proj
    r_cut= args.r_cut
    r_cut=r_cut/3600 # convert to deg
    method= args.method
    threshold = args.threshold
    range3 = args.range3
    
    if os.path.exists(outname):
        if args.force:
            print(f"will overwrite the existing output file {outname}")
        else:
            print(f"File exists {outname}")
            print("exit... Use ' -f ' to overwrite it or change outname.")
            sys.exit(0)
            
    #make sure bwidth not to samll
    if (np.array(bwidth)<1).any():
        raise(ValueError('input --bwidth is too small, its unit is arc second'))
    
    
    files=[]
    for thefname in fname:
        files+= glob(thefname,recursive=True)
    files= list(set(files))
    files.sort()
    print('processing files:')
    [print(i) for i in files]
    
    with h5py.File(files[0],'r') as f:
        try:
            vel_type= f['Header'].attrs['vel_type']
            frame= f['Header'].attrs['frame']
        except:
            print("No header found, using arbitrary value for velcity type and frame")
            vel_type= 'VRAD'
            frame='LSRK'
        if key is None:
            if 'flux' in f['S'].keys():
                key= 'flux'
            elif 'Ta' in f['S'].keys():
                key= 'Ta'
            else:
                raise(ValueError('no flux or temperature information'))
            
    ra=[]
    dec=[]
    Ta=[]
    vel=[]
    fs=[]
    
    if args.w_on_t:
        t_sample = []
    
    for file in files:
        f = h5py.File(file,'r')
        S = f['S']
        ra += [S['ra'],]
        dec += [S['dec'],]
        Ta_ = PolarMjdChan_to_MjdChanPolar(S[key][:])
        # merge polar
        if Ta_.ndim == 3:
            Ta_ = np.mean(Ta_, axis=2, dtype='float64')
        vel_ = S[type3]
        if range3 is not None:
            is_ = (vel_[:] >= range3[0]) & (vel_[:] <= range3[1])
            vel_ = vel_[:][is_]
            Ta_ = Ta_[:][:, is_]
        Ta += [Ta_, ]
        vel += [vel_,]
        fs += [f,]
        
        if args.w_on_t:
            t_sample += [S['mjd'][1] - S['mjd'][0]]
    if args.w_on_t:
        t_sample = np.array(t_sample)
        t_sample /= np.max(t_sample)
        wi = np.hstack([np.full(len(ra[i]), t_sample[i]) for i in range(len(t_sample))])
    else:
        wi = None
        
    ra= np.hstack(ra)
    dec= np.hstack(dec)
    #     Ta=np.vstack(Ta)
    #     vel=np.vstack(vel)
    vel, Ta= _stack_spec(vel, Ta) #vel are in descending order.
    [f.close() for f in fs]
    #vel, Ta= remove_nan(vel, Ta)
    Ta= _preprocess(Ta, threshold)
    # use some vel sample for each specta
    std_vel= np.std(vel,axis=0,dtype=np.float64) # single precision can be inaccurate
    if np.nanmax(std_vel)>0.2:
        raise(ValueError('vel dispersion (std) at same Ta order is two large'))
    else:
        print(f'vel dispersion (std) at same Ta order is between {np.nanmin(std_vel)} and {np.nanmax(std_vel)}.')
    vel_refine= np.mean(vel,axis=0,dtype=np.float64) # single precision can be inaccurate
    vel_refine= np.append(vel_refine- (vel_refine[1]-vel_refine[0])/2, vel_refine[-1]+(vel_refine[1]-vel_refine[0])/2)

    if ra_range is None:
        ra_range= get_ra_range(ra)
    if dec_range is None:
        dec_range= [np.nanmin(dec), np.nanmax(dec)]
    
    bwidth= [i/60/60 for i in bwidth]
    if len(bwidth)==1:
        bwidth=bwidth*2
    print('ra range:',ra_range,' deg')
    print('dec range:',dec_range,' deg')
    print('bwidth:',bwidth,' deg')
    
    #record history
    histories = []
    try:
        from ._version import get_versions
        histories = ['version: ' + get_versions()['version']]
        histories += [f"{k}: {v}" if k !='fname' else "args" for k,v in args.__dict__.items()]
    except:
        pass
    histories += files
    
    header= gen_header(ra_range,dec_range,*bwidth,vel[0],proj, vel_type, frame, histories=histories)
    ra_grid, dec_grid = gen_grid_radec(header)
    #print('ra range in generated cube fits file', np.min(ra_grid), np.max(ra_grid))
    #print('dec range in generated cube fits file', np.min(dec_grid), np.max(dec_grid))
    out, nums= grid.gridding(ra, dec, Ta, ra_grid, dec_grid, wi=wi, r=r_cut, method=method, 
                            frac_finite_min=args.frac_finite_min) 
    hdu = fits.PrimaryHDU(out.astype('float32'), header=header)
    print(f'Saving to {outname}.')
    overwrite = True if args.force else False
    hdu.writeto(outname,overwrite=overwrite)
    # save the spec count in each grid
    outname_c = '.'.join(outname.split('.')[:-1]) + '-count.fits'
    hdu = fits.PrimaryHDU(nums, header=header_3to2(header))
    hdu.writeto(outname_c, overwrite=True)
