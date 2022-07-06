#!/usr/bin/env python
# coding: utf-8

import numpy as np
from glob import glob
import os
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
import h5py
from hifast.core import grid
from hifast.core import conf
from hifast.core.corr_vel import freq2vel, vel2freq
from hifast.core.wcs import gen_header, gen_grid_radec, header_3to2

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
    print(np.diff(vel))
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
    from hifast.utils.io import get_nB
    parser = argparse.ArgumentParser(allow_abbrev=False, formatter_class=formatter_class)
    parser.add_argument('fname', nargs='+',
                        help='file name')
    parser.add_argument('--ra_range', type=float,nargs=2,
                       help='ra range; unit: deg')
    parser.add_argument('--dec_range', type=float,nargs=2,
                       help='dec_range; unit: deg')
    parser.add_argument('--range3', type=float,nargs=2,
                       help='freq or vel range')
    parser.add_argument('--type3', choices=['vopt', 'vrad', 'freq'], default='vrad',
                       help='third axis in cube,  ``vopt, ``vrad`` or ``freq``')
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
    parser.add_argument('-m','--method', default='gaussian', choices=['mean','median', 'gaussian', 'bessel_gaussian', 'sinc_gaussian'],
                       help='method to process the spec in r_cut')
    parser.add_argument('--frac_finite_min', '--frac-finite-min', type=float, default=1, help='For a grid point having ``n``` spectra in ``r_cut``, if the number of ``finite value`` in a channel is zero or smaller than ``frac_finite_min*n``, the output value in the channel will be set as ``nan``')
    parser.add_argument('-t','--threshold', type=float,
                       help='vals less than threshold will be masked as nan')
    parser.add_argument('--w_on_t', action='store_true',
                       help='weight on sample time')
    parser.add_argument('--wcs_from',
                       help='use the wcs parameters from input fits')
    

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
            #vel_type= f['Header'].attrs['vel_type']
            frame= f['Header'].attrs['frame']
        except:
            print("No header found, using arbitrary value for velcity type and frame")
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
    nBs = []
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
        _freq = S['freq'][()]
        if type3.upper() == 'FREQ':
            vel_ = _freq
        elif type3.upper() == 'VOPT':
            vel_ = freq2vel(_freq, vtype='optical')
        elif type3.upper() == 'VRAD':
            vel_ = freq2vel(_freq, vtype='radio')
        if range3 is not None:
            is_ = (vel_[:] >= range3[0]) & (vel_[:] <= range3[1])
            vel_ = vel_[:][is_]
            Ta_ = Ta_[:][:, is_]
        Ta += [Ta_, ]
        vel += [vel_,]
        nBs += [np.full(len(S['ra'][:]), get_nB(file))]
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
    nBs = np.hstack(nBs)
    #     Ta=np.vstack(Ta)
    #     vel=np.vstack(vel)
    vel, Ta= _stack_spec(vel, Ta) #vel are in descending order.
    [f.close() for f in fs]
    #vel, Ta= remove_nan(vel, Ta)
    Ta= _preprocess(Ta, threshold)
    # use some vel sample for each specta
    std_vel= np.std(vel,axis=0,dtype=np.float64) # single precision can be inaccurate
    print(f'vel dispersion (std) at same Ta order is between {np.nanmin(std_vel)} and {np.nanmax(std_vel)}.')
    if np.nanmax(std_vel)>0.2:
        raise(ValueError('vel dispersion (std) at same Ta order is two large'))
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
    
    header = gen_header(ra_range, dec_range, *bwidth,vel[0], proj, type3, frame, histories=histories)
    if key == 'flux':
        header["BUNIT"] = 'Jy/beam'
    elif key == 'Ta':
        header["BUNIT"] = 'K'
    
    if args.wcs_from is not None:
        print(f'use the wcs sky coordinates parameters from {args.wcs_from}')
        fa = fits.open(args.wcs_from)
        header2 = fa[0].header
        for key in header.keys():
            if key[-1:] in ['1', '2'] or key == 'LONPOLE' or key == 'LATPOLE':
                print(f'replacing {key}')
                try:
                    header[key] = header2[key]
                except:
                    print(f'replace {key} fail')
    
    ra_grid, dec_grid = gen_grid_radec(header)
    #print('ra range in generated cube fits file', np.min(ra_grid), np.max(ra_grid))
    #print('dec range in generated cube fits file', np.min(dec_grid), np.max(dec_grid))
    out, nums= grid.gridding(ra, dec, Ta, ra_grid, dec_grid, wi=wi, r=r_cut, method=method, frac_finite_min=args.frac_finite_min) 
    hdu = fits.PrimaryHDU(out.astype('float32'), header=header)
    print(f'Saving to {outname}.')
    overwrite = True if args.force else False
    hdu.writeto(outname,overwrite=overwrite)
    # save the spec count in each grid
    outname_c = '.'.join(outname.split('.')[:-1]) + '-count.fits'
    hdu = fits.PrimaryHDU(nums, header=header_3to2(header))
    hdu.writeto(outname_c, overwrite=True)
