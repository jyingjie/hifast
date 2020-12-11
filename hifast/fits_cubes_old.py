#!/usr/bin/env python
# coding: utf-8

import numpy as np
from glob import glob
from scipy.stats import binned_statistic_2d
from astropy.io import fits
from astropy import wcs
import h5py

def bin_stat(ra,dec,Ta,ranges,bwidth=5, valmin=-np.inf):
    
    x,y,values= ra, dec, Ta.T # Ta.T to fit the input of binned_statistic_2d
    if np.isscalar(bwidth):
        bwidth1=bwidth2=bwidth
    else:
        bwidth1,bwidth2=bwidth
    
    bins=[np.arange(ranges[0][0],ranges[0][1]+bwidth1,bwidth1), np.arange(ranges[1][0],ranges[1][1]+bwidth2,bwidth2)]

    val,xedge,yedge,_= binned_statistic_2d(x,y,values,statistic='mean',bins=bins)
    val[val<valmin]=np.nan
    ra_refine= xedge # (xedge[:-1]+xedge[1:])/2
    dec_refine= yedge #(yedge[:-1]+yedge[1:])/2
    return val,ra_refine,dec_refine


def fits_cubes(x,y,z,vals,fitsname):
    # Create a new WCS object.  The number of axes must be set
    # from the start
    w = wcs.WCS(naxis=3)
    # center pixel of the XY grid.
    w.wcs.crpix = list(map(lambda x:len(x)/2, [x,y,z]))
    # coordinate and z value of that pixel.
    w.wcs.crval = list(map(lambda x:(x[-1]+x[0])/2, [x,y,z]))
    # the pixel scale in (ra,dec, z)
    w.wcs.cdelt = list(map(lambda x:(x[-1]-x[0])/(len(x)-1), [x,y,z]))
    # projection? 
    w.wcs.ctype = ["RA---CAR", "DEC--CAR", "VELO-LSR"]

    # write the HDU object WITH THE HEADER
    header = w.to_header()
    hdu = fits.PrimaryHDU(vals, header=header,)

    hdu.header["EQUINOX"] = 2000.0                                                  
    hdu.header["LINE"]    = 'HI'
    hdu.header["CUNIT3"]  =  'km/s' 
    hdu.header["RESTFRQ"]  =   1.420405751E+9 
    hdu.writeto(fitsname,overwrite=False)

def _stack_spec(vlsr, Ta, vrange=None):
    """
    vlrs, Ta: list; vlsr in descending order
    """
    if vrange is None:
        ranges= np.array([np.min([i[0]for i in vlsr]), np.max([i[-1] for i in vlsr])])
        vmax,vmin= ranges.max(), ranges.min()
    else:
        vmax,vmin=vrange[1],vrange[0]
    _vlsr = []
    _Ta = []
    for i,j in zip(vlsr,Ta):
        is_use= (i>=vmin)&(i<=vmax)
        _vlsr+= [i[is_use]]
        _Ta+= [j[:,is_use]]

    vlsr_len= np.array([len(i) for i in _vlsr])
    len_use= vlsr_len.min()

    vlsr= np.vstack([i[:len_use] for i in _vlsr])
    Ta= np.vstack([i[:,:len_use] for i in _Ta])
    return vlsr,Ta

def remove_nan(vel, Ta):
    is_use= np.sum(np.isnan(Ta),axis=0) < Ta.shape[0]/2
    return vel[:,is_use],Ta[:,is_use]

def _tight_ra_old(ra):
    
    is_= (ra>180)&(ra<=360)
    ind_= np.where(is_)[0]
    if len(ind_)==0 or len(ind_)==len(ra):
        return ra
    vmin1,vmax1= ra[~is_].min(), ra[~is_].max()
    vmin2,vmax2= ra[is_].min(), ra[is_].max()
    if vmin2-vmax1 > vmin1+360-vmax2:
        ra= np.copy(ra)
        ra[is_]= ra[is_]-360
    return ra

def _tight_ra(ra):
    ra_s= np.sort(ra)
    diff_s= np.diff(ra_s)
    ind_max= np.argmax(diff_s)
    if diff_s[ind_max] < ra_s[0]+360-ra_s[-1]:
        return ra
    else:
        ra= np.copy(ra)
        is_c= ra>=ra_s[ind_max+1]
        ra[is_c]= ra[is_c]-360
        return ra
    
if __name__ == '__main__':
    import os
    import argparse
    parser = argparse.ArgumentParser(allow_abbrev=False)
    parser.add_argument('fname',nargs='+',
                        help='file name')
    parser.add_argument('--ra_range', type=float,nargs=2,
                       help='ra range; unit: deg')
    parser.add_argument('--dec_range', type=float,nargs=2,
                       help='dec_range; unit: deg')
    parser.add_argument('--bwidth', type=float, default=[90.],nargs='+',
                       help='unit: arc second')
    parser.add_argument('--outname', required=True,
                       help='output file name, full path')
    parser.add_argument('--key', default='Ta',
                       help='output file name, full path')

    args = parser.parse_args()
    fname= args.fname
    ra_range= args.ra_range
    dec_range= args.dec_range
    bwidth= args.bwidth
    outname= args.outname
    key= args.key
    if os.path.exists(outname):
        raise(OSError(f"{outname} already exists."))
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
    
    
    ra=[]
    dec=[]
    Ta=[]
    vlsr=[]

    for file in files:
        tdata= h5py.File(file,'r')
        ra+= [tdata['ra'],]
        dec+= [tdata['dec'],]
        Ta+= [tdata[key],]
        vlsr+= [tdata['vel']]

    ra= np.hstack(ra)
    ra= _tight_ra(ra)
    dec= np.hstack(dec)
    #     Ta=np.vstack(Ta)
    #     vlsr=np.vstack(vlsr)
    vlsr, Ta= _stack_spec(vlsr, Ta) #vlsr are in descending order.
    #vlsr, Ta= remove_nan(vlsr, Ta)
    # use some vlsr sample for each specta
    std_vlsr= np.std(vlsr,axis=0,dtype=np.float64) # single precision can be inaccurate
    if np.nanmax(std_vlsr)>0.2:
        raise(ValueError('vlsr dispersion (std) at same Ta order is two large'))
    else:
        print(f'vlsr dispersion (std) at same Ta order is between {np.nanmin(std_vlsr)} and {np.nanmax(std_vlsr)}.')
    vlsr_refine= np.mean(vlsr,axis=0,dtype=np.float64) # single precision can be inaccurate
    vlsr_refine= np.append(vlsr_refine- (vlsr_refine[1]-vlsr_refine[0])/2, vlsr_refine[-1]+(vlsr_refine[1]-vlsr_refine[0])/2)

    if ra_range is None:
        ra_range= [np.nanmin(ra), np.nanmax(ra)]
    if dec_range is None:
        dec_range= [np.nanmin(dec), np.nanmax(dec)]
    
    bwidth= [i/60/60 for i in bwidth]
    if len(bwidth)==1:
        bwidth=bwidth[0]
        
    print('ra range:',ra_range,' deg')
    print('dec range:',dec_range,' deg')
    print('bwidth:',bwidth,' deg')
    Ta_refine, ra_refine,dec_refine= bin_stat(ra,dec,Ta,[ra_range,dec_range],bwidth)
    
    #reverse ra and the value at ra; .T
    ra_refine=ra_refine[::-1] # ra invert
    Ta_refine= np.array([ii[::-1,:].T for ii in Ta_refine])
    
    print(f'Saving to {outname}.')
    fits_cubes(ra_refine, dec_refine, vlsr_refine,Ta_refine, outname)
