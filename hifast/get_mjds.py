#!/usr/bin/env python
# coding: utf-8

from .core.tcal_onoff import tcal_onoff
from .utils.io import save_specs_hdf5
from glob import glob
import numpy as np
import os
import re


# In[13]:


# nB=1
# fname_part= f'/data/inspur_disk06/fast_data/3047/M31_45_0/20191227/M31_45_0_MultiBeamOTF-M{nB:02d}_F_'


# In[432]:


def get_mjds(fname_part, outdir=None, read_all=False, t_step=1):
    mjd_delta= t_step*np.array(1.16508454084396362304687500e-05, dtype=np.float64)
    end_file= len(glob(fname_part+'*.fits'))
    
    if read_all:
        mjds=tcal_onoff(fname_part,1,end_file).mjds
    else:
        mjds_1=tcal_onoff(fname_part,1,1).mjds
        mjds_2=tcal_onoff(fname_part,end_file,end_file).mjds
        mjds= np.arange(mjds_1[0],mjds_2[-1]+mjd_delta,mjd_delta)
    
    if outdir is None:
        return {'mjd':mjds}
    else:
        outname= os.path.join(outdir, 
                f"{os.path.basename(fname_part)[:-1]}-specs_T-mjd-{os.path.basename(os.path.dirname(fname_part))}.hdf5")
        save_specs_hdf5(outname,{'mjd':mjds})


# In[478]:


def main(basedirs,outdir='./', band=None, **kargs):
    if band is None:
        band = ['F','W','N']
    band = ','.join(band)
    fnames=[]
    for basedir in basedirs:
        fnames+= glob(basedir+f'/**/*M01_[{band}]*0001.fits',recursive=True)
    fnames=list(set(fnames))
    fnames.sort()
    [print(i) for i in fnames]
    for fname in fnames:
        fname_part= re.sub("[0-9]{4}\.fits\Z",'',fname)
        print(fname_part)
        get_mjds(fname_part,outdir, **kargs)


# In[467]:


if __name__ == '__main__':
    import os
    import argparse
    parser = argparse.ArgumentParser(allow_abbrev=False)
    parser.add_argument('basedirs',nargs='+',
                        help='basedirs')
    parser.add_argument('--outdir', required=False, default='./',
                       help='output file name, full path')
    parser.add_argument('-a', '--read_all', action='store_true',
                       help='read all files to obtain mjd')
    parser.add_argument('--t_step', type=float, default=1,
                       help='time step; s')
    parser.add_argument('-b', '--band', nargs='+', choices=['F', 'W', 'N'],
                       help='freq band type, F, W, N')
    args = parser.parse_args()
    main(args.basedirs, args.outdir, args.band, read_all=args.read_all, t_step=args.t_step)
