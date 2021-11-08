#!/usr/bin/env python
# coding: utf-8

import numpy as np
from astropy import log
from hifast.utils.misc import down_sample
import os
import sys

def down(properties,keys,ra,dec,mjd,freq,vel,chanfactor=None,specfactor=None):
    
    if specfactor is not None:
        dec=down_sample(dec[None,:,None],specfactor)[0,:,0]
        ra=down_sample(ra[None,:,None],specfactor)[0,:,0]
        mjd=down_sample(mjd[None,:,None],specfactor)[0,:,0]
        for key in keys:
            properties[key]=down_sample(properties[key].T[:,:,None],specfactor)[:,:,0].T
        print(f"After spec down={properties[key].shape}")
    
    if chanfactor is not None:
        vel= down_sample(vel[None,:,None],chanfactor)[0,:,0]
        freq=down_sample(freq[None,:,None],chanfactor)[0,:,0]
        for key in keys:
            properties[key]=down_sample(properties[key][:,:,None],chanfactor)[:,:,0]
        print(f"After channel down={properties[key].shape}")
    return properties,ra,dec,mjd,freq,vel
    
if __name__ == '__main__':
    
    import argparse
    parser = argparse.ArgumentParser(allow_abbrev=False)
    parser.add_argument('fname',
                        help='file name')
    parser.add_argument('--outdir',
                       help='default is same with the input file')

    parser.add_argument('--key', nargs='+',
                       help='properties to down sample, flux, Ta or both')

    parser.add_argument('--chanfactor',type=int,
                       help=' down sample on channel')

    parser.add_argument('--specfactor',type=int,
                       help=' down sample on spec')

    parser.add_argument('-f', '--force', action='store_true',
                       help='overwriting file if out file exists')

    args = parser.parse_args()
    fname= args.fname
    outdir= args.outdir
    keys= args.key
    chanfactor= args.chanfactor
    specfactor= args.specfactor

    #save
    if outdir is None:
        outdir= os.path.dirname(fname)
    fileout= os.path.join(outdir, '.'.join(os.path.basename(fname).split('.')[:-1]) +'-ds.hdf5' )

    if os.path.exists(fileout):
        if args.force:
            print(f"will overwrite the existing out file {fileout}")
        else:
            print(f"File exists {fileout}")
            print('exit... Using -f to overwrite it.')
            sys.exit()

    import h5py
    f= h5py.File(fname,'r')
    properties={}
    if keys is None:
        keys=[]
        for key in ['Ta', 'flux', ]:
            if key in f.keys():
                keys+=[key]
    for key in keys:
        T = f[key][()]
        if len(T.shape) == 3:
            T = np.mean(T, axis=2, dtype='float64') 
        properties[key]= T
    print(f"Origin shape={properties[key].shape}")

    freq= f['freq'][()]
    mjd= f['mjd'][()]
    ra= f['ra'][()]
    dec= f['dec'][()]
    vel= f['vel'][()]

    if 'Header' in f.keys():
        from collections import OrderedDict
        header_in= OrderedDict(f['Header'].attrs.items())
    else:
        header_in=None


    #down sample
    properties,ra,dec,mjd,freq,vel = down(properties,keys,ra,dec,mjd,freq,vel,chanfactor,specfactor)

    dict_out={}
    dict_out['freq']=freq
    dict_out['vel']=vel
    dict_out['ra']=ra
    dict_out['dec']=dec
    dict_out['mjd']=mjd
    for key in properties.keys():
        dict_out[key]= properties[key].astype('float32')

    from hifast.utils.io import  rec_his, save_dict_hdf5
    f.close()
    header=rec_his(args=args)
    
    if header_in is not None: header.update(header_in)
    save_dict_hdf5(fileout, dict_out, header=header)
    log.info(f"Saved to {fileout}")
