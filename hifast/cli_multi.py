#!/usr/bin/env python
# coding: utf-8

import sys
import os
import re

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(allow_abbrev=False)
    parser.add_argument('fname',
                        help='file name')
    parser.add_argument('--outdir',
                       help='default is same with the input file')
    parser.add_argument('-f', '--force', action='store_true',
                        help='overwriting file if out file exists')
    ## rfi
    parser.add_argument('--pr', '--polar_rfi', action='store_true',
                       help='')
    parser.add_argument('--pr_s_sigma', type=float, default=5,
                       help='gaussian smooth size for spectra smoothing along time axis')
    parser.add_argument('--pr_times', type=float, default=6,
                       help='')
    parser.add_argument('--pr_times_s', type=float, default=1,
                       help='')
    parser.add_argument('--tr', '--time_rfi', action='store_true',
                       help='')
    parser.add_argument('--tr_method', default='smooth',
                       help=" 'folder': merge spectra every 'folder'; smooth: gaussian smooth spectra along time axis")
    parser.add_argument('--tr_folder', type=int, default=40, 
                       help='')
    parser.add_argument('--tr_s_sigma', type=float, default=5,
                       help='gaussian smooth size for spectra smoothing along time axis')
    parser.add_argument('--tr_n_continue', type=float, default=40,
                       help='gaussian smooth size for spectra smoothing along time axis')
    parser.add_argument('--tr_times', type=float, default=6.,
                       help='')
    parser.add_argument('--tr_times_s', type=float, default=1.5,
                       help='')
    parser.add_argument('--ext_add', type=int, default=0,
                       help='extend rfi range')
    parser.add_argument('--ext_frac', type=float, default=0.,
                       help='between 0 and 1, extend rfi range')
    
    ## flux calibration
    parser.add_argument('--flux', action='store_true',
                       help='flux calibration')
    parser.add_argument('-c', '--cali_fname',
                       help='quasar calibration file name')
    
    ## frame correct
    parser.add_argument('--fc', action='store_true',
                       help='frame correct')
    parser.add_argument('--frame', choices=['HELIOCENT', 'LSRK'], default='LSRK',
                       help='Velocity Rest Frames, HELIOCENT or LSRK')
    parser.add_argument('--keep_rfi', action='store_true',
                       help='keep rfi')
    parser.add_argument('--keep_polar', action='store_true',
                       help='keep two polarizations')
#     parser.add_argument('--key', nargs='+',
#                        help='properties to correct, flux, Ta or both')
    
    args = parser.parse_args()
    fname= args.fname
    outdir= args.outdir
    # rfi
    pr = args.pr
    tr = args.tr
    if pr:
        pr_kwargs = {}
        pr_kwargs['s_sigma'] = args.pr_s_sigma
        pr_kwargs['times'] = args.pr_times
        pr_kwargs['times_s'] = args.pr_times_s
        pr_kwargs['ext_add'] = args.ext_add
        pr_kwargs['ext_frac']= args.ext_frac
    if tr:
        tr_kwargs = {}
        tr_kwargs['method'] = args.tr_method
        if tr_kwargs['method'] == 'folder':
            tr_kwargs['folder'] = args.tr_folder
        if tr_kwargs['method'] == 'smooth':
            tr_kwargs['n_continue'] = args.tr_n_continue
            tr_kwargs['times_s'] = args.tr_times
        tr_kwargs['times'] = args.tr_times_s
        tr_kwargs['ext_add'] = args.ext_add
        tr_kwargs['ext_frac']= args.ext_frac
    #keys= args.key
    # flux
    flux = args.flux
    cali_fname=args.cali_fname
    # frame
    fc = args.fc
    frame = args.frame
    keep_rfi = args.keep_rfi
    keep_polar = args.keep_polar
    
    outparts = []
    if pr or tr: outparts += ['rfi']
    if flux: outparts += ['flux']
    if fc: outparts += ['fc']
    outpart = '-' + '_'.join(outparts)
    if outdir is None: outdir = os.path.dirname(fname)
    fileout = os.path.join(outdir, '.'.join(os.path.basename(fname).split('.')[:-1]) + f'{outpart}.hdf5')
    if os.path.exists(fileout):
        if args.force:
            print(f"will overwrite the existing out file {fileout}")
        else:
            print(f"File exists {fileout}")
            print('exit... Using -f to overwrite it.')
            sys.exit()
    

    # load data
    import h5py
    import numpy as np
    f = h5py.File(fname, 'r')
    if 'T' in f.keys():
        T = f['T'][()]
        outfield = 'T'
    elif 'Ta' in f.keys():
        T = f['Ta'][()]
        outfield = 'Ta'
    elif 'flux' in f.keys():
        T = f['flux'][()]
        outfield = 'flux'
        if flux:
            raise(ValueError(f'flux already exists, please remove \"--flux\"'))
    if flux:
        outfield = 'flux'
            
    freq= f['freq'][()]
    mjd= f['mjd'][()]
    ra= f['ra'][()]
    dec= f['dec'][()]
    ra[ra<0] += 360.
    
    if 'frame' in f['Header'].attrs.keys():
            raise(ValueError(f'rest frame already corrected, please remove \"--fc\"'))
    if 'Header' in f.keys():
        from collections import OrderedDict
        header_in= OrderedDict(f['Header'].attrs.items())
    else:
        header_in=None
    # f.close()
    
    #do correct
    is_rfi = f['is_rfi'][()] if 'is_rfi' in f.keys() else None    
    if pr or tr:
        if is_rfi is None:
            is_rfi = np.full(T.shape[:2], False, dtype=bool)
    if pr:
        from .rfi_polar import mask_rfi_p
        is_rfi |= mask_rfi_p(T, **pr_kwargs)
    if tr:
        from .rfi_t import mask_rfi_t
        is_rfi |= mask_rfi_t(freq, T, **tr_kwargs)
    if flux:
        from .flux import cali_src
        nB = int(re.findall(r'-M[0-1][0-9]',fname)[-1][2:])
        T = cali_src(T, nB, freq, cali_fname, ra=ra, dec=dec, mjd=mjd)
    if fc:
        from .corr_vel import freq2vel, frame_correct
        if is_rfi is not None:
            if keep_rfi:
                is_rfi, _ = frame_correct(is_rfi, freq, mjd, ra, dec, frame=frame, interp_kind='nearest')
                is_rfi = np.array(is_rfi, dtype=bool)
            else:
                T[is_rfi] = np.nan
        if not keep_polar:
            T = np.mean(T, axis=2, dtype='float64')
        T, freq = frame_correct(T, freq, mjd, ra, dec, frame=frame)
        vel = freq2vel(freq)
    # out dict
    dict_out={}
    if fc:
        dict_out['vel'] = vel
    if is_rfi is not None and not fc:
        dict_out['is_rfi'] = is_rfi
    dict_out['freq'] = freq
    dict_out['ra'] = ra
    dict_out['dec'] = dec
    dict_out['mjd'] = mjd
    dict_out[outfield]= T.astype('float32')
    
    #save file
    print('Saving...')
    from .util import add_extra
    add_extra(f, dict_out)
    f.close()
    from .util import rec_his, save_dict_hdf5
    header=rec_his(args=args)
    if fc:
        header['frame']= frame
        header['vel_type']= 'VRAD'
    if header_in is not None: header.update(header_in)
    save_dict_hdf5(fileout, dict_out, header=header)
    print(f"Saved to {fileout}")
    