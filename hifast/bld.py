#!/usr/bin/env python
# coding: utf-8

import sys
import os
import re

if __name__ == '__main__':
    import warnings 
    warnings.filterwarnings("ignore",r'overflow encountered in exp')
    import argparse
    parser = argparse.ArgumentParser(allow_abbrev=False)
    parser.add_argument('fname',
                        help='file name')
    parser.add_argument('-n', '--nproc', type=int, default=1,
                        help='number process')
    parser.add_argument('-f', '--force', action='store_true',
                        help='overwriting file if out file exists')
    parser.add_argument('--frange', type=float, nargs=2,
                       help='freq range')
    # interaction
    parser.add_argument('-i', '--interact', action='store_true',
                       help='interaction')
    parser.add_argument('--ylim', nargs='+', default=['auto'],
                       help='ylim')
    parser.add_argument('--figsize', type=float, nargs=2, default=(10,7),
                       help='figsize')
    parser.add_argument('--length', type=int, nargs=1, default=20,
                       help='spetra numbe used to test')
    
    parser.add_argument('--nB_radec', type=int, default=1, 
                       help='Beam number of the radec file name')
    parser.add_argument('--outdir',
                       help='default is same with the input file')
    parser.add_argument('--flux', action='store_true',
                        help='convert temperature to flux before subtract baseline')
    parser.add_argument('-c', '--cali_fname',
                       help='quasar calibration file name')
    
    parser.add_argument('--method', default='arPLS', choices=['arPLS', 'srPLS', 'Chebyshev', 'poly', 'sin_poly', 'sin_poly_2', 'S', 'SP', 'original', 'asPLS'],
                       help='method to fit baseline')
    
    parser.add_argument('--lam', type=float, default=1.0e8, 
                       help='baseline fit parameters')
    parser.add_argument('--deg', type=int, default=2, 
                       help='baseline fit parameters')
    parser.add_argument('--offset', type=float, default=2, 
                       help='baseline fit parameters')
    parser.add_argument('--ratio', type=float, default=0.01, 
                       help='baseline fit parameters')
    parser.add_argument('--niter', type=int, default=100, 
                       help='baseline fit parameters')
    parser.add_argument('--sin_f', type=float, nargs='+',
                       help='sin freq')
    parser.add_argument('--s_method_freq', choices=['median', 'gaussian', 'boxcar', 'PLS', 'fft'],
                       help='')
    parser.add_argument('--s_sigma_freq', type=float,
                       help='')
    parser.add_argument('--average_every_freq', type=int,
                       help='')
    
    parser.add_argument('--lam_2', type=float, default=1.0e12, 
                       help='baseline fit parameters')
    parser.add_argument('--deg_2', type=int, default=2, 
                       help='baseline fit parameters')
    parser.add_argument('--offset_2', type=float, default=2, 
                       help='baseline fit parameters')
    parser.add_argument('--ratio_2', type=float, default=0.01, 
                       help='baseline fit parameters')
    parser.add_argument('--niter_2', type=int, default=100, 
                       help='baseline fit parameters')
    parser.add_argument('--sin_f_2', type=float, nargs='+',
                       help='sin freq')
    parser.add_argument('--s_method_freq_2', choices=['median', 'gaussian', 'boxcar', 'PLS'],
                       help='')
    parser.add_argument('--s_sigma_freq_2', type=int,
                       help='')
    parser.add_argument('--average_every_freq_2', type=int,
                       help='')
    
    parser.add_argument('--njoin', type=int,
                       help='join nspec along t')
    parser.add_argument('--njoin_2', type=int,
                       help='join nspec along t ')
    parser.add_argument('--s_method_t', choices=['gaussian', 'boxcar', 'median'],
                       help='')
    parser.add_argument('--s_sigma_t', type=int, default=5,
                       help='')
    parser.add_argument('--s_method_t_2', choices=['gaussian', 'boxcar', 'median'],
                       help='')
    parser.add_argument('--s_sigma_t_2', type=int, default=5,
                       help='')
    parser.add_argument('-T', '--trans', action='store_true',
                       help='')
    parser.add_argument('--exclude_m', type=int, default=0,
                       help='')
    parser.add_argument('--no_radec', action='store_true',
                       help='')
    
    args = parser.parse_args()
    
    if args.interact:
        import h5py
        from . import _bld_i
        fs = h5py.File(args.fname,'r')
        if 'T' in fs.keys():
            T = fs['T']
        elif 'Ta' in fs.keys():
            T = fs['Ta']
        elif 'flux' in fs.keys():
            T = fs['flux']
        _bld_i.T2p = T
        _bld_i.freq = fs['freq'][:]
        _bld_i.frange = args.frange
        _bld_i.nproc = args.nproc
        _bld_i.length = args.length
        _bld_i.figsize = args.figsize
        _bld_i.ylim = args.ylim[0] if len(args.ylim) ==1 else args.ylim
        _bld_i.main()
        sys.exit()
    
    nproc = args.nproc
    file_spec = args.fname
    outdir = args.outdir
    trans= args.trans
    flux = args.flux
    # determine output filename and check it  
    fpart = '-flux_' if args.flux else '-'
    fpart += 'bld'
    if args.trans: fpart += '_T'
    if outdir is None: outdir = os.path.dirname(file_spec)
    fileout = os.path.join(outdir, '.'.join(os.path.basename(file_spec).split('.')[:-1]) + f'{fpart}.hdf5')
    if os.path.exists(fileout):
        if args.force:
            print(f"will overwrite the existing out file {fileout}")
        else:
            print(f"File exists {fileout}")
            print('exit... Using -f to overwrite it.')
            sys.exit()
    # import 
    import numpy as np
    import h5py
    from scipy import ndimage
    from .core.baseline import get_baseline, get_baseline_mp, sub_baseline
    from .utils.misc import extend_Trues, average_every_n, smooth1d, smooth1d_fft
    exclude_m = args.exclude_m
    no_radec = args.no_radec
    # 
    frange = args.frange
    nB_radec = args.nB_radec
    cali_fname = args.cali_fname
    
    method = args.method
    njoin = args.njoin
    njoin_2 = args.njoin_2
    s_method_t = args.s_method_t
    s_sigma_t = args.s_sigma_t
    s_method_t_2 = args.s_method_t_2
    s_sigma_t_2 = args.s_sigma_t_2
    #smooth before fit baseline
    s_method_freq = args.s_method_freq
    s_sigma_freq = args.s_sigma_freq
    average_every_freq = args.average_every_freq
    s_method_freq_2 = args.s_method_freq_2
    s_sigma_freq_2 = args.s_sigma_freq_2
    average_every_freq_2 = args.average_every_freq_2
    
    #fit args
    fit_args={}
    if method in ['arPLS', 'srPLS', 'S', 'SP']:
        fit_args['lam'] = args.lam
    if method == 'sin_poly' or method == '1':
        fit_args['f'] = args.sin_f
        fit_args['rew'] = False  # only fit once
    fit_args['offset'] = args.offset
    fit_args['deg'] = args.deg
    fit_args['ratio'] = args.ratio
    fit_args['niter'] = args.niter
    
    if method in ['S', 'SP']:
        fit_args_2 = {}
        fit_args_2['sin_f'] = args.sin_f_2
        fit_args_2['offset'] = args.offset_2
        fit_args_2['deg'] = args.deg_2
        fit_args_2['ratio'] = args.ratio_2
        fit_args_2['niter'] = args.niter_2
        fit_args_2['rew'] = False  # only fit once
    
     
    nB = int(re.findall(r'-M[0-1][0-9]',file_spec)[-1][2:]) # used in flux cali
    fs = h5py.File(file_spec,'r')
    mjd = fs['mjd'][()]
    freq = fs['freq'][:]
    if not no_radec:
        if 'ra' not in fs.keys():
            from add_radec import get_radec
            ra, dec, is_extrapo = get_radec(fs, nB_radec)
        else:
            ra = fs['ra'][()]
            dec = fs['dec'][()]
            is_extrapo = None
    if 'T' in fs.keys():
        T = fs['T'][()]
        outfield = 'Ta'
    elif 'Ta' in fs.keys():
        T = fs['Ta'][()]
        outfield = 'Ta'
    elif 'flux' in fs.keys():
        T = fs['flux'][()]
        outfield = 'flux'
        if flux:
            raise(ValueError(f'flux already exists, please remove \"--flux\"'))
    if flux:
        outfield = 'flux'
    if frange is not None:
        is_ = (freq >= frange[0]) & (freq <= frange[1])
        freq = freq[is_]
        T = T[:,is_]

    # try to load header
    if 'Header' in fs.keys():
        from collections import OrderedDict
        header_in= OrderedDict(fs['Header'].attrs.items())
        try:
            header_in.update(OrderedDict(fs['Header'].attrs.items()))
        except:
            print('radec file has no header')
    else:
        header_in = None
    #close files
    # fs.close() # don't close, load extra later
    # flux cali
    if flux:
        from .core.flux import cali_src
        print('Flux calibrating ...')
        T = cali_src(T, nB, freq, cali_fname, ra=ra, dec=dec, mjd=mjd)
    
    if trans:
        T = T.transpose((1,0,2))
        freq_bak = freq
        freq = np.arange(T.shape[1])
    ## subtract baseline
    verbose = True
    if method == 'S' or method == 'SP':
        method_a = 'arPLS'
        para = {
        's_sigma_freq': s_sigma_freq,
        's_method_freq': s_method_freq,
        'average_every_freq': average_every_freq,
        }
        para.update(fit_args)
        
        method_b = 'sin_poly'
        bounds = [(0.,1.), (0.899, 0.961), (0, 2*np.pi), (-1,1)]
        fit_args_2['opt_para'] = {'bounds':bounds,}
        para2 = {
        's_sigma_freq': s_sigma_freq_2,
        's_method_freq': s_method_freq_2,
        'average_every_freq': average_every_freq_2,
        }
        para2.update(fit_args_2)
        T_bld1 = sub_baseline(freq, T, njoin=njoin, nproc=nproc, s_method_t=s_method_t, s_sigma_t=s_sigma_t, method=method_a, **para)
        #exclude_fun = lambda x:abs(x) > 1.2*np.diff(np.percentile(x, [16,84], axis=1), axis=0)[0][:,None,:]
        if exclude_m == 0:
            exclude_fun = lambda x: extend_Trues(abs(x) > 1.*np.diff(np.percentile(x, [16,84], axis=1), axis=0)[0][:,None,:], axis=1, ext_frac=1/3, ext_add=3)
        elif exclude_m == 1:
            exclude_fun = lambda x: extend_Trues(abs(x) > 2.5*np.min(np.diff(np.percentile(x, [10, 50, 90], axis=1), axis=0), axis=0)[:,None,:], axis=1, ext_frac=1/3, ext_add=3)
        bl_sin = T_bld1 - sub_baseline(freq, T_bld1, njoin=njoin_2, nproc=nproc, s_method_t=s_method_t_2, s_sigma_t=s_sigma_t_2, exclude_fun=exclude_fun, method=method_b, **para2)
        del(T_bld1)
        T = T - bl_sin
        del(bl_sin)
        if method == 'SP':
            T = sub_baseline(freq, T, njoin=njoin, method=method_a, **para)
    elif method == '1':
        pass
    else:
        para = {
        's_sigma_freq': s_sigma_freq,
        's_method_freq': s_method_freq,
        'average_every_freq': average_every_freq,
        }
        para.update(fit_args)
        T = sub_baseline(freq, T, njoin=njoin, nproc=nproc, s_method_t=s_method_t, s_sigma_t=s_sigma_t, method=method, **para)
        
    if trans:
        T = T.transpose((1,0,2))
        freq = freq_bak
    print(f"Saving...")
    dict_out= {}
    dict_out['mjd'] = mjd
    if not no_radec:
        dict_out['ra'] = ra
        dict_out['dec'] = dec
        if is_extrapo is not None:
            dict_out['is_extrapo'] = is_extrapo
    dict_out[outfield] = T.astype('float32')
    dict_out['freq'] = freq
    
    #save file
    from .utils.io import rec_his, save_dict_hdf5, add_extra
    add_extra(fs, dict_out)
    fs.close()
    header=rec_his(args=args)
    if header_in is not None: header.update(header_in)
    save_dict_hdf5(fileout, dict_out, header=header)
    print(f"Saved to {fileout}")
        
