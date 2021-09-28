#!/usr/bin/env python
# coding: utf-8

# author: Xu Chen, Li Fujia, 2021.06
# code：Xu Chen
# original code: Jing Yingjie

import sys
import os
import re

from astropy import log

import numpy as np
import h5py
from copy import deepcopy
from tqdm import tqdm
from glob import glob


if __name__ == '__main__':
    import warnings 
    warnings.filterwarnings("ignore",r'overflow encountered in exp')
    warnings.filterwarnings("ignore",r'Polyfit may be poorly conditioned')
    
    import argparse
    parser = argparse.ArgumentParser(allow_abbrev=False)
    
    parser.add_argument('fname',
                        help='file name')
    parser.add_argument('--outdir',
                       help='default is same with the input file')
    parser.add_argument('-f', '--force', action='store_true',
                        help='overwriting file if out file exists')
    parser.add_argument('--frange', type=float, nargs=2,
                       help='freq range')
    
    parser.add_argument('--rms_frange', type=float, nargs=2,
                       help='freq range to compute rms')
    parser.add_argument('--rms_sigma', type=float, default =6,
                       help='gauss filter sigma to compute real rms')
    parser.add_argument('--mw_frange', type=float, nargs=2,
                       help='milky way freq range')  
    ## time rfi
    parser.add_argument('--time_rfi', action='store_true',
                        help='find time rfi')
    parser.add_argument('--lf', action='store_true',
                        help='find long time time rfi')
    parser.add_argument('--sf', action='store_true',
                        help='find short time time rfi')
    # short freq
    parser.add_argument('--sf_frange', type=float, nargs=2,
                       help='freq range exists short-freq time rfi ')
    parser.add_argument('--sf_file',         
                        help='freq range exists short-freq time rfi npy filename')
    parser.add_argument('--sf_frange_step',type = int,default=20,
                        help='if sf_frange is None and sf_file is None, cycle in whole freq band.')
    
    parser.add_argument('--sf_times', type=float, default=3,
                       help='first threhold, rfi is this times of median value')
    parser.add_argument('--sf_thr', type=float, default=1,
                       help='sharp edge on time axis. diff above this times of next point will be recognized.')
    parser.add_argument('--sf_rfi_last',type=float, default=10, 
                       help='rfi lasts at least 20 spec numbers')
    parser.add_argument('--sf_T_thr_times',type=float, default=3, 
                       help='T above thr will be masked (default 2 times RMS)')
    parser.add_argument('--sf_ext',type = int,default=0, 
                       help='extend edge')
    # long freq
    parser.add_argument('--lf_beams', 
                        help='beam numbers which has long-freq time rfi')
    parser.add_argument('--lf_sepname', 
                        help='search lf in sep data? If it is None, search in sub data input.')
    parser.add_argument('--lf_frange', type=float, nargs=2,
                       help='freq range exists long-freq time rfi')
    parser.add_argument('--lf_times', type=float, default=1.5,
                       help='first threhold, rfi is this times of median value')
    parser.add_argument('--lf_thr', type=float, default=0,
                       help='set 0 and do not change')
    parser.add_argument('--lf_rfi_last',type=float, default=50, 
                       help='rfi lasts at least 20 spec numbers')
    parser.add_argument('--lf_ext',type = int,default=1, 
                       help='extend edge')

    
    ## flux calibration
    parser.add_argument('--flux', action='store_true',
                       help='flux calibration')
    parser.add_argument('-c', '--cali_fname',
                       help='quasar calibration file name')
    parser.add_argument('--keep_polar',action='store_true',
                       help='keep polar')
    
    parser.add_argument('--keep_rfi', action='store_true',
                       help='keep rfi')
    parser.add_argument('--plot', action= 'store_true',
                       help='plot')
    parser.add_argument('--ylim', type=float, nargs=2,
                        help='set ylim in plot')
    parser.add_argument('--vmin_max', type=float, nargs=2,
                        help='plot waterfall')
    
    
    args = parser.parse_args()
    fname = args.fname
    outdir = args.outdir
    frange = args.frange
    # flux
    flux = args.flux
    cali_fname=args.cali_fname
    
    keep_polar = args.keep_polar
    keep_rfi = args.keep_rfi
    plot = args.plot
    
    time_rfi = args.time_rfi
    longf_rfi = args.lf
    shortf_rfi = args.sf
    
    if (longf_rfi == True) or (shortf_rfi == True):
        time_rfi = True
    
    outparts = []
    if time_rfi:outparts += ['tr']
    if flux: outparts += ['flux']
    
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
        
    if plot:
        from matplotlib import pyplot as plt
        plt.switch_backend('agg')
        from matplotlib.backends.backend_pdf import PdfPages
        if not os.path.exists(outdir+'/fig/'):
            os.mkdir(outdir+'/fig/')
        pdfname = os.path.join(outdir+'/fig/','.'.join(os.path.basename(fileout).split('.')[:-1])+ '.pdf')
        pdf = PdfPages(pdfname)
    else:
        pdf = None

    # load data
    import h5py
    import numpy as np
    f = h5py.File(fname, 'r')
    mjd= f['mjd'][()]
    ind_sort = np.argsort(f['mjd'])
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
        
    if len(T.shape) == 3:
        if keep_polar:
            T3 = deepcopy(T)
        T = np.mean(T, axis=2, dtype='float64')    
    
    T = T[ind_sort]           
    freq= f['freq'][()]
    
    # fdelta = freq[1] - freq[0]
    if frange is not None:
        is_ = (freq >= frange[0]) & (freq <= frange[1])
        freq = freq[is_]
        T = T[:,is_]
            
    if 'ra' in f.keys():
        ra= f['ra'][()][ind_sort]
        dec= f['dec'][()][ind_sort]
        ra[ra<0] += 360.
    
    if 'Header' in f.keys():
        from collections import OrderedDict
        header_in= OrderedDict(f['Header'].attrs.items())
    else:
        header_in=None
        
    T_ori = deepcopy(T)

    
    # load data
    sep_fname = args.lf_sepname 
    if sep_fname == 'input_subname':
        sep_fname = fname
        
    elif sep_fname is None:
        # find higher class folder
        sep_dir_default = os.path.dirname(fname).split('/')[:-1]
        sep_dir_default.append('sep')
        sep_dir_default = '/'.join(sep_dir_default)
        sep_dir = os.path.join(sep_dir_default, '.'.join(os.path.basename(fname).split('-bld',1)[:-1]) + '.hdf5')
        
        sep_fnames = glob(sep_dir)
        if len(sep_fnames) == 1:
            sep_fname = sep_fnames[0]
        elif len(sep_fnames) == 0:
            #find in same class folder
            sep_fname = os.path.join(os.path.dirname(fname),'.'.join(os.path.basename(fname).split('-bld',1)[:-1]) +'.hdf5')
        else:
            raise FileNotFoundError("Which sep file do you want? ")
    else:
        sep_fnames = glob(sep_fname)
        if len(sep_fnames) == 1:
            sep_fname = sep_fnames[0]
        else:
            filedir = os.path.dirname(sep_fname)
            sep_fname = os.path.join(filedir, os.path.basename(fname).split('-bld')[0] + '.hdf5')
    
    if sep_fname != fname:
        sep = h5py.File(sep_fname,'r')
        print(f"Find sep file {sep_fname}to find long-freq time rfi.") 
        from util import get_data
        T_sep = get_data(sep,polar = 'merged',xrange = frange)
        
        if len(T_sep.shape) == 3:
            T_sep = np.mean(T_sep, axis=2, dtype='float64')  
        if T_sep.shape != T.shape:
            raise ValueError(f"sep data shape {T_sep.shape} is not matched with sub data shape {T.shape}.")
    else:
        T_sep = deepcopy(T_ori)
        print(f"Use the input sub file to find long-freq time rfi.") 
        
    from markRFI import real_rms
    rms_frange = args.rms_frange
    rms_sigma = args.rms_sigma
    
    mw_frange = args.mw_frange

    if mw_frange is None:
        protect_use = np.zeros_like(freq,dtype='bool')
    else:
        protect_use = (freq>mw_frange[0])&(freq<mw_frange[1])
        
    t_rfi = np.isnan(T)
    Tt = deepcopy(T)
    Tt[:,protect_use] = 0  
    T_sep[:,protect_use] = 0

    ####################### time RFI ################################
    
    if time_rfi:
        lf_beams = args.lf_beams
        if lf_beams is not None:
            m_pos = fname.find('-M')
            nB = fname[m_pos+2:m_pos+4]
            if nB not in lf_beams:
                longf_rfi = False
            else:
                print(f"M{nB} in beam {lf_beams}")
  
        from markRFI import mask_time_rfi
        if longf_rfi:
            longf_args = {}
            longf_args['frange'] = args.lf_frange
            longf_args['times'] = args.lf_times
            longf_args['thr'] = args.lf_thr
            longf_args['rfi_width_lim'] = args.lf_rfi_last
            longf_args['ext_add'] = args.lf_ext
            print("long freq args:",longf_args)
            l_rfi = mask_time_rfi(T_sep,freq, rtype = 'long-freq',plot = plot,pdf = pdf,
                                  **longf_args)
            t_rfi = t_rfi | l_rfi
            
    whole_rfi = np.all(t_rfi,axis = 1)
    not_rfi_num = np.arange(T.shape[0])[~whole_rfi]
    is_rfi_num = np.arange(T.shape[0])[whole_rfi]
    
    if time_rfi:
        frange_step = args.sf_frange_step
        
        if shortf_rfi:
            tn = np.argmin(np.abs(np.nansum(data[not_rfi_num,:],axis = 1)))
            spec = deepcopy(T[tn,:])
            RMS = real_rms(spec,freq,rms_sigma,rms_frange)
            
            shortf_args = {}
            shortf_args['file'] = args.sf_file
            shortf_args['frange'] = args.sf_frange
            shortf_args['times'] = args.sf_times
            shortf_args['thr'] = args.sf_thr
            shortf_args['rfi_width_lim'] = args.sf_rfi_last
            shortf_args['T_thr_times'] = args.sf_T_thr_times
            shortf_args['ext_add'] = args.sf_ext
            print("short freq args:",shortf_args)
            s_rfi = mask_time_rfi(Tt,freq, rtype = 'short-freq',plot = plot,pdf = pdf,RMS = RMS,
                                  frange_step =frange_step, **shortf_args)
            
            t_rfi = t_rfi | s_rfi

    if time_rfi:
        rfi_mask = deepcopy(t_rfi)
    else:
        raise ValueError("Are you really want to mask time RFI?")
    
    
    is_rfi_ori = f['is_rfi'][()] if 'is_rfi' in f.keys() else None 
    if is_rfi_ori is not None:
        if time_rfi:
            rfi_mask = rfi_mask | is_rfi_ori
            log.warning("'is_rfi' alredy exists, found again and use union set.")
        else:
            rfi_mask = is_rfi_ori
            print("Use 'is_rfi' in itself.")

    if np.sum(rfi_mask) > 0:
        if not keep_rfi:
            if keep_polar: 
                T_ret = deepcopy(T3)
                T_ret[rfi_mask,:] = np.nan
            else:
                T_ret = deepcopy(T_ori)
                T_ret[rfi_mask] = np.nan  
        else:
            if keep_polar: 
                T_ret = deepcopy(T3)
            else:
                T_ret = deepcopy(T_ori)

        if plot:    
            print(" 'Wait for plotting patiently, you must.' Master Yoda said.")
            vmin_max = args.vmin_max
            from util import plot_waterfall
            if keep_polar:
                plot_waterfall(f,data = T_ret[:,:,0], vmin_max=vmin_max,cmap='plasma',figsize=(20,5),
                               title = os.path.basename(fileout).split('.')[:-1][0] + '_polar xx',pdf = pdf)

                plot_waterfall(f,data = T_ret[:,:,1], vmin_max=vmin_max,cmap='plasma',figsize=(20,5),
                               title = os.path.basename(fileout).split('.')[:-1][0] + '_polar yy',pdf = pdf)
            else:
                plot_waterfall(f,data = T_ret, vmin_max=vmin_max,cmap='plasma',figsize=(20,5),
                               title = os.path.basename(fileout).split('.')[:-1][0] + '_polar merged',pdf = pdf)              

            pdf.close()
            log.info(f"Plot to {pdfname}")
    else:
        rfi_mask = np.zeros_like(T,dtype='bool')
        log.info("Nothing is masked.")
        if plot:
            pdf.close()
            os.remove(pdfname)
        T_ret = np.array(0.)

    if flux:
        from hifast.core.flux import cali_src
        nB = int(re.findall(r'-M[0-1][0-9]',fname)[-1][2:])
        T_ret = cali_src(T_ret, nB, freq, cali_fname, ra=ra, dec=dec, mjd=mjd)
    
    # out dict
    dict_out={}
    dict_out['is_rfi'] = rfi_mask
    print("saved 'is_rfi'")
    dict_out['freq'] = freq
    if 'ra' in f.keys():
        dict_out['ra'] = ra
        dict_out['dec'] = dec
        
    dict_out['mjd'] = mjd
    dict_out[outfield] = T_ret.astype('float32')
    
    #save file
    print('Saving...')
    
    from hifast.utils.io import rec_his, save_dict_hdf5,add_extra
    add_extra(f, dict_out)
    f.close()
    if sep_fname != fname:
        sep.close()
    header=rec_his(args=args)

    if header_in is not None: header.update(header_in)
    save_dict_hdf5(fileout, dict_out, header=header)
    log.info(f"Saved to {fileout}")
    
