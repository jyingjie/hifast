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
    parser.add_argument('-rfi', '--rfi_fname',
                       help='mask rfi file name')
    
    parser.add_argument('--rms_frange', type=float, nargs=2,
                       help='freq range to compute rms')
    parser.add_argument('--rms_sigma', type=float, default =6,
                       help='gauss filter sigma to compute real rms')
    # find M31
    parser.add_argument('--m31_frange', type=float, nargs=2,
                       help='M31 freq range(estimate)')
    parser.add_argument('--m31_times', type=float, default=1.5,
                       help='first threhold, M31 is this times of median value')
    parser.add_argument('--m31_thr', type=float, default=0,
                       help='sharp edge on time axis. diff above this times of next point will be recognized.')
    parser.add_argument('--m31_rfi_last',type=float, default=50, 
                       help='m31 lasts at least 50 spec numbers')
    parser.add_argument('--m31_ext',type = int,default=10, 
                       help='extend edge')
    
    #smooth
    parser.add_argument('--s_method_freq', default='None', choices=['gaussian', 'boxcar', 'median', 'None'],
                       help='smooth spec along freq')
    parser.add_argument('--s_sigma_freq', type=int, default=3,
                       help='smooth spec along freq sigma')
    parser.add_argument('--s_method_t', default='None', choices=['gaussian', 'boxcar', 'median', 'None'],
                       help='smooth spec along time')
    parser.add_argument('--s_sigma_t', type=int, default=10,
                       help='smooth spec along time sigma')
   
    # find rfi
    parser.add_argument('--rfi_thr', type=float, default=3, 
                       help='default 3 times of rms threshold')
    parser.add_argument('--mw_frange', type=float, nargs=2,
                       help='protect milky way freq range(estimate)')
    parser.add_argument('--no_m31_frange', type=float, nargs=2,
                       help='no M31 freq range(estimate)')
    
    parser.add_argument('--rfi_width_lim', type=float, default=10, 
                       help='rfi should contain more channels than limit')
    parser.add_argument('--ext_sec', type=int, default=30, 
                       help='extend channel number of start and end of each section')
    parser.add_argument('--freq_thr', type=float, default=.5, 
                       help='estimate error when looking for peaks to polyfit')
    parser.add_argument('--freq_step', type=float, default = 8.1, 
                       help='rfi period MHz')
    parser.add_argument('--rfi_groups', default = 'two groups', choices=['two groups','three groups','all'],
                        help='divide rfi into 2 or 3 groups')
    
    ## mask rfi
    parser.add_argument('--freq_from_theory', type=float,default = 3,
                        help='mask width from theory center')
    parser.add_argument('--mask_thr', type=float, default=3, 
                       help='default below 3 times of rms threshold will be masked')
    parser.add_argument('--ext_edge', type=int, default= 0, 
                       help=' extend result channel on freq axis')
    parser.add_argument('--mask_RFI_method',default='fixed freq',choices=['2 sides','fixed freq'],
                         help='from center to two sides, or use a fixed freq width')
    # 2 sides  
    parser.add_argument('--small_rfi_times',type=float, default = 2,
                        help='small rfi below 2*RMS will not be masked')
    parser.add_argument('--chan_step' ,type=int, default= 5,
                        help='channel step when walk from center to two sides')
    # fixed freq
    parser.add_argument('--mask_all_theory', action= 'store_true',
                       help='mask_all_theory')
    ## 
    parser.add_argument('--only_M31', action='store_true',
                       help='only subtract M31 shape?')
    
    ## flux calibration
    parser.add_argument('--flux', action='store_true',
                       help='flux calibration')
    parser.add_argument('-c', '--cali_fname',
                       help='quasar calibration file name')
    #parser.add_argument('--keep_polar',action='store_true',
    #                   help='keep polar')
    
    
    parser.add_argument('--plot', action= 'store_true',
                       help='plot')
    parser.add_argument('--ylim', type=float, nargs=2,
                        help='set ylim in plot')

    
    args = parser.parse_args()
    fname = args.fname
    outdir = args.outdir

    # flux
    flux = args.flux
    cali_fname=args.cali_fname
    rfi_fname = args.rfi_fname
    #keep_polar = args.keep_polar

    plot = args.plot

    find_args = {}
    find_args['rfi_width_lim'] = args.rfi_width_lim
    find_args['ext_sec'] = args.ext_sec
    
    rfi_thr = args.rfi_thr
    freq_thr = args.freq_thr
    freq_step = args.freq_step
    ext_edge = args.ext_edge
    rfi_groups = args.rfi_groups

    mask_RFI_method = args.mask_RFI_method
    small_rfi_times = args.small_rfi_times
    chan_step = args.chan_step
    mask_all_theory = args.mask_all_theory
    freq_from_theory = args.freq_from_theory
    mask_thr = args.mask_thr
    
    from markRFI import real_rms
    rms_frange = args.rms_frange
    rms_sigma = args.rms_sigma
    
    outpart = '-shape'
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
        #if keep_polar:
        #    T3 = deepcopy(T)
        T = np.mean(T, axis=2, dtype='float64')    
    
    T = T[ind_sort]           
    freq= f['freq'][()]
    if 'ra' in f.keys():
        ra= f['ra'][()][ind_sort]
        dec= f['dec'][()][ind_sort]
        ra[ra<0] += 360.
    
    if 'Header' in f.keys():
        from collections import OrderedDict
        header_in= OrderedDict(f['Header'].attrs.items())
    else:
        header_in=None

    ###### read RFI
    if rfi_fname is not None:
        if rfi_fname == 'none':
            print("Don't use rfi file.")
    else:
        rfi_dir_default = os.path.dirname(fname).split('/')[:-1]
        rfi_dir_default.append('rfi')
        
        rfi_dir_default = '/'.join(rfi_dir_default)
        outpart = '*-tr*.hdf5'
        rfi_dir = os.path.join(rfi_dir_default, '.'.join(os.path.basename(fname).split('T')[:-1]) + f'{outpart}')
        from glob import glob
        rfi_fnames = glob(rfi_dir)
        if len(rfi_fnames) == 1:
            rfi_fname = rfi_fnames[0]
        elif len(rfi_fnames) == 0:
            rfi_dir = os.path.join(os.path.dirname(fname),'.'.join(os.path.basename(fname).split('T')[:-1]) + f'{outpart}')
            rfi_fnames = glob(rfi_dir)
            if len(rfi_fnames) == 1:
                rfi_fname = rfi_fnames[0]
            else:
                raise FileNotFoundError("Where is rfi mask?")
        else:
            raise FileNotFoundError("Which rfi mask do you want? ")
    
    if rfi_fname != 'none':
        rfi = h5py.File(rfi_fname,'r')
        print(f"Find rfi file {rfi_fname}")

        if 'short_rfi' in rfi.keys():
            s_rfi = rfi['short_rfi'][()]
        else:
            s_rfi = np.full(T.shape,False)
            log.warning("rfi file doesn't contain key 'short_rfi',so ignore short time rfi.")
        rfi.close()
        if len(s_rfi.shape) == 3:
            # use 2D rfi mask
            s_rfi = s_rfi[:,:,0]|s_rfi[:,:,1]

    T[s_rfi] = 0
    T_ori = deepcopy(T)
    ############################## find m31 #####################################
    m31_args = {}
    m31_args['frange'] = args.m31_frange
    m31_args['times'] = args.m31_times
    m31_args['thr'] = args.m31_thr
    m31_args['rfi_width_lim'] = args.m31_rfi_last
    m31_args['ext_add'] = args.m31_ext
    from markRFI import find_t
    m31_tuse = find_t(T,freq,plot = plot,pdf = pdf,**m31_args)
    
    m31_trange = np.where(m31_tuse==True)[0][0],np.where(m31_tuse==True)[0][-1]
    
    ########################### smooth ###############################3
    s_method_freq = args.s_method_freq
    s_sigma_freq = args.s_sigma_freq
    s_method_t = args.s_method_t
    s_sigma_t = args.s_sigma_t
    
    from hifast.utils.misc import smooth1d
    if  s_method_t in ['gaussian', 'boxcar', 'median']:
        print('Smooth ing ...')
        T = smooth1d(T,axis = 0,sigma = s_sigma_t, method = s_method_t)
    
    if  s_method_freq in ['gaussian', 'boxcar', 'median']:
        print('Smooth ing ...')
        T = smooth1d(T,axis = 1,sigma = s_sigma_freq, method = s_method_freq)
    
    ############################## subtract #####################################
    no_m31_frange = args.no_m31_frange
    
    mw_frange = args.mw_frange
    only_M31 = args.only_M31
    
    T_ret = deepcopy(T)
    
    mw_use = (freq>mw_frange[0])&(freq<mw_frange[1])
    m31_use = (freq>no_m31_frange[0])&(freq<no_m31_frange[1])
    
    from RFI_shape import rfis_in_one_spec,common_shape
    from markRFI import find_RFI
    
    if plot:    
        print(" 'Wait for plotting patiently, you must.' Master Yoda said.")
        tn = (m31_trange[0]+m31_trange[1])//2
        spec = deepcopy(T[tn,:]) 
        RMS = real_rms(spec,freq,rms_sigma,rms_frange)
        is_rfi = (spec > RMS * rfi_thr) & (~mw_use)
        pd_rfi,theorys = find_RFI(
            spec,freq,is_rfi,freq_step = freq_step,RMS = RMS,
            freq_thr = freq_thr,ext_edge = ext_edge,rfi_fit_use = rfi_groups ,plot = plot,pdf=pdf,
            mask_RFI_method = mask_RFI_method, small_rfi_times = small_rfi_times,chan_step = chan_step,
            mask_all_theory = mask_all_theory,freq_from_theory = freq_from_theory,mask_thr = mask_thr,
            ret_type='theory 123',**find_args)

        syn_rfis,syn_freqs = rfis_in_one_spec(spec,freq,pd_rfi,mw_use,freq_step=8.1,RMS = RMS,
                 rfi_fit_use = rfi_groups,step=0,plot = plot,pdf=pdf,freq_thr = freq_thr,ylim=[-1,5],
                  **find_args)

        _,com_rfi = common_shape(spec,T_ori[tn],freq,pd_rfi,theorys,syn_rfis,syn_freqs,
                 m31_use = m31_use,RMS=RMS,freq_thr=freq_thr,plot = plot,pdf=pdf,
                 only_M31=only_M31,ylim=[-1,5])

        
    log.info("Dealing with RFI near M31...")
    for tn in tqdm(range(m31_trange[0],m31_trange[1])):
        spec = deepcopy(T[tn,:]) 
        RMS = real_rms(spec,freq,rms_sigma,rms_frange)
        is_rfi = (spec > RMS * rfi_thr) & (~mw_use)
        try:
            pd_rfi,theorys = find_RFI(
                spec,freq,is_rfi,freq_step = freq_step,RMS = RMS,
                freq_thr = freq_thr,ext_edge = ext_edge,rfi_fit_use = rfi_groups ,plot = False,
                mask_RFI_method = mask_RFI_method, small_rfi_times = small_rfi_times,chan_step = chan_step,
                mask_all_theory = mask_all_theory,freq_from_theory = freq_from_theory,mask_thr = mask_thr,
                ret_type='theory 123',**find_args)
            
            syn_rfis,syn_freqs = rfis_in_one_spec(spec,freq,pd_rfi,mw_use,freq_step=8.1,RMS = RMS,
                     rfi_fit_use = rfi_groups,step=0,plot = False,freq_thr = freq_thr,ylim=[-1,5],
                                                  **find_args)
            
            T_ret[tn,:],com_rfi = common_shape(spec,T_ori[tn],freq,pd_rfi,theorys,syn_rfis,syn_freqs,
                      m31_use = m31_use,RMS=RMS,freq_thr=freq_thr,plot = False,
                     only_M31=only_M31,ylim=[-1,5])

        except ValueError:
            T_ret[tn,:] = np.nan#spec
            log.warning(f"tn={tn} has a ValueError !")
            import traceback
            traceback.print_exc()  
        except IndexError:
            log.warning(f"tn={tn} has a IndexError !")
            import traceback
            traceback.print_exc()  
            import sys
            sys.exit()

    log.info("Finish")

    if plot:    

        from util import plot_waterfall
        
        plot_waterfall(f,data = T_ret, vmin_max=[-.05,.05],cmap='plasma',figsize=(18,5),
                           title = os.path.basename(fileout).split('.')[:-1][0] + '_polar merged',pdf = pdf)              

        pdf.close()
        log.info(f"Plot to {pdfname}")

    if flux:
        from hifast.core.flux import cali_src
        nB = int(re.findall(r'-M[0-1][0-9]',fname)[-1][2:])
        T_ret = cali_src(T_ret, nB, freq, cali_fname, ra=ra, dec=dec, mjd=mjd)
    
    # out dict
    dict_out={}
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
    rfi.close()
    header=rec_his(args=args)

    if header_in is not None: header.update(header_in)
    save_dict_hdf5(fileout, dict_out, header=header)
    log.info(f"Saved to {fileout}")
    
