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
    
    parser.add_argument('--rms_frange', type=float, nargs=2,
                       help='freq range to compute rms')
    parser.add_argument('--rms_sigma', type=float, default =6,
                       help='gauss filter sigma to compute real rms')
    
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
    parser.add_argument('--sf_times', type=float, default=10,
                       help='first threhold, rfi is this times of median value')
    parser.add_argument('--sf_thr', type=float, default=10,
                       help='sharp edge on time axis. diff above this times of next point will be recognized.')
    parser.add_argument('--sf_rfi_last',type=float, default=20, 
                       help='rfi lasts at least 20 spec numbers')
    parser.add_argument('--sf_T_thr_times',type=float, default=2, 
                       help='T above thr will be masked (default 2 times RMS)')
    parser.add_argument('--sf_ext',type = int,default=0, 
                       help='extend edge')
    # long freq
    parser.add_argument('--lf_beams', 
                        help='beam numbers which has long-freq time rfi')
    parser.add_argument('--lf_frange', type=float, nargs=2,
                       help='freq range exists long-freq time rfi')
    parser.add_argument('--lf_times', type=float, default=1.5,
                       help='first threhold, rfi is this times of median value')
    parser.add_argument('--lf_thr', type=float, default=0,
                       help='sharp edge on time axis. diff above this times of next point will be recognized.')
    parser.add_argument('--lf_rfi_last',type=float, default=20, 
                       help='rfi lasts at least 20 spec numbers')
    parser.add_argument('--lf_ext',type = int,default=10, 
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
    
    ## freq period rfi
    parser.add_argument('--period_rfi', action='store_true',
                        help='find freq period rfi')
    # find rfi
    parser.add_argument('--rfi_thr', type=float, default=3, 
                       help='default 3 times of rms threshold')
    parser.add_argument('--mw_frange', type=float, nargs=2,
                       help='protect milky way freq range(estimate)')
    
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
    
    # time coherent
    parser.add_argument('--time_coherent_per', type=float, default = 1,
                       help='rfi in one freq appears more than emmm, maybe 70%, mask them all on time axis.')
    
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
    parser.add_argument('--save_rfi_list', action= 'store_true',
                       help='save rfi freq list in hdf5')
    parser.add_argument('--save_sf', action= 'store_true',
                       help='save short time rfi in hdf5')
    
    
    args = parser.parse_args()
    fname = args.fname
    outdir = args.outdir
    
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
    
    period_rfi = args.period_rfi
    if period_rfi:
        find_args = {}
        find_args['rfi_width_lim'] = args.rfi_width_lim
        find_args['ext_sec'] = args.ext_sec
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
        save_rfi_list = args.save_rfi_list
        if save_rfi_list:
            plot = True
    
    outparts = []
    if time_rfi:outparts += ['tr']
    if period_rfi:
        outparts += ['pdr']
        if mask_all_theory: outparts += ['strict']
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
    
    
    from markRFI import real_rms
    rms_frange = args.rms_frange
    rms_sigma = args.rms_sigma
    
    ####################### time RFI ################################
    t_rfi = np.zeros_like(T,dtype = 'bool')
    
    if time_rfi:
        lf_beams = args.lf_beams
        if lf_beams is not None:
            m_pos = fname.find('-M')
            nB = fname[m_pos+2:m_pos+4]
            if nB not in lf_beams:
                longf_rfi = False
            else:
                print(f"M{nB} in beam {lf_beams}")
  
        t_rfi = np.zeros_like(T,dtype = 'bool')
        from markRFI import mask_time_rfi
        if longf_rfi:
            longf_args = {}
            longf_args['frange'] = args.lf_frange
            longf_args['times'] = args.lf_times
            longf_args['thr'] = args.lf_thr
            longf_args['rfi_width_lim'] = args.lf_rfi_last
            longf_args['ext_add'] = args.lf_ext
            
            l_rfi = mask_time_rfi(T,freq, rtype = 'long-freq',plot = plot,pdf = pdf,**longf_args)
            t_rfi = l_rfi
            
    whole_rfi = np.all(t_rfi,axis = 1)
    not_rfi_num = np.arange(T.shape[0])[~whole_rfi]
    is_rfi_num = np.arange(T.shape[0])[whole_rfi]
    
    if time_rfi:
        if shortf_rfi:
            tn = not_rfi_num[0]
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
            
            s_rfi = mask_time_rfi(T,freq, rtype = 'short-freq',plot = plot,pdf = pdf,RMS = RMS,**shortf_args)
            
            t_rfi = t_rfi | s_rfi
    
        if shortf_rfi:
            T[s_rfi] = 0

    if len(is_rfi_num)>0:
        T = np.delete(T,is_rfi_num,axis = 0)
        
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
    
    if len(is_rfi_num)>0:
        T_ = np.full(T_ori.shape,np.nan)
        T_[not_rfi_num,:] = T
        T = deepcopy(T_)
        
    ###################### freq period RFI ##########################    
    is_rfio = f['is_rfi'][()] if 'is_rfi' in f.keys() else None 
    
    if period_rfi:
        #mark RFI 
        from markRFI import find_RFI
        rfi_thr = args.rfi_thr
        mw_frange = args.mw_frange

        if mw_frange is None:
            protect_use = np.zeros_like(freq,dtype='bool')
        else:
            protect_use = (freq>mw_frange[0])&(freq<mw_frange[1])
        
        if plot:
            ylim = args.ylim

            tn = not_rfi_num[0]
            spec = deepcopy(T[tn,:])
            RMS = real_rms(spec,freq,rms_sigma,rms_frange)
            is_rfi_mw = (spec > RMS * rfi_thr)
            is_rfi_ = is_rfi_mw & (~protect_use)
            _,theory = find_RFI(
                spec,freq,is_rfi_,is_rfi_mw,freq_step = freq_step,RMS = RMS,freq_thr = freq_thr,
                ext_edge = ext_edge,rfi_fit_use = 'two groups',plot = plot,pdf=pdf,ylim=ylim,
                mask_RFI_method = mask_RFI_method, small_rfi_times = small_rfi_times,chan_step = chan_step,
                mask_all_theory = mask_all_theory,freq_from_theory = freq_from_theory,mask_thr = mask_thr,
                **find_args)
        ###    
        pd_rfi = np.full(T.shape[:2], False, dtype=bool)
        if save_rfi_list:
            rfi_freq_list = np.full((T.shape[0],len(theory)+50),np.nan)

        log.info("Looking for period RFI...")
        for tn in tqdm(range(T.shape[0])):
            if tn in not_rfi_num:
                spec = deepcopy(T[tn,:]) 
                #RMS = rms(spec,freq,rms_frange)
                RMS = real_rms(spec,freq,rms_sigma,rms_frange)
                is_rfi_mw = (spec > RMS * rfi_thr)
                is_rfi = (spec > RMS * rfi_thr) & (~protect_use)
                try:
                    pd_rfi[tn,:],rfi_theory = find_RFI(
                        spec,freq,is_rfi,is_rfi_mw,freq_step = freq_step,RMS = RMS,
                        freq_thr = freq_thr,ext_edge = ext_edge,rfi_fit_use = rfi_groups ,plot = False,
                        mask_RFI_method = mask_RFI_method, small_rfi_times = small_rfi_times,chan_step = chan_step,
                        mask_all_theory = mask_all_theory,freq_from_theory = freq_from_theory,mask_thr = mask_thr,
                        **find_args)

                    if save_rfi_list:
                        if len(rfi_theory) < rfi_freq_list.shape[1]:
                            rfi_theory = np.pad(rfi_theory,(0,rfi_freq_list.shape[1] - len(rfi_theory)))
                        elif len(rfi_theory) > rfi_freq_list.shape[1]:
                            rfi_theory = rfi_theory[:rfi_freq_list.shape[1]]

                        rfi_freq_list[tn,:] = rfi_theory

                except ValueError:
                    pd_rfi[tn,:] = True
                    log.warning(f"tn={tn} has a ValueError !")
                    import traceback
                    traceback.print_exc()    


        log.info("Finish finding period RFI...")
    
    if period_rfi & time_rfi:
        rfi_mask = pd_rfi | t_rfi
    else:
        if period_rfi:
            rfi_mask = deepcopy(pd_rfi)
        if time_rfi:
            rfi_mask = deepcopy(t_rfi)
            
    if is_rfio is not None:
        if period_rfi | time_rfi:
            rfi_mask = rfi_mask | is_rfio
            log.warning("'is_rfi' alredy exists, found again and use union set.")
        else:
            rfi_mask = is_rfio
            print("Use 'is_rfi' in itself.")

    
    if np.sum(rfi_mask) > 0:

        time_coherent_per = args.time_coherent_per
        if (time_coherent_per > 0)&(time_coherent_per <1):
            per_use = (np.sum(rfi_mask,axis = 0)/rfi_mask.shape[0] > time_coherent_per)
            rfi_mask[:,per_use] = True

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

            from util import plot_waterfall
            if keep_polar:
                plot_waterfall(f,data = T_ret[:,:,0], vmin_max=[-.05,.05],cmap='plasma',figsize=(18,5),
                               title = os.path.basename(fileout).split('.')[:-1][0] + '_polar xx',pdf = pdf)

                plot_waterfall(f,data = T_ret[:,:,1], vmin_max=[-.05,.05],cmap='plasma',figsize=(18,5),
                               title = os.path.basename(fileout).split('.')[:-1][0] + '_polar yy',pdf = pdf)
            else:
                plot_waterfall(f,data = T_ret, vmin_max=[-.05,.05],cmap='plasma',figsize=(18,5),
                               title = os.path.basename(fileout).split('.')[:-1][0] + '_polar merged',pdf = pdf)              

            pdf.close()
            log.info(f"Plot to {pdfname}")
    else:
        rfi_mask = np.zeros_like(T,dtype='bool')
        log.info("Nothing is masked.")
        if plot:
            pdf.close()

    if flux:
        from hifast.core.flux import cali_src
        nB = int(re.findall(r'-M[0-1][0-9]',fname)[-1][2:])
        T_ret = cali_src(T_ret, nB, freq, cali_fname, ra=ra, dec=dec, mjd=mjd)
    
    # out dict
    dict_out={}
    dict_out['is_rfi'] = rfi_mask
    print("saved 'is_rfi'")
    if shortf_rfi:
        save_sf = args.save_sf
        if save_sf and (np.sum(s_rfi) > 0):
            dict_out['short_rfi'] = s_rfi
            print("saved 'short_rfi'")
    dict_out['freq'] = freq
    if 'ra' in f.keys():
        dict_out['ra'] = ra
        dict_out['dec'] = dec
        
    dict_out['mjd'] = mjd
    dict_out[outfield] = T_ret.astype('float32')
    if period_rfi:
        if save_rfi_list: dict_out['rfi_list'] = rfi_freq_list
    #save file
    print('Saving...')
    
    from hifast.utils.io import rec_his, save_dict_hdf5,add_extra
    add_extra(f, dict_out)
    f.close()
    
    header=rec_his(args=args)

    if header_in is not None: header.update(header_in)
    save_dict_hdf5(fileout, dict_out, header=header)
    log.info(f"Saved to {fileout}")
    
