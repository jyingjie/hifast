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
    import argparse
    parser = argparse.ArgumentParser(allow_abbrev=False)
    parser.add_argument('fname',
                        help='file name to sub standing waves')
    parser.add_argument('-rfi', '--rfi_fname',
                       help='mask rfi file name')
    parser.add_argument('-sep', '--sep_fname',
                       help='sep file name')
    parser.add_argument('-f', '--force', action='store_true',
                        help='overwriting file if out file exists')
    parser.add_argument('--frange', type=float, nargs=2,
                       help='freq range')
    parser.add_argument('--nB_radec', type=int, default=1, 
                       help='Beam number in the radec file name')
    parser.add_argument('--outdir',
                       help='default is same with the input file')
    parser.add_argument('--flux', action='store_true',
                        help='convert temperature to flux before subtract baseline')
    parser.add_argument('-c', '--cali_fname',
                       help='quasar calibration file name')
    
    #smooth
    parser.add_argument('--s_method_freq', default='None', choices=['gaussian', 'boxcar', 'median', 'None'],
                       help='smooth spec along freq')
    parser.add_argument('--s_sigma_freq', type=int, default=3,
                       help='smooth spec along freq sigma')
    parser.add_argument('--s_method_t', default='None', choices=['gaussian', 'boxcar', 'median', 'None'],
                       help='smooth spec along time')
    parser.add_argument('--s_sigma_t', type=int, default=10,
                       help='smooth spec along time sigma')
    
    ## replace big RFI
    parser.add_argument('--rfi_method', default='near ripple', choices=['subtract trpdr',
                      'subtract tr','subtract rfi', 'lower','set zeros','set noise','near ripple'],
                        help='method to replace big RFI')
    parser.add_argument('--mw_frange', type=float, nargs=2,
                       help='milky way freq range')  
    parser.add_argument('--rms_sigma', type=float, default = 6,
                       help='gauss filter sigma to compute real rms')
    parser.add_argument('--rms_frange', type=float, nargs=2,
                       help='freq range to compute rms')
    # subtract methods use sg filter and gauss smooth
    parser.add_argument('--sg_window', type=float, default=1.0, 
                       help='savgol_filter window_length (MHz)')
    parser.add_argument('--sg_polyorder', type=int, default=7, 
                       help='savgol_filter polyorder')
    parser.add_argument('--gauss_sigma',type=float, default=1.0,
                       help='gauss smooth sigma')
    # other methods need
    parser.add_argument('--times_lower_thr', type=float, default=2, 
                       help='above * times of rms will be lowered')
    parser.add_argument('--times_lower', type=float, default=1.0e4, 
                       help='lower mw area when fft, default mw/1e4')
    parser.add_argument('--ext_freq',type=float, 
                        help='extend freq range to replace (mhz)')
    parser.add_argument('--rfi_width_lim', type=float, 
                       help='rfi should contain more channels than limit')
    parser.add_argument('--ext_sec', type=int,
                       help='extend channel number of start and end of each section')
        
    ## fft remove ripple
    parser.add_argument('--fft_method', default='rfft', choices=['rfft'],
                       help='method to remove ripples')
    # fft method
    parser.add_argument('--amp_thr_mean_factor', type=float, default=1.05, 
                       help='above mean amptitude threshold will be chosed')
    parser.add_argument('--amp_thr_factor', type=float, default=1.4, 
                       help='above amptitude threshold will be chosed in every spec')
    parser.add_argument('--chan_wide', type=int, default=5, 
                       help='channel numbers near 1mhz to be chosed (wide)')
    parser.add_argument('--chan_narr', type=int, default=3, 
                       help='channel numbers near 1mhz to be chosed (narrow)')
    parser.add_argument('--choose_method', default='all', choices=['all','interpolate'],
                       help='method to choose components in fft')

    # remove which component in fft ?
    parser.add_argument('--rip_base', action='store_true',
                       help='remove constant components')
    parser.add_argument('--rip_1mhz', action='store_true',
                       help='remove 1.08mhz ripple')
    parser.add_argument('--rip_2mhz', action='store_true',
                       help='remove 1.92mhz ripple')
    parser.add_argument('--rip_0_04mhz', action='store_true',
                       help='remove 0.039 mhz ripple')
    # DELETE
    #parser.add_argument('--rfi_8mhz', action='store_true',
    #                   help='remove 8.1 mhz components')
    #parser.add_argument('--rfi_8mhz_step', type=float, 
    #                   help='big RFI linspace step, freq(\mu s) in Fourier space, nearly 1/16')
    parser.add_argument('--fft_ylim', type=float, nargs=2,
                        help='set ylim in plotting fft components')
    
    parser.add_argument('--ylim', type=float, nargs=2,
                        help='set ylim in plotting spec')
    parser.add_argument('--vmin_max', type=float, nargs=2,
                        help='vmin vmax in plotting waterfall')
    parser.add_argument('--one_spec', action='store_true',
                       help='plot only one spec or mean specs')
    
    parser.add_argument('-T', '--trans', action='store_true',
                       help='trans')
    parser.add_argument('--keep_polar', action='store_true',
                       help='keep two polarizations')
    parser.add_argument('--fill_rfi', default='rfi', choices=['nan','rfi'],
                       help='keep rfi')
    parser.add_argument('--plot', action= 'store_true',
                       help='plot')
    parser.add_argument('--no_radec', action='store_true',
                       help='do not check radec')
    
    
    args = parser.parse_args()
    file_spec = args.fname
    outdir = args.outdir
    trans= args.trans
    flux = args.flux
    
    fill_rfi = args.fill_rfi
    keep_polar = args.keep_polar
    plot = args.plot
    
    frange = args.frange
    nB_radec = args.nB_radec
    cali_fname = args.cali_fname
    rfi_fname = args.rfi_fname
    #save_rfi_fname =  args.save_rfi_fname
    
    fft_method = args.fft_method
    rfi_method = args.rfi_method
    
    # replace args
    rep_args = {}
    if rfi_method in ['subtract trpdr','subtract tr','subtract rfi',
                       'lower','set zeros','set noise','near ripple']:
        rep_args['rms_sigma'] = args.rms_sigma
        rep_args['rms_frange'] = args.rms_frange
        if 'subtract' in rfi_method:
            rep_args['sg_window'] = args.sg_window
            rep_args['sg_polyorder'] = args.sg_polyorder
            rep_args['times_lower_thr'] = args.times_lower_thr
            if 'subtract tr' in rfi_method:
                rep_args['s_sigma'] = args.gauss_sigma
            if rfi_method == 'subtract tr':
                rep_args['times_low_thr'] = args.times_lower_thr
            elif rfi_method == 'subtract rfi':
                rep_args['times_lower'] = args.times_lower
                
        else:
            rep_args['times_lower_thr'] = args.times_lower_thr
            if rfi_method == 'lower':
                rep_args['times_lower'] = args.times_lower
            elif rfi_method == 'near ripple':
                rep_args['rfi_width_lim'] = args.rfi_width_lim
                rep_args['ext_sec'] = args.ext_sec
                rep_args['ext_freq'] = args.ext_freq
            if rfi_method == 'lower' or 'set' in rfi_method :
                rfi_fname = 'none'
            
        print("rep_args:",rep_args)    
    else:
        raise ValueError("Unsupport replace RFI method.")

    #fit args
    fit_args = {}
    cw = args.chan_wide
    cn = args.chan_narr
    choose_method = args.choose_method
    if fft_method =='rfft':
        fit_args['chan_wide'] = cw
        fit_args['chan_narr'] = cn
        fit_args['amp_thr_mean_factor'] = args.amp_thr_mean_factor
        fit_args['amp_thr_factor'] = args.amp_thr_factor
        fit_args['choose_method'] = choose_method
        fit_args['rip_base'] = args.rip_base
        fit_args['rip_1mhz'] = args.rip_1mhz
        fit_args['rip_2mhz'] = args.rip_2mhz
        fit_args['rip_0_04mhz'] = args.rip_0_04mhz
        fit_args['fft_ylim'] = args.fft_ylim
        #if args.rfi_8mhz_step is not None:
        #    fit_args['rfi_8mhz'] = args.rfi_8mhz
        #    fit_args['rfi_8mhz_step'] = args.rfi_8mhz_step
        print("fit_args:",fit_args) 
    else:
        raise ValueError("Unsupport fit baseline ripple method.")
    
    fpart = '-flux_' if args.flux else '-'
    if fft_method == 'rfft': 
        fpart += 'fft_bld'
        
    if rfi_method =='subtract trpdr':
        fpart += 'p'
    elif rfi_method =='subtract tr':
        fpart += 't'
    elif rfi_method =='subtract rfi':
        fpart += 'r'
    elif rfi_method =='lower':
        fpart += 'l'
    elif rfi_method =='set zeros':
        fpart += 'z'
    elif rfi_method =='set noise':
        fpart += 'n'
    elif rfi_method =='near ripple':
        fpart += 'e'
    
    if choose_method == 'all':
        fpart += 'a'
    elif choose_method == 'interpolate':
        fpart += 'i'
    # test
    #fpart += f'{cw}_{cn}'    
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
 
    nB = int(re.findall(r'-M[0-1][0-9]',file_spec)[-1][2:])

    ################################ load data ###############################
    fs = h5py.File(file_spec,'r')
    mjd = fs['mjd'][()]
    
    no_radec = args.no_radec
    if not no_radec:
        if 'ra' not in fs.keys():
            from hifast.cli_baseline import get_radec 
            ra, dec, is_extrapo = get_radec(file_spec, nB, nB_radec, mjd)
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
            
    if 'is_on' in fs.keys():
        is_on = fs['is_on'][()]
    else:
        raise(ValueError(f'Can not find is_on'))

    ind_sort = np.argsort(mjd)
    mjd = mjd[ind_sort]
    ra = ra[ind_sort]
    dec = dec[ind_sort]
    T = T[ind_sort]
    is_on = is_on[ind_sort]
    
    # read RFI
    from glob import glob
    
    if rfi_fname is not None:
        if rfi_fname == 'none':
            print("Don't use rfi file.")
    else:
        rfi_dir_default = os.path.dirname(file_spec).split('/')[:-1]
        rfi_dir_default.append('rfi')
        
        rfi_dir_default = '/'.join(rfi_dir_default)
        if rfi_method =='subtract trpdr' or rfi_method =='subtract rfi':
            outpart = '*-tr_pdr.hdf5'
        if rfi_method =='subtract tr' or rfi_method =='near ripple':
            outpart = '*-tr*.hdf5'
        rfi_dir = os.path.join(rfi_dir_default, '.'.join(os.path.basename(file_spec).split('.')[:-1]) + f'{outpart}')

        rfi_fnames = glob(rfi_dir)
        if len(rfi_fnames) == 1:
            rfi_fname = rfi_fnames[0]
        elif len(rfi_fnames) == 0:
            rfi_dir = os.path.join(os.path.dirname(file_spec),'.'.join(os.path.basename(file_spec).split('.')[:-1]) + f'{outpart}')
            rfi_fnames = glob(rfi_dir)
            if len(rfi_fnames) < 3:
                rfi_fname = rfi_fnames[0]
            else:
                raise FileNotFoundError("Where is rfi mask?")
        else:
            raise FileNotFoundError("Which rfi mask do you want? ")
    
    if rfi_fname != 'none':
        rfi = h5py.File(rfi_fname,'r')
        print(f"Find rfi file {rfi_fname}")
        
        if 'is_rfi' in rfi.keys():
            is_rfi = rfi['is_rfi'][()]
            if 'time_rfi' in rfi.keys():
                t_rfi = rfi['time_rfi'][()]
            else:
                if rfi_method =='subtract tr' or rfi_method =='near ripple':
                    t_rfi = is_rfi
                else:
                    t_rfi = np.full(is_rfi.shape,False)
                    log.warning("rfi file doesn't contain key 'time_rfi',so don't use short time rfi.")
            rfi.close()
            if len(is_rfi.shape) == 3:
                # use 2D rfi mask
                is_rfi = is_rfi[:,:,0]|is_rfi[:,:,1]
        else:
            if 'T' in rfi.keys():
                Tr = rfi['T'][()]
            elif 'Ta' in rfi.keys():
                Tr = rfi['Ta'][()]
            elif 'flux' in rfi.keys():
                Tr = rfi['flux'][()]
            else:
                print(f"{rfi.keys()}")
                raise ValueError(f"rfi keys don't have 'T','Ta' or 'flux'.")
            if len(Tr.shape) == 3:
                Tr = np.mean(Tr, axis=2, dtype='float64')
            is_rfi = np.isnan(Tr)
        
        whole_rfi = np.all(is_rfi,axis = 1)
        not_rfi_num = np.arange(T.shape[0])[~whole_rfi]
        is_rfi_num = np.arange(T.shape[0])[whole_rfi]
    else:
        is_rfi = None;t_rfi = None
        not_rfi_num = np.arange(T.shape[0])
        is_rfi_num = np.array([])

    freq = fs['freq'][:]
    fdelta = freq[1] - freq[0]
    if frange is not None:
        is_ = (freq >= frange[0]) & (freq <= frange[1])
        freq = freq[is_]
        T = T[:,is_,:]
        if rfi_fname != 'none':
            is_rfi = is_rfi[:,is_];t_rfi = is_rfi[:,is_]
    
    # load data
    sep_fname = args.sep_fname    
    if sep_fname is None:
        # find higher class folder
        sep_dir_default = os.path.dirname(file_spec).split('/')[:-1]
        sep_dir_default.append('sep')
        sep_dir_default = '/'.join(sep_dir_default)
        sep_dir = os.path.join(sep_dir_default, '.'.join(os.path.basename(file_spec).split('-bld',1)[:-1]) + '.hdf5')
        
        sep_fnames = glob(sep_dir)
        if len(sep_fnames) == 1:
            sep_fname = sep_fnames[0]
        elif len(sep_fnames) == 0:
            #find in same class folder
            sep_fname = os.path.join(os.path.dirname(file_spec),'.'.join(os.path.basename(file_spec).split('-bld',1)[:-1]) +'.hdf5')
        else:
            raise FileNotFoundError("Which sep file do you want? ")
    else:
        sep_fnames = glob(sep_fname)
        if len(sep_fnames) == 1:
            sep_fname = sep_fnames[0]
        else:
            filedir = os.path.dirname(sep_fname)
            sep_fname = os.path.join(filedir, os.path.basename(file_spec).split('-bld')[0] + '.hdf5')

    sep = h5py.File(sep_fname,'r')
    print(f"Find sep file {sep_fname}") 
    from util import get_data
    T_sep = get_data(sep,polar = 'merged',xrange = frange)

    if T_sep.shape != T.shape:
        raise ValueError(f"sep data shape {T_sep.shape} is not matched with sub data shape {T.shape}.")
    
    mw_frange = args.mw_frange
    if mw_frange is None:
        if (max(freq) <= 1419)|(min (freq)>= 1422):
            print("don't contain MW")
            mw_use = np.zeros_like(freq,dtype='bool')
        else:
            if rfi_method == 'subtract':
                raise ValueError("--mw_range is None!")
            else:
                mw_use = None
    else:
        mw_use = (freq>=mw_frange[0])&(freq<=mw_frange[1])
    
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

    if trans:
        T = T.transpose((1,0,2))
        if flux:
            raise()     
  
    if flux:
        ra = comm.scatter(ra, root=0)
        dec = comm.scatter(dec, root=0)
        from hifase.flux import cali_src
        print('Flux calibrating ...')
        T = cali_src(T, nB, freq, cali_fname, ra=ra, dec=dec, mjd=mjd)
    
    ########################### smooth ###############################3
    s_method_freq = args.s_method_freq
    s_sigma_freq = args.s_sigma_freq
    s_method_t = args.s_method_t
    s_sigma_t = args.s_sigma_t
    def do_smooth(T,s_method_t,s_sigma_t,s_method_freq,s_sigma_freq):
        from hifast.utils.misc import smooth1d
        if  s_method_t in ['gaussian', 'boxcar', 'median']:
            print('Smooth ing ...')
            T = smooth1d(T,axis = 0,sigma = s_sigma_t, method = s_method_t)

        if  s_method_freq in ['gaussian', 'boxcar', 'median']:
            print('Smooth ing ...')
            T = smooth1d(T,axis = 1,sigma = s_sigma_freq, method = s_method_freq)
        return T
    
    ######################### fft ##############################
    
    from baseline_2 import replace_rfi,fit_ripple
    if len(T.shape) == 3:
        if keep_polar:
            T3 = deepcopy(T)
            T_xx = T[:,:,0]
            T_yy = T[:,:,1]
            ori_shape = T_xx.shape
            print("polar xx ...")
            # replace big rfi
            data_rep_xx = replace_rfi(T_xx,freq,is_rfi,time_rfi=t_rfi,
                                               method = rfi_method,mw_use =mw_use,**rep_args)
            data_rep_xx = do_smooth(data_rep_xx,
                                             s_method_t,s_sigma_t,s_method_freq,s_sigma_freq)
            # get standing waves
            sw_fit_xx = fit_ripple(data_rep_xx, freq,fft_method,is_rfi_num,not_rfi_num,ori_shape,
                                   is_on = is_on, plot = plot,pdf=pdf,title='polar xx',**fit_args)
            print("polar yy ...")
            data_rep_yy = replace_rfi(T_yy,freq,is_rfi,time_rfi=t_rfi,
                                               method = rfi_method,mw_use =mw_use,**rep_args)
            data_rep_yy = do_smooth(data_rep_yy,
                                             s_method_t,s_sigma_t,s_method_freq,s_sigma_freq)
            
            sw_fit_yy = fit_ripple(data_rep_yy, freq,fft_method,is_rfi_num,not_rfi_num,ori_shape,
                                   is_on = is_on, plot = plot,pdf=pdf,title='polar yy',**fit_args)
            
        else:
            T = np.mean(T, axis=2, dtype='float64')   
            ori_shape = T.shape
            
            # replace big rfi
            data_rep, is_rfi_no_mw = replace_rfi(T,freq,is_rfi,time_rfi=t_rfi,
                                                          method = rfi_method,mw_use =mw_use,**rep_args)
            # get standing waves
            sw_fit = fit_ripple(data_rep, freq,fft_method,is_rfi_num,not_rfi_num,ori_shape,
                                   plot = plot,pdf=pdf,title='polar merged',**fit_args)
    else:
        raise ValueError('data should be 3D, and has 2 polars') 

    
    #sub remove standing waves 
    if keep_polar: 
        sw_fit = np.append(sw_fit_xx,sw_fit_yy)
        sw_fit = sw_fit.reshape((2,ori_shape[0],ori_shape[1])).transpose((1,2,0))
        rmsw_data = T3 - sw_fit
    else:
        rmsw_data = T_ori - sw_fit
    
    # fill rfi with ?
    if fill_rfi == 'nan':
        if rfi_fname != 'none':
            rmsw_data[is_rfi] = np.nan
    elif fill_rfi == 'rfi':
        pass
    
           
    if plot:
        ylim = args.ylim
        vmin_max = args.vmin_max
        one_spec = args.one_spec
        print(" 'Wait for plotting patiently, you must.' Master Yoda said.")
        tn = not_rfi_num[10]
        def plot_in_pdf(data_rep,T,sw_fit,rmsw_data,polar,pdf = None,one_spec = False,
                        ylim = None,vmin_max=None):
            global tn, freq
            fig = plt.figure(figsize=(40,4))
            ax = fig.add_subplot(111)
            if not one_spec:
                ax.hlines([0,-.5],freq[0],freq[-1],alpha = .8)
                ax.plot(freq,np.mean(data_rep[tn-5:tn+5,:],axis = 0),'b',label='rm rfi',alpha = .5)
                ax.plot(freq,np.mean(T[tn-5:tn+5,:],axis = 0),label='original')
                ax.plot(freq,np.mean(sw_fit[tn-5:tn+5,:],axis = 0),label='ripple')
                ax.plot(freq,np.mean(rmsw_data[tn-5:tn+5,:],axis = 0) - .5,label='result')
                ax.set_title(f'ten specs mean, polar {polar}')
                if ylim is not None:
                    ax.set_ylim(ylim[0],ylim[1])
            else:
                ax.hlines([0,-2],freq[0],freq[-1],alpha = .8)
                ax.plot(freq,data_rep[tn,:],'b',label='rm rfi',alpha = .5)
                ax.plot(freq,T[tn,:],label='original')
                ax.plot(freq,sw_fit[tn,:],label='ripple')
                ax.plot(freq,rmsw_data[tn,:] - 2,label='result')
                ax.set_title(f'single spec, polar {polar}')
                if ylim is not None:
                    ax.set_ylim(ylim[0],ylim[1])
            ax.grid();ax.legend();
            ax.set_xlim(freq[0],freq[-1])
            pdf.savefig();plt.close()

            from util import plot_waterfall
            plot_waterfall(fs,data = T, vmin_max=vmin_max,cmap='plasma',figsize=(18,5),
                           title = os.path.basename(file_spec).split('.')[:-1][0],pdf = pdf)

            plot_waterfall(fs,data = rmsw_data, vmin_max=vmin_max,cmap='plasma',figsize=(18,5),pdf = pdf,
                           title = f'remove standing waves, polar {polar}')
            
        if keep_polar:
            plot_in_pdf(data_rep_xx,T_xx,sw_fit_xx,rmsw_data[:,:,0],polar='xx',pdf = pdf,
                        one_spec = one_spec, ylim = ylim,vmin_max=vmin_max)
            plot_in_pdf(data_rep_yy,T_yy,sw_fit_yy,rmsw_data[:,:,1],polar='yy',pdf = pdf,
                       one_spec = one_spec, ylim = ylim,vmin_max=vmin_max)
        else:
            plot_in_pdf(data_rep,T,sw_fit,rmsw_data,polar='merged',pdf = pdf)
        
        pdf.close()
        log.info(f"Plot to {pdfname}")
        
        
    #sep remove standing waves 
    if keep_polar: 
        rmsw_data = T_sep - sw_fit
    else:
        rmsw_data = np.mean(T_sep,axis = 2) - sw_fit
        
    # fill rfi with ?
    if fill_rfi == 'nan':
        if rfi_fname != 'none':
            rmsw_data[is_srfi] = np.nan
    elif fill_rfi == 'rfi':
        pass
    
    print(f"Saving...")
    dict_out= {}
    dict_out['mjd'] = mjd
    dict_out['ra'] = ra
    dict_out['dec'] = dec
    dict_out[outfield] = rmsw_data.astype('float32')
    #dict_out['ripple'] = sw_fit.astype('float32')
    dict_out['freq'] = freq
    if is_rfi is not None:
        dict_out['is_rfi'] = is_rfi
    if rfi_fname != 'none':      
        rfi.close()
    if is_extrapo is not None:
        dict_out['is_extrapo'] = is_extrapo
    #save file
    from hifast.utils.io import rec_his, save_dict_hdf5,add_extra
    add_extra(fs, dict_out)
    fs.close()
    header=rec_his(args=args);
    if header_in is not None: header.update(header_in)
    save_dict_hdf5(fileout, dict_out, header=header)
    log.info(f"Saved to {fileout}")
        
