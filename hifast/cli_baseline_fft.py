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
from fft import FFT
#from .util import extend_Trues#, boxcar_smooth1d, median_filter_1d
 
    
def replace_rfi(data, is_rfi, method,**rep_args):
    """
    replace big RFI to reduce the impact in FFT
    
    Parameters:
    data: np array
    sg_window: scipy.signal.savgol_filter window length (unit MHz)
    sg_polyorder: scipy.signal.savgol_filter polyorder
    mw_lower: lower milky way area by divide a number
    """
    global mw_use,fdelta,is_rfi_num,not_rfi_num
    
    sg_window = rep_args['sg_window']
    sg_polyorder = rep_args['sg_polyorder']
    mw_lower = rep_args['mw_lower']
    times_lower = rep_args['times_lower']

    if method == 'subtract':

        is_rfi_no_mw = deepcopy(is_rfi)
        is_rfi_no_mw[:,mw_use] = False
        is_rfi_no_mw[is_rfi_num,:] = True
        data_no_mw = deepcopy(data)
        data_no_mw[:,mw_use] = 0
        data_rmrfi = deepcopy(data)
        data_rmrfi[is_rfi_num,:] = np.nan
        
        from .util import _round_up_to_odd_integer
        window_length = _round_up_to_odd_integer(sg_window/fdelta)
        #print(window_length,sg_polyorder,data_no_mw.shape)
        
        log.info(f"Replace RFI with {method} method ...")
        from scipy.signal import savgol_filter
        for tn in tqdm(range(data.shape[0])):
            if tn in not_rfi_num:
                sg = savgol_filter(data_no_mw[tn,:],window_length = window_length ,polyorder=sg_polyorder,)
                data_rmrfi[tn,is_rfi_no_mw[tn]] = data[tn,is_rfi_no_mw[tn]] - sg[is_rfi_no_mw[tn]]

        from .markRFI import rms
        RMS = rms(data_rmrfi[0,:],freq,[freq[0],freq[0]+20])
        thr_lower = RMS*times_lower
        low_use = (data_rmrfi > thr_lower) 
        from .util import extend_Trues
        low_use = extend_Trues(low_use,ext_add =10,leng_lim = 20,axis = -1)
          
        data_rmrfi_low_mw = deepcopy(data_rmrfi)
        data_rmrfi_low_mw[low_use] = data[low_use]/mw_lower
        data_rmrfi_low_mw[is_rfi_num,:] = np.nan      

    return data_rmrfi, data_rmrfi_low_mw, is_rfi_no_mw
    

def fit_ripple(data_rmrfi_low_mw, method, plot = False,**fit_args): 
    """
    fit baseline ripple (standing wave) by FFT
    Parameter:
    data_rmrfi_low_mw: data array (after lower mw)
    sw_freq: standing wave 'frequence' in Fourier space, unit \mu s
    rfi_freq_step: big RFI residual influence in Fourier space, nearly 1/16.2 \mu s.
                    (period is 16.2 MHz)
    amp_thr: above amptitude threshold will be chosed.
    sw_n: channel number next to sw_freq will be chosed all.
    
    """
    
    global freq,fdelta,is_rfi_num,not_rfi_num,ori_shape
    sw_freq = fit_args['sw_freq']
    rfi_freq_step = fit_args['rfi_freq_step'] 
    amp_thr = fit_args['amp_thr'] 
    sw_n = fit_args['sw_n'] 

    if method == 'rfft':
        log.info(f"Fit baseline ripple with {method} method ...")
        
        if len(is_rfi_num) > 0:
            data_rmrfi_low_mw = np.delete(data_rmrfi_low_mw,is_rfi_num,axis = 0)
        
        fftf = FFT(data_rmrfi_low_mw, freq)
        x = fftf.x
        amp_data = fftf.amp
        loc1 = np.argmin(np.abs((x - sw_freq)))
        x_loc16 = np.arange(0,x[-1],rfi_freq_step)

        N = freq.shape[0]
        
        fs = 1/fdelta
        loc16 = np.around(x_loc16 / (fs/N)).astype('int')
        # two sides
        loc16_ = np.hstack((loc16 - 1,loc16,loc16 + 1))
        loc16_.sort()

        use16 = np.full(len(x),False)
        use16[loc16_] = True;use16[loc1] = False ; use16[0] = False
        use16 = use16 &  (amp_data >= amp_thr) 

        use = (np.abs(np.arange(len(x))-loc1) < sw_n) & (amp_data >= amp_thr) 

        from scipy.interpolate import interp1d
        x_ = deepcopy(x)
        # interplate
        amp_data_inpd = deepcopy(amp_data)
        for ti in tqdm(range(amp_data.shape[0])):
            x_mask = x_[~use16[ti]]
            amp_data_mask = amp_data[ti][~use16[ti]]
            amp_data_interp = interp1d(x_mask,amp_data_mask,kind='linear')#,fill_value="extrapolate")
            amp_data_inpd[ti][use16[ti]] = amp_data_interp(x_[use16[ti]])
        amp_data_inpd = amp_data - amp_data_inpd
        amp_data_inpd[:,0] = amp_data[:,0]
        amp_data_inpd[use] = amp_data[use]
        amp_data_inpd[amp_data_inpd < 0] = 0
        
        if plot:
            tn = not_rfi_num[0]
    
            is_ = (x >= 0) & (x <= 3)
            fig = plt.figure(figsize=(22,5))
            ax = fig.add_subplot(111) 
            ax.stem(x[is_],amp_data[tn,is_],linefmt='--',markerfmt ='C0o',label = 'whole fft')
            ax.stem(x[is_],amp_data_inpd[tn,is_],linefmt='--',markerfmt ='C1o',label = 'sw fft')
            ax.text(x[loc1],.1,f'{loc1}')
            ax.plot(x[use[tn]],np.zeros_like(x[use[tn]]),'co',label = f'near {x[loc1]:.4f} $\mu$s ')
            ax.plot(x[loc16],np.zeros_like(loc16),'gs',label = f'step {rfi_freq_step:.4f} $\mu$s ')
            ax.axhline(amp_thr,color = 'k',linestyle='--')
            ax.grid();ax.set_xlim(0,3)
            ax.set_xlabel('k [$\mu$s]')
            plt.legend()
            plt.tight_layout()

        A_data_inpd =  amp_data_inpd * np.exp(1j*fftf.phi)
        A_data_ifft = np.real(np.fft.irfft(A_data_inpd,n=data_rmrfi_low_mw.shape[1]))
        
        if len(is_rfi_num)>0:
            A_data_ifft_ = np.full(ori_shape,np.nan)
            A_data_ifft_[not_rfi_num,:] = A_data_ifft
            A_data_ifft = deepcopy(A_data_ifft_)
        
        return A_data_ifft


if __name__ == '__main__':
    import warnings 
    warnings.filterwarnings("ignore",r'overflow encountered in exp')
    import argparse
    parser = argparse.ArgumentParser(allow_abbrev=False)
    parser.add_argument('fname',
                        help='file name to sub standing waves')
    parser.add_argument('-rfi', '--rfi_fname',
                       help='mask rfi file name')
    #parser.add_argument('--nproc', type=int,
    #                    help='number process')
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
    
    
    parser.add_argument('--rfi_method', default='subtract', choices=['subtract'],
                       help='method to replace big RFI')
    parser.add_argument('--mw_frange', type=float, nargs=2,
                       help='milky way freq range')
    parser.add_argument('--protect_mw', action= 'store_true',
                       help='protect mw')    
    parser.add_argument('--sg_window', type=float, default=1.0, 
                       help='savgol_filter window_length (MHz)')
    parser.add_argument('--sg_polyorder', type=int, default=7, 
                       help='savgol_filter polyorder')
    parser.add_argument('--mw_lower', type=float, default=1.0e4, 
                       help='lower mw area when fft, default mw/1e4')
    parser.add_argument('--times_lower', type=float, default=5, 
                       help='above 5 times of rms will be lowered')
    
    
    parser.add_argument('--fft_method', default='rfft', choices=['rfft'],
                       help='method to remove ripples')
    parser.add_argument('--sw_freq', type=float, default=0.9254, 
                       help='standing waves freq(\mu s) in Fourier space')
    parser.add_argument('--rfi_freq_step', type=float, default=0.0617283950617284, 
                       help='big RFI linspace step, freq(\mu s) in Fourier space, nearly 1/16')
    parser.add_argument('--amp_thr', type=float, default=35, 
                       help='above amptitude threshold will be chosed')
    parser.add_argument('--sw_n', type=int, default=5, 
                       help='channel numbers near 1mhz to be chosed')
    
    
    parser.add_argument('-T', '--trans', action='store_true',
                       help='')
    #parser.add_argument('--keep_polar', action='store_true',
    #                   help='keep two polarizations')
    parser.add_argument('--fill_rfi', default='nan', choices=['nan','noise','rfi'],
                       help='keep rfi')
    parser.add_argument('--plot', action= 'store_true',
                       help='plot')
    
    
    args = parser.parse_args()
    #nproc = args.nproc
    file_spec = args.fname
    outdir = args.outdir
    trans= args.trans
    flux = args.flux
    
    fill_rfi = args.fill_rfi
    #keep_polar = args.keep_polar
    plot = args.plot
    
    fpart = '-flux_' if args.flux else '-'
    fpart += 'fft_rfi_bld'
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
        

    frange = args.frange
    nB_radec = args.nB_radec
    cali_fname = args.cali_fname
    rfi_fname = args.rfi_fname
    
    fft_method = args.fft_method
    rfi_method = args.rfi_method

    # replace args
    rep_args = {}
    if rfi_method in ['subtract']:
        rep_args['sg_window'] = args.sg_window
        rep_args['sg_polyorder'] = args.sg_polyorder
        rep_args['mw_lower'] = args.mw_lower
        rep_args['times_lower'] = args.times_lower
    else:
        raise ValueError("Unsupport replace RFI method.")
    
    
    #fit args
    fit_args = {}
    if fft_method in ['rfft']:
        fit_args['sw_freq'] = args.sw_freq
        fit_args['rfi_freq_step'] = args.rfi_freq_step
        fit_args['amp_thr'] = args.amp_thr
        fit_args['sw_n'] = args.sw_n
    else:
        raise ValueError("Unsupport fit baseline ripple method.")
        
    nB = int(re.findall(r'-M[0-1][0-9]',file_spec)[-1][2:])

    # load data
    fs = h5py.File(file_spec,'r')
    mjd = fs['mjd'][()]
    if 'ra' not in fs.keys():
        from .cli_baseline import get_radec 
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
    
    if flux:
        outfield = 'flux'

    ind_sort = np.argsort(mjd)
    mjd = mjd[ind_sort]
    ra = ra[ind_sort]
    dec = dec[ind_sort]
    T = T[ind_sort]
    
    # read RFI
    if rfi_fname is not None:
        rfi = h5py.File(rfi_fname)
        print(f"Use rfi file {rfi_fname}")
    else:
        rfi_dir_default = os.path.dirname(file_spec).split('/')[:-1]
        rfi_dir_default.append('cor_vel_once')
        rfi_dir_default = '/'.join(rfi_dir_default)
        outpart = '*-pdrfi*'
        rfi_dir = os.path.join(rfi_dir_default, '.'.join(os.path.basename(file_spec).split('.')[:-1]) + f'{outpart}')
        from glob import glob
        rfi_fnames = glob(rfi_dir)
        if len(rfi_fnames) == 1:
            print(f"Find rfi file {rfi_fnames[0]}")
            rfi = h5py.File(rfi_fnames[0],'r')
        elif len(rfi_fnames) == 0:
            raise FileNotFoundError("Run cli_multi or cli_markRFI first!")
        else:
            raise FileNotFoundError("Which rfi mask do you want? ")
    
    
    if 'is_rfi' in rfi.keys():
        is_rfi = rfi['is_rfi'][()]
        rfi.close()
        if len(is_rfi.shape) == 3:
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
    
    
    freq = fs['freq'][:]
    fdelta = freq[1] - freq[0]
    if frange is not None:
        is_ = (freq >= frange[0]) & (freq <= frange[1])
        freq = freq[is_]
        T = T[:,is_]
        is_rfi = is_rfi[:,is_]
    
    mw_frange = args.mw_frange
    if mw_frange is None:
        if (max(freq) < 1419)|(min(freq)>1422):
            print("don't contain MW")
        else:
            raise ValueError("--mw_frange is empty!")
        #mw_use = np.full(len(freq),False)
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

    #close files
    # fs.close() # don't close, load extra later

    if trans:
        T = T.transpose((1,0,2))
        if flux:
            raise()     
  
    if flux:
        ra = comm.scatter(ra, root=0)
        dec = comm.scatter(dec, root=0)
        from .flux import cali_src
        print('Flux calibrating ...')
        T = cali_src(T, nB, freq, cali_fname, ra=ra, dec=dec, mjd=mjd)
    
    #if not keep_polar:
    if len(T.shape) == 3:
        T = np.mean(T, axis=2, dtype='float64')
        
    ori_shape = T.shape
        
    # replace big rfi
    data_rmrfi , data_rmrfi_low_mw, is_rfi_no_mw = replace_rfi(T,is_rfi,method = rfi_method,**rep_args)
    
    # get standing waves
    sw_fit = fit_ripple(data_rmrfi_low_mw,method = fft_method,plot = plot, **fit_args)
    if plot:
        pdf.savefig();plt.close()
    
    # remove standing waves
    rmsw_data = T - sw_fit
    
    # protect_mw
    protect_mw = args.protect_mw
    if protect_mw:
        is_rfi = deecopy(is_rfi_no_mw)
    
    # fill rfi with ?
    if fill_rfi == 'nan':
        rmsw_data[is_rfi] = np.nan
    elif fill_rfi == 'noise':
        rmsw_data[is_rfi] = data_rmrfi[is_rfi]
    elif fill_rfi == 'rfi':
        pass
        
    if plot:
        print(" 'Wait for plotting patiently, you must.' Master Yoda said")
        tn = not_rfi_num[10]
        fig = plt.figure(figsize=(40,4))
        ax = fig.add_subplot(111)
        ax.hlines([0,-.5],1320,1440,alpha = .8)
        ax.plot(freq,np.mean(data_rmrfi_low_mw[tn-5:tn+5,:],axis = 0),'b',label='rm rfi',alpha = .5)
        ax.plot(freq,np.mean(T[tn-5:tn+5,:],axis = 0),label='original')
        ax.plot(freq,np.mean(sw_fit[tn-5:tn+5,:],axis = 0),label='ripple')
        ax.plot(freq,np.mean(rmsw_data[tn-5:tn+5,:],axis = 0) - .5,label='result')
        ax.grid();ax.legend();ax.set_title('ten specs mean')
        ax.set_xlim(1320,1440)
        ax.set_ylim(-1,.5)
        pdf.savefig();plt.close()
        
        from .util import plot_waterfall
        plot_waterfall(fs,data = T, vmin_max=[-.05,.05],cmap='plasma',figsize=(18,5),
                       title = os.path.basename(file_spec).split('.')[:-1],pdf = pdf)
        
        plot_waterfall(fs,data = rmsw_data, vmin_max=[-.05,.05],cmap='plasma',figsize=(18,5),pdf = pdf,
                       title = 'remove standing waves')
        
        pdf.close()
        log.info(f"Plot to {pdfname}")
        
        
    print(f"Saving...")
    dict_out= {}
    dict_out['mjd'] = mjd
    dict_out['ra'] = ra
    dict_out['dec'] = dec
    dict_out[outfield] = rmsw_data.astype('float32')
    dict_out['freq'] = freq
    dict_out['is_rfi'] = is_rfi
    if is_extrapo is not None:
        dict_out['is_extrapo'] = is_extrapo
    #save file
    from .util import add_extra
    add_extra(fs, dict_out)
    fs.close()
    rfi.close()
    from .util import rec_his, save_dict_hdf5
    header=rec_his(args=args)
    if header_in is not None: header.update(header_in)
    save_dict_hdf5(fileout, dict_out, header=header)
    log.info(f"Saved to {fileout}")
        
