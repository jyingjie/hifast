#!/usr/bin/env python
# coding: utf-8

# author: Xu Chen, 2021.07
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
                        help='sep file name to sub baseline')
    parser.add_argument('-rfi', '--rfi_fname',
                       help='mask rfi file name')
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
 

    parser.add_argument('--sub_method', default='mean', choices=['mean','median'],
                       help='method to remove baseline')
    parser.add_argument('--nspec',type=int, default=10,
                        help='average how many specs to fit baseline.')
    parser.add_argument('--func', default='iter', choices=['iter','smooth'],
                       help='function to mean or median')
    
    parser.add_argument('-T', '--trans', action='store_true',
                       help='trans')
    parser.add_argument('--keep_polar', action='store_true',
                       help='keep two polarizations')
    parser.add_argument('--plot', action= 'store_true',
                       help='plot')
    parser.add_argument('--no_radec', action='store_true',
                       help='do not check radec')
    
    parser.add_argument('--ylim', type=float, nargs=2,
                        help='set ylim in plotting spec')
    parser.add_argument('--vmin_max', type=float, nargs=2,
                        help='vmin vmax in plotting waterfall')
    parser.add_argument('--one_spec', action='store_true',
                       help='plot only one spec or mean specs')
    parser.add_argument('--fill_rfi', default='rfi', choices=['nan','rfi'],
                       help='keep rfi')
    
    args = parser.parse_args()
    file_spec = args.fname
    outdir = args.outdir
    trans= args.trans
    flux = args.flux

    sub_method = args.sub_method
    keep_polar = args.keep_polar

    if sub_method in ['mean','median']:
        fit_args = {}
        nspec = args.nspec
        fit_args['nspec'] = nspec
        fit_args['func'] = args.func
        print("fit args:",fit_args)
    else:
        raise ValueError("Unsupport fit baseline ripple method.")
    
    fpart = '-flux_' if args.flux else '-'
    if sub_method == 'mean': 
        fpart += f'mean{nspec}_bld'
    elif sub_method == 'median': 
        fpart += f'med{nspec}_bld'

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
            
    plot = args.plot
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
        

    frange = args.frange
    nB_radec = args.nB_radec
    cali_fname = args.cali_fname
 
    nB = int(re.findall(r'-M[0-1][0-9]',file_spec)[-1][2:])

    # load data
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

    ind_sort = np.argsort(mjd)
    mjd = mjd[ind_sort]
    ra = ra[ind_sort]
    dec = dec[ind_sort]
    T = T[ind_sort]
        
    # read RFI
    from glob import glob
    rfi_fname = args.rfi_fname
    
    if rfi_fname is not None:
        if rfi_fname == 'none':
            print("Don't use rfi file.")
    else:
        rfi_dir_default = os.path.dirname(file_spec).split('/')[:-1]
        rfi_dir_default.append('rfi')
        
        rfi_dir_default = '/'.join(rfi_dir_default)
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
    else:
        is_rfi = None 
    
    freq = fs['freq'][:]
    if frange is not None:
        is_ = (freq >= frange[0]) & (freq <= frange[1])
        freq = freq[is_]
        T = T[:,is_,:]
        is_rfi = is_rfi[:,is_]
        
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
    T_ori = deepcopy(T)
    if flux:
        ra = comm.scatter(ra, root=0)
        dec = comm.scatter(dec, root=0)
        from hifase.flux import cali_src
        print('Flux calibrating ...')
        T = cali_src(T, nB, freq, cali_fname, ra=ra, dec=dec, mjd=mjd)
    
    if rfi_fname != 'none':
        if keep_polar:
            T[is_rfi,:] = np.nan
        else:
            T[is_rfi] = np.nan
    
    from baseline_2 import fit_ripple
    if len(T.shape) == 3:
        if keep_polar:
            #T3 = deepcopy(T)
            T_xx = T[:,:,0]
            T_yy = T[:,:,1]
            ori_shape = T_xx.shape
            print("polar xx ...")
            # get standing waves
            sw_fit_xx = fit_ripple(T_xx,method = sub_method, **fit_args)
            print("polar yy ...")
            sw_fit_yy = fit_ripple(T_yy,method = sub_method,  **fit_args)
            
        else:
            T = np.mean(T, axis=2, dtype='float64')   
            sw_fit = fit_ripple(T,method = sub_method,  **fit_args)
    else:
        raise ValueError('data should be 3D, and has 2 polars') 

    #sub remove standing waves 
    if keep_polar: 
        sw_fit = np.append(sw_fit_xx,sw_fit_yy)
        sw_fit = sw_fit.reshape((2,ori_shape[0],ori_shape[1])).transpose((1,2,0))

    rmsw_data = T_ori - sw_fit
    fill_rfi = args.fill_rfi
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
        print(" 'Wait for plotting patiently, you must.' Master Yoda said")
        tn = 10
        def plot_in_pdf(T,sw_fit,rmsw_data,polar,pdf = None,one_spec = False,frange = None,
                        ylim = None,vmin_max=None):
            global tn, freq
            fig = plt.figure(figsize=(40,4))
            ax = fig.add_subplot(111)
            if not one_spec:
                ax.hlines([0,-.5],freq[0],freq[-1],alpha = .8)
                ax.plot(freq,np.mean(T[tn-5:tn+5,:],axis = 0),label='original')
                ax.plot(freq,np.mean(sw_fit[tn-5:tn+5,:],axis = 0),label='ripple')
                ax.plot(freq,np.mean(rmsw_data[tn-5:tn+5,:],axis = 0) - .5,label='result')
                ax.set_title(f'ten specs mean, polar {polar}')
                if ylim is not None:
                    ax.set_ylim(ylim[0],ylim[1])
            else:
                ax.hlines([0,-2],freq[0],freq[-1],alpha = .8)
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
            plot_waterfall(fs,data = T, vmin_max=vmin_max,cmap='plasma',figsize=(18,5),xrange = frange,
                           title = os.path.basename(file_spec).split('.')[:-1][0],pdf = pdf)

            plot_waterfall(fs,data = rmsw_data, vmin_max=vmin_max,cmap='plasma',figsize=(18,5),pdf = pdf,xrange = frange,
                           title = f'remove baseline, polar {polar}')
            
        if keep_polar:
            plot_in_pdf(T_xx,sw_fit_xx,rmsw_data[:,:,0],polar='xx',pdf = pdf,frange = frange,
                        one_spec = one_spec, ylim = ylim,vmin_max=vmin_max)
            plot_in_pdf(T_yy,sw_fit_yy,rmsw_data[:,:,1],polar='yy',pdf = pdf,frange = frange,
                       one_spec = one_spec, ylim = ylim,vmin_max=vmin_max)
        else:
            plot_in_pdf(T,sw_fit,rmsw_data,polar='merged',pdf = pdf,frange = frange)

        
        pdf.close()
        log.info(f"Plot to {pdfname}")
        

    print(f"Saving...")
    dict_out= {}
    dict_out['mjd'] = mjd
    dict_out['ra'] = ra
    dict_out['dec'] = dec
    dict_out[outfield] = rmsw_data.astype('float32')
    #dict_out['ripple'] = sw_fit.astype('float32')
    dict_out['freq'] = freq
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
        
