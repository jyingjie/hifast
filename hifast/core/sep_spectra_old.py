#!/usr/bin/env python
# coding: utf-8

# In[1]:

import numpy as np
import h5py
from .tcal_onoff_old import tcal_onoff
from ..utils.tcal import read_tcal
from ..utils.misc import smooth_axis1_d3, down_sample, median_filter_axis1_d3
from ..utils.io import save_dict_hdf5

from astropy import constants as const
from astropy.time import Time
from astropy import units as u
import scipy.interpolate as interp

import os
import re
from glob import glob
import copy


global n1, n2
global nB
global freql,freqh
global tc,tc_freq
global fname_part

global smooth #"mean", "gaussian",'poly'
global s_sigma #for freq (MHz)
global s_deg # for ploy
global med_filter # median filter
global med_size # median filter kernel size; odd number
global noise_mode

def plot_sep(nspec,outdir=None):
    '''
    nspec: class tcal_onoff, after run nspec.sep_onoff(...)
    '''
    from matplotlib import pyplot as plt
    plt.switch_backend('agg')
    
    fig= plt.figure(figsize=(15,4*3))
    data= nspec.data[:,nspec.use,:]
    len_freq= data.shape[1]
    wbin= np.min([200,len_freq//4])
    for i,ind_s in enumerate((len_freq-wbin)//4*np.array([1,2,3])):
        ax= fig.add_subplot(3,1,i+1)
        xx=np.mean(data[:200, ind_s:ind_s+wbin,0],axis=1)
        mjds= nspec.mjds[:200]
        plt.plot(mjds,xx,'k',label='ori')
        is_use= (nspec.on_mjds >= mjds.min()) & (nspec.on_mjds <= mjds.max())
        plt.plot(nspec.on_mjds[is_use], np.mean(nspec.on[is_use,ind_s:ind_s+wbin,0],axis=1),'r.',label='on')
        is_use= (nspec.off_mjds >= mjds.min()) & (nspec.off_mjds <= mjds.max())
        plt.plot(nspec.off_mjds[is_use], np.mean(nspec.off[is_use,ind_s:ind_s+wbin,0],axis=1),'b.',label='off')
        plt.ylabel('xx')
        plt.legend()
        plt.minorticks_on()
        #plt.title(f"{nspec.freq[ind_s]:.3f}|{nspec.vlsr[ind_s]:.3f}")
    plt.xlabel('mjd')
    plt.tight_layout()
    
    if outdir is None:
        outdir= './' + (Time(nspec.on_mjds[0],format='mjd')+8*u.hour).to_datetime().strftime('%Y%m%d') # mjd is UTC
        os.makedirs(outdir, exist_ok=True)
    plt.savefig(outdir+'/'+f"{os.path.basename(fname_part)[:-1]}-sep.pdf")


def test(outdir,chunk=None):
    if chunk is None:
        chunk=1
    try:
        nspec= tcal_onoff(fname_part,chunk,chunk+1)
    except:
        nspec= tcal_onoff(fname_part,chunk,chunk)
    nspec.sep_onoff(n1,n2,[freql,freqh],ave_cycle=True)
    nspec.sep_mjds_onoff(n1,n2,ave_cycle=True)
    plot_sep(nspec,outdir)

def load(start, stop, verbose= True, plot= False, outdir=None, dfactor=None, ave_cycle=False):

    if n1!=n2:
        ave_cycle=False
    nspec= tcal_onoff(fname_part,start,stop,verbose)
    nspec.sep_onoff(n1,n2,[freql,freqh],ave_cycle=ave_cycle)
    nspec.sep_mjds_onoff(n1,n2,ave_cycle=ave_cycle)
    freq= nspec.freq[nspec.use]
    on,off= nspec.on, nspec.off
    if med_filter:
        on= median_filter_axis1_d3(on, med_size)
        off= median_filter_axis1_d3(off, med_size)
    if dfactor is not None:
        freq= down_sample(freq[None,:,None],dfactor)[0,:,0]
        on= down_sample(on,dfactor)
        off= down_sample(off,dfactor)
    if n1!=n2:
        #gather spec in same period
        on_cyc= on.reshape(tuple([-1,n1]+list(on.shape[1:])))
        off_cyc= off.reshape(tuple([-1,n2]+list(off.shape[1:])))
        # smooth on-off each period
        n= min(n1,n2)
        if smooth=='gaussian':
            sigma= s_sigma/(np.nanmax(freq)-np.nanmin(freq))*len(freq)
            s_cal= smooth_axis1_d3(np.mean(on_cyc[:,:n]-off_cyc[:,:n],axis=1), method=smooth, sigma=sigma)
        if smooth=='poly':
            s_cal= smooth_axis1_d3(np.mean(on_cyc[:,:n]-off_cyc[:,:n],axis=1), method=smooth, x=freq, deg=s_deg)
        s_cal= s_cal[:,None,:,:]
        #T from tcal
        tc_inter= interp.interp1d(tc_freq,tc, axis=1, kind='linear',fill_value ='extrapolate')(freq)
        T_off= (off_cyc/s_cal).reshape([-1,]+list(off_cyc.shape[2:]))*tc_inter
        T_on= (on_cyc/s_cal).reshape([-1,]+list(on_cyc.shape[2:]))*tc_inter- tc_inter
    else:
        if smooth=='gaussian':
            sigma = s_sigma/(np.nanmax(freq)-np.nanmin(freq))*len(freq)
            s_cal= smooth_axis1_d3(on-off, method=smooth, sigma=sigma)
        if smooth=='poly':
            s_cal= smooth_axis1_d3(on-off, method=smooth, x=freq, deg=s_deg)
        tc_inter= interp.interp1d(tc_freq,tc, axis=1, kind='linear',fill_value ='extrapolate')(freq)
        t_per= tc_inter/s_cal
        T_off= off*t_per
        T_on= on*t_per-tc_inter
    
    if plot:
        plot_sep(nspec,outdir)
    res= [freq, np.vstack([T_on,T_off]).astype('float32'), np.hstack([nspec.on_mjds, nspec.off_mjds]),nspec.drop_num,len(nspec.tdata[0]), tc_inter]
    #res= [copy.deepcopy(freq), copy.deepcopy(np.vstack([T_on,T_off])), copy.deepcopy(np.hstack([nspec.on_mjds, nspec.off_mjds]))]
    #del(nspec)
    return res

def load_tc():
    global tc
    global tc_freq
    #global tc_  #del after
    tc_freq, tc_= read_tcal(nB, mode= noise_mode)
    is_use= (tc_freq<freqh) & (tc_freq>freql)
    tc_= tc_[:,is_use]
    tc_freq= tc_freq[is_use]
    if smooth=='gaussian':
        sigma = s_sigma/(np.nanmax(tc_freq)-np.nanmin(tc_freq))*len(tc_freq)
        tc= smooth_axis1_d3(tc_.T[None,:,:],method=smooth, sigma=sigma)
    if smooth=='poly':
        tc= smooth_axis1_d3(tc_.T[None,:,:],method=smooth, x=tc_freq, deg=s_deg)
def load_all(step, start_all=None, stop_all=None, sep_save=False, outdir=None, dfactor=None, ave_cycle=False, header=None):

    if start_all is None:
        start_all=1
    if stop_all is None:
        stop_all= len(glob(fname_part+'*.fits'))
    starts= list(range(start_all,stop_all+1,step))
    stops= list(range(start_all+step-1,stop_all,step))+[stop_all]
    
    res= []
    plot=True
    for i,j in zip(starts,stops):
        res_tmp= load(i,j,plot=plot,outdir=outdir,dfactor=dfactor,ave_cycle=ave_cycle)
        plot=False
        specs= {}
        specs['freq']= res_tmp[0]
        specs['T']= res_tmp[1]
        specs['mjd']= res_tmp[2]
        specs['Tcal']= res_tmp[5]
        drop_num=res_tmp[3]
        if len(starts)>1 and drop_num>0 and i==starts[0]:
            num_per= res_tmp[4]
            raise(ValueError(f"--step should be an integral multiple of {np.lcm(n1+n2,num_per%(n1+n2))//(num_per%(n1+n2))}."))
        
        if sep_save:
            if outdir is None:
                outdir= './' + (Time(specs['mjd'][0],format='mjd')+8*u.hour).to_datetime().strftime('%Y%m%d') # mjd is UTC
                os.makedirs(outdir, exist_ok=True)
            outname= outdir+'/'+f"{os.path.basename(fname_part)[:-1]}-specs_T_{i:04d}_{j:04d}.hdf5"
            print(f"saving to {outname}")
            save_dict_hdf5(outname,specs, header=header)
        else:
            res+= [res_tmp,]
    if not sep_save:
        res=np.asarray(res)
        specs= {}
        specs['freq']=res[0][0]
        specs['T']= np.vstack(res[:,1])
        specs['mjd']= np.hstack(res[:,2])
        specs['Tcal']= res[0][5]

        if outdir is None:
            outdir= './' + (Time(specs['mjd'][0],format='mjd')+8*u.hour).to_datetime().strftime('%Y%m%d') # mjd is UTC
            os.makedirs(outdir, exist_ok=True)
        outname= outdir+'/'+f"{os.path.basename(fname_part)[:-1]}-specs_T.hdf5"
        print(f"saving to {outname}")
        save_dict_hdf5(outname,specs, header=header)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(allow_abbrev=False)
    parser.add_argument('fname',
                        help='file name')
    parser.add_argument('-m', type=int,
                       help='time of Tcal_on divided by sampling time')
    parser.add_argument('-n', type=int,
                        help='time of Tcal_off divided by sampling time')
    parser.add_argument('--freql', type=float, default=1415.4,
                        help='lowest freq, default 1415.4')
    parser.add_argument('--freqh', type=float, default=1425.4,
                        help='highest freq, default 1425.4')
    parser.add_argument('--test', action='store_true',
                       help='only test the separation of Tcal on off')
    parser.add_argument('--step', type=int, default=5,
                       help='number of files loaded into memory every time, default 5')
    parser.add_argument('--start', type=int,
                       help='chunk number of start')
    parser.add_argument('--stop', type=int,
                       help='chunk number of stop')
    parser.add_argument('--sep_save', action='store_true',
                       help='every step save to a file')
    parser.add_argument('--outdir',
                       help='the directory to store output files.')
    
    parser.add_argument('--smooth', default='mean',
                        help="smooth method: 'mean','poly',gaussian',.. default:mean")
    parser.add_argument('--s_sigma', type=float, default=10,
                        help='sigma for gaussian smooth, default 10MHz')
    parser.add_argument('--s_deg', type=int, default=1,
                        help='Degree of the fitting polynomial, 0 equals mean; 1 is linear fit. default 1')
    parser.add_argument('--dfactor', type=int,
                        help='Down-sample the data before cal Teff')
    parser.add_argument('--med_filter', action='store_true',
                        help='median filter before Down-sample the data')
    parser.add_argument('--med_size', type=int, default=5,
                        help='median filter kernel size; odd number; default 5')
    
    parser.add_argument('--ave_cycle', action='store_true', 
                        help='average on and off spec in the same period')
    parser.add_argument('--noise_mode', default='high',
                        help='noise_mode, high or low')
   
    args = parser.parse_args()
    
    fname_part= re.sub('[0-9]{4}\.fits\Z','',args.fname)
    nB= int(re.findall(r'-M[0-1][0-9]',fname_part)[-1][2:])# used to read tcal
    n1, n2= args.m, args.n
    #freql, freqh= 1415.4, 1425.4
    freql, freqh= args.freql, args.freqh
    step= args.step
    sep_save=args.sep_save
    outdir= args.outdir
    if outdir is not None and not os.path.exists(outdir):
        raise(NotADirectoryError(outdir))
    start_all, stop_all= args.start, args.stop
    
    smooth= args.smooth #"mean", "gaussian","poly"
    s_sigma= args.s_sigma
    s_deg= args.s_deg
    if smooth== "gaussian" and (freqh- freql)/s_sigma <3:
        raise(ValueError('s_sigma is too larger for the input freq range'))
    if smooth=='mean':
        smooth='poly'
        s_deg=0
    dfactor= args.dfactor
    med_filter= args.med_filter
    med_size= args.med_size
    
    ave_cycle= args.ave_cycle
    noise_mode= args.noise_mode
    print('processing %s'%fname_part)
    print("Beam number:",nB)
    #print('beam: %d'%nB)
    print('m,n: %d, %d'%(n1,n2))
    print('freq range',freql, freqh)
    #record history
    from ..utils.io import rec_his
    header=rec_his(args=args)
    if args.test:
        test(outdir,start_all)
    else:
        load_tc()
        load_all(step,start_all, stop_all, sep_save, outdir, dfactor, ave_cycle,header=header)
