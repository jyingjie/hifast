#!/usr/bin/env python
# coding: utf-8
# Author: Nekomata  zmliu@nao.cas.cn
#修改nBs获取方式

from ast import arg
import numpy as np

import h5py
import os
from glob import glob
import copy
import re
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy import units as u
import scipy.interpolate as interp
from scipy.optimize import curve_fit

#for hifast@dev
from hifast.core.cal import CalOnOff, FastRawSpec
from hifast.core.radec import get_radec
from hifast.utils.io import MjdChanPolar_to_PolarMjdChan
from hifast.utils.io import PolarMjdChan_to_MjdChanPolar
from hifast.utils.misc import smooth1d
from hifast.utils.io import save_dict_hdf5
from hifast.utils.io import replace_nB
from hifast.utils.io import ArgumentParser, bool_fun, formatter_class
from hifast.utils.io import sub_patten

import sys


global nB
global radec
global fit_k
global rlim_src # max allowed radius[arcsec] to src
global freq
global obsdate
global outname
global flux
global nproc
global outdir
nproc = 1
 

def load_crd(obj,crds):
    if crds==None:
        from astroquery.simbad import Simbad
        try:
            sbo = Simbad.query_object(obj)
            crd = SkyCoord((sbo['RA'][0]+sbo['DEC'][0]),unit=(u.hourangle,u.deg))
        except:
            print('Calibrator can not find in Simbad.query_object, please input coordinates using ``--crd``.')
            sys.exit(1)
    else:
        crd= SkyCoord(*crds,unit=(u.deg))
    return crd

import csv

def load_flux_profile(calname,freq_key,fpparas):
    """
    freq_key should in MHz
    calname should be string
    
    """
    freq_key= freq_key/1000
    if fpparas is not None:
        a0,a1,a2,a3= fpparas
    else:
        # import pandas as pd
        # fpf= pd.read_csv(os.path.dirname(__file__) + '/data/FluxProfiles.csv')
        try:
            fpath = os.path.dirname(__file__) + '/data/FluxProfiles.csv'
            found = False
            with open(fpath, 'r') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    if row['name'] == calname:
                        a0 = float(row['a0'])
                        a1 = float(row['a1'])
                        a2 = float(row['a2'])
                        a3 = float(row['a3'])
                        found = True
                        break
            if not found:
                 raise ValueError("Calibrator name not found in FluxProfiles.csv")

            # narg= np.where(fpf['name']==calname)[0][0]
            # a0,a1,a2,a3= fpf['a0'][narg],fpf['a1'][narg],fpf['a2'][narg],fpf['a3'][narg]
        except:
            print('Can not find calibrator name in FluxProfile.csv. Please check calname or input flux profile parameters.')
            sys.exit(1)
    logS = a0+a1*np.log10(freq_key)+a2*(np.log10(freq_key))**2+a3*(np.log10(freq_key))**3
    flux= 10**(logS)
    return flux

#fitting method of Scan modes
def fit_curve(T,ra,figname=None):
    rcon= (ra>-0.025)&(ra<0.025)
    rcoff= ~((ra>-0.1)&(ra<0.1))
    tmax= np.max(T)
    tsys= np.mean(T[rcoff])
    if figname is not None:
        from matplotlib import pyplot as plt
        plt.figure(figsize=(8,4))
        plt.plot(ra, T, color='k')
        plt.scatter(ra[rcon], T[rcon], marker='.',color='r')
        plt.scatter(ra[rcoff],T[rcoff], marker='.',color='b')
        plt.xlabel('RA [deg]')
        plt.ylabel('T')
        plt.grid()
        plt.tight_layout()
        plt.savefig(figname)
    
    onbound=([tmax-0.05*np.abs(tmax)-tsys,-0.01,0,tsys-0.05*np.abs(tsys),],[tmax+0.05*np.abs(tmax)-tsys,0.01,1,tsys+0.05*np.abs(tsys),])
    onp,oncon = curve_fit(fit_scon,ra[rcon],T[rcon],bounds=onbound)
    offbound= ([-1,tsys-0.1*np.abs(tsys)],[1,tsys+0.1*np.abs(tsys)])
    offp,offcon = curve_fit(fit_scoff,ra[rcoff],T[rcoff],bounds=offbound)
    #ton= onp[0]+onp[3]
    #toff= offp[1]
    #return np.array([ton,toff,tmax,onp[1]])
    return tmax, onp, offp 
def fit_scon(x,a,b,c,d):
    return a*np.exp(-(x-b)**2/(2*c**2))+d
def fit_scoff(x,k,b):
    return k*x+b

#separate On Off time of Track modes
def load_track_time(mjd,t_src,nB,obsmode,n_cir):
    if obsmode=='OnOff':
        t_change = 30/60/60/24
        t_cir= t_src+t_change
        mjds= mjd[0]+np.arange(n_cir)*t_cir
        is_src,is_ref = np.full(len(mjd), False),np.full(len(mjd), False)
        for st in mjds:
            is_src = is_src | ((mjd > st) & (mjd < st+t_src))
            is_ref = is_ref | ((mjd > st+t_src+t_change) & (mjd < st+2*t_src+t_change))
        if nB==1:
            mjd_src,mjd_ref= is_src,is_ref
        else:
            mjd_src,mjd_ref= is_ref,is_src
    elif obsmode=='MultiBeamCalibration':
        t_change = 40/60/60/24
        t_cir= t_src+t_change
        mjds= mjd[0]+np.arange(20)*t_cir
        obslist= np.array([1,2,3,4,5,6,7,19,8,9,10,11,12,13,14,15,16,17,18])
        mjdon= mjds[np.where(obslist==nB)[0][0]+1]
        mjd_src= (mjdon<mjd)&(mjdon+t_src>mjd)                  
        if nB==1:
            mjd_src= mjd_src|(mjd<(mjds[0]+t_src))
        mjd_ref= (mjdon-t_cir>mjd)|(mjdon+t_cir<mjd)
    return mjd_src,mjd_ref

#mask rfi roughly
def rfi_mask(T,plot=False):
    Tsmt= smooth1d(T,'gaussian_fft',sigma_nchan,axis=0)
    Tdiff= np.abs(np.diff(Tsmt,axis=0))
    is_rfi= np.full(T.shape,False)
    dlimit= 0.015
    k= 500
    T_tmp= []
    for i in range(2):
        Tdiff_use=Tdiff[:,i]
        is_rfi_use= is_rfi[:,i]
        nlist= np.where(Tdiff_use>dlimit)[0]
        for j in range(len(nlist)):
            begin= nlist[j]-k
            end= nlist[j]+k
            if nlist[j]+1 in nlist:
                end= nlist[j]+k+1
            is_rfi_use[begin:end]= True
        Ttmp= T[~is_rfi_use,i]
        Ttmpitp= interp.interp1d(freq[~is_rfi_use],Ttmp, kind='linear', axis=0,fill_value ='extrapolate')(freq)
        T_tmp.append(Ttmpitp[:,np.newaxis])
        is_rfi[:,i]= is_rfi_use
    T_itp= np.concatenate(T_tmp,axis=1)
    T_itpsmt= smooth1d(T_itp,'gaussian_fft',sigma_nchan*10,axis=0)
    return T_itpsmt
def rfi_mask_pcal(p_cal,plot=False):
    pdiff= np.abs(np.diff(p_cal,axis=0))
    is_rfi= np.full(p_cal.shape,False)
    pdiffm= np.median(pdiff,axis=0)
    k= 500
    p_tmp= []
    for i in range(2):
        pdiff_use=pdiff[:,i]
        is_rfi_use= is_rfi[:,i]
        nlist= np.where(pdiff_use>pdiffm[i]*10)[0]
        for j in range(len(nlist)):
            begin= nlist[j]-k
            end= nlist[j]+k
            if nlist[j]+1 in nlist:
                end= nlist[j]+k+1
            is_rfi_use[begin:end]= True
        ptmp= p_cal[~is_rfi_use,i]
        ptmpitp= interp.interp1d(freq[~is_rfi_use],ptmp, kind='linear', axis=0,fill_value ='extrapolate')(freq)
        p_tmp.append(ptmpitp[:,np.newaxis])
        is_rfi[:,i]= is_rfi_use
    p_itp= np.concatenate(p_tmp,axis=1)
    return p_itp

# separate spectrum of calibrator
def sep_spe(fname,d,m,n):
    fname_part= fname
    para = {}
    para['n_delay'] = d
    para['n_on'] = m
    para['n_off'] = n
    para['start'] = 1
    para['stop'] = None
    para['frange'] = (frange[0], frange[1])
    para['verbose'] = True
    para['smooth'] ='gaussian'
    para['s_para'] = {'s_sigma':sigma}
    para['dfactor'] = None
    para['med_filter_size'] = None
    para['noise_mode'] = noise_mode
    para['noise_date'] = noise_date
    spec = CalOnOff(fname_part, **para)
    return spec

def load_pcal(pcal_s,mjd):
    mjd= mjd[:(pcal_s.shape[0])]
    if obsmode in ['Drift','DriftWithAngle','DecDriftWithAngle','MultiBeamOTF']:
        pcal_smt= np.mean(pcal_s,axis=0)
    else:
        ptime= load_track_time(mjd,t_src,nB,obsmode,n_cir=ncir)
        pcal_smt= np.mean(pcal_s[ptime[1]],axis=0)
    pcal_smt= smooth1d(pcal_smt,'gaussian_fft',sigma_nchan,axis=0)
    return pcal_smt

def plot_radec_near_cbr(ra, dec, ra_cbr, dec_cbr, figname, r_max=0.15):

    from matplotlib import pyplot as plt
    plt.figure(figsize=(5,5))
    
    is_use = (ra < ra_cbr+r_max ) & (ra > ra_cbr-r_max)
    is_use &= (dec < dec_cbr+r_max ) & (dec > dec_cbr-r_max)
    
    plt.scatter(ra[is_use], dec[is_use], marker='.', s=1, color='k')
    plt.scatter(ra_cbr, dec_cbr, marker='.', s=7, color='r')
    plt.xlabel('RA [deg]')
    plt.ylabel('DEC [deg]')
    plt.grid()
    plt.tight_layout()
    plt.savefig(figname)

def load_Tsc(spec,T,mjd):
    global radec
    global fit_k
    global rlim_src # max allowed radius[arcsec] to src
    if nB==1:
        print(f'Calculating RA-DEC for Beam {list(nBs)}')
        radec = get_radec(mjd, guess_str=fname,tol=5,nBs=list(nBs),nproc=nproc)
        save_dict_hdf5(outname+'-radec.hdf5',radec) 
        print('Saved RA-DEC to '+outname+'-radec.hdf5')
    else:
        try:
            _,_ = radec[f'ra{nB}'], radec[f'dec{nB}']
        except:
            print(f'Calculating RA-DEC for Beam {nB:02d}')
            radec = get_radec(mjd, guess_str=fname,tol=5,nBs=[nB,])
    ra0,dec0 = radec[f'ra{nB}'], radec[f'dec{nB}']
    
    radec_figname = outname + f'-M{nB:02d}-radec_near_src.png'
    plot_radec_near_cbr(ra0, dec0, ccrd.ra.deg, ccrd.dec.deg, figname=radec_figname, r_max=0.15)
    
    if obsmode in ['Drift','DriftWithAngle','DecDriftWithAngle','MultiBeamOTF']:
        T= smooth1d(T,'gaussian_fft',sigma_nchan,axis=1)
        ra= ra0-ccrd.ra.deg
        dec= dec0-ccrd.dec.deg
        data_cut= (ra>-0.5)&(ra<0.5)&(dec<0.05)&(dec>-0.05)
        Tuse= T[data_cut]
        rause= ra[data_cut]
        ONp= []
        OFFp= []
        TMAX= []
        freq_key= np.arange(frange[0]+1,frange[1],1)
        fit_k= 0
        for i in range(len(freq_key)):
            f_cut= (freq>freq_key[i]-1)&(freq<freq_key[i]+1)
            T_fit= np.mean(Tuse[:,f_cut,:],axis=1)
            try:
                if i==0:
                    figname = outname + f'-M{nB:02d}-freq_{freq_key[i]:.2f}.pdf'
                else:
                    figname = None
                tmaxXX, onpXX, offpXX  = fit_curve(T_fit[:,0],rause,figname=figname)
                tmaxYY, onpYY, offpYY  = fit_curve(T_fit[:,1],rause,figname=figname)
            except:                
                print('Fit failed in '+str(int(freq_key[i])) + f'MHz in beam {nB}, please check.')
                tmaxXX, onpXX, offpXX= np.nan,[np.nan]*4,[np.nan]*2
                tmaxYY, onpYY, offpYY= np.nan,[np.nan]*4,[np.nan]*2
                fit_k+=1
                if fit_k>=int(0.1*len(freq_key)):
                    print(f'Fit failed at too many Freq in beam {nB}. Please check')
                    sys.exit(1)
            ONp.append(np.array([onpXX,onpYY],dtype='float64'))
            OFFp.append(np.array([offpXX,offpYY],dtype='float64'))
            TMAX.append(np.array([tmaxXX,tmaxYY],dtype='float64'))
        ONp= np.array(ONp,dtype='float64')
        OFFp= np.array(OFFp,dtype='float64')
        TMAX= np.array(TMAX,dtype='float64')
        Tsrc= interp.interp1d(freq_key, ONp[:,:,0]+ONp[:,:,3], kind='linear', axis=0,fill_value ='extrapolate')(freq)
        Tref= interp.interp1d(freq_key, OFFp[:,:,1], kind='linear', axis=0,fill_value ='extrapolate')(freq)
        fit_k= fit_k/len(freq_key)
    else:
        crds= SkyCoord(ra0,dec0,unit=u.deg)          
        mjd_src,mjd_ref= load_track_time(mjd,t_src,nB,obsmode,n_cir=ncir)
        T_src,T_ref= T[mjd_src],T[mjd_ref]
        T_select = True
        if T_select:
            print(f'Further restrict the on-source spectra within {args.src_drange}')
            Tsrc,_= radec_select(T_src,crds[mjd_src], args.src_drange)
            fit_k= 1- Tsrc.shape[0]/T_src.shape[0]
            print(f'Further restrict the off-source spectra within {args.ref_drange}')
            Tref,_= radec_select(T_ref,crds[mjd_ref], args.ref_drange)
        else:
            Tsrc,Tref= T_src,T_ref
            fit_k= 1
        Tsrc= np.mean(Tsrc,axis=0)
        Tref= np.mean(Tref,axis=0)
        ONp,OFFp,TMAX= Tsrc,Tref,np.nan
    return Tsrc,Tref,ONp,OFFp,TMAX

def radec_select(T,crds,seprange):
    seps= ccrd.separation(crds).to_value(u.arcsec)
    is_use= (seps>=seprange[0])&(seps<=seprange[1])
    num_used= T[is_use].shape[0]
    if num_used==0:
        print(f'Beam {nB:02d}: No spectra inside {seprange}, please check input parameters.')
        sys.exit(1)
    else:
        print(f'Beam {nB:02d}: Input {T.shape[0]} spectra, {num_used} spectra is in {seprange}')
    return T[is_use],is_use

def load_data(fname):
    global nB
    dataname = replace_nB(fname, nB)
    spec= sep_spe(dataname,d,m,n)
    global freq
    global sigma_nchan
    freq= spec.freq_use
    sigma_nchan =  int(sigma/(np.max(freq) - np.min(freq))*(len(freq)-1))
    mjd= spec.get_mjds()
    global obsdate
    obsdate= (Time(np.median(mjd),format='mjd').to_value('iso','date_hm'))
    obsdate= re.sub('\D','',obsdate)
    
    global outdir
    outdir = sub_patten(outdir, date=obsdate+'UTC', project=calname)
    outdir = os.path.expanduser(outdir)
    if outdir!='' and (not os.path.exists(outdir)):
        print(f'outdir {outdir} not exists. Create it now')
        os.makedirs(outdir, exist_ok=True)
    global outname
    outname= f"{outdir}/{calname}-{obsmode}-{obsdate}UTC"
    
    mjd_on = mjd[spec.inds_on]
    mjd_off = mjd[spec.inds_off]
    
    Tcal_s= spec.get_Tcal_s()
    spec.gen_out_name_base(outdir)
    spec.plot = True
    spec.sep_on_off_inds()
    p_on = spec.get_field(spec.inds_on, 'DATA',)
    p_off= spec.get_field(spec.inds_off, 'DATA', close_file=False)
    #c_on, c_off, pcal_s,_= spec.get_count_tcal(spec.inds_on,spec.inds_off)
    # used power of cal not smoothed
    pcal_s = spec._get_cal_power(spec.inds_ton, spec.inds_toff_bef, spec.inds_toff_aft)
    p_cal= load_pcal(pcal_s,mjd_on[::m])
    p_cal= rfi_mask_pcal(p_cal)
    T_off= p_off*Tcal_s/p_cal

#Save T files for checking
    if saveT==True:
        odata= {}
        odata['T']= T_off
        odata['freq']= freq
        odata['mjd']= mjd_off
        odata['pcal']= p_cal
        odata['Tcal']= Tcal_s
        odata['Tcal_file']= spec.tcal_file
        save_dict_hdf5(outname+f'-M{nB:02d}_T.hdf5',odata) 
        print('Saved the temperature data to '+outname+f'-M{nB:02d}_T.hdf5')

    Tsrc,Tref,ONp,OFFp,TMAX= load_Tsc(spec,T_off,mjd_off,)
    Tsc= Tsrc-Tref
    Tsc_rfi= rfi_mask(Tsc)
    global flux
    flux= load_flux_profile(calname,freq,fpparas)
    flux= flux[:,np.newaxis]
    K_Jy= Tsc_rfi/flux
    K_Jy= K_Jy[np.newaxis,:]
    return K_Jy,Tcal_s,spec.tcal_file,mjd_off,ONp,OFFp,TMAX



sep_line = '##'+'#'*70+'##'
parser = ArgumentParser(prog=f"python -m hifast.{os.path.basename(sys.argv[0])[:-3]}", formatter_class=formatter_class, 
                        allow_abbrev=False,
                        description='Processing the calibrator data', )
parser.add_argument('fname',
                    help="File path of raw data, to be used with the `--nBs` parameter: '*M01*-0001.fits'.")
parser.add_argument('--nBs', type=int, nargs='*',
                    help='beam numbers. Input int number, seperated by space, such as `1 3 11`. \
                    If no numbers are provided, process all beams when the `--obsmode` is `MultiBeamCalibration` or `MultiBeamOTF`, otherwise only beam 1 will be processed.', default=None)
parser.add_argument('--outdir', type=str,
                   help='output directory',default=None)

group = parser.add_argument_group(f'*Calibrator setting\n{sep_line}')
group.add_argument('--obsmode',type=str, choices=['MultiBeamCalibration', 'OnOff', 'MultiBeamOTF','DriftWithAngle','DecDriftWithAngle','Drift'],
                    help='Observation mode. Supports `MultiBeamCalibration`, `OnOff`, and some special scanning.',)
group.add_argument('--calname', type=str, required=True,
                   help='calibrator name. e.g. 3C48',)
group.add_argument('--crd', type=float, nargs=2,
                   help='Calibrator coordinates in deg. Two float numbers, for RA and DEC. \
                   If None, use Simbad.query_object to query by `--calname`.', default=None)

calnames = "3C48,3C286,..."
try:
    calnames_list = []
    with open(os.path.dirname(__file__) + '/data/FluxProfiles.csv', 'r') as f:
         reader = csv.DictReader(f)
         for row in reader:
             calnames_list.append(row['name'])
    if calnames_list:
        calnames = ','.join(calnames_list)
except:
    pass

# fpf= pd.read_csv(os.path.dirname(__file__) + '/data/FluxProfiles.csv')
# calnames = ','.join(fpf['name'])
group.add_argument('--fluxProfilePara', type=float, nargs=4,
                   help=f'The function of calibrator flux with respect to frequency. Function of {calnames} are already embeded in code.\
                   Others must be specified and include four floating-point numbers: a0, a1, a2, and a3. \
                   The function will be represented by the equation `a0 + a1*log10(f) + a2*log10(f)**2 + a3*log10(f)**3`.',default= None)


group = parser.add_argument_group(f'*Parameters same with hifast.sep\n{sep_line}')
group.add_argument('-d','--d', '--n_delay', type=int,
                   help='time of delay divided by sampling time')
group.add_argument('-m','--m', '--n_on', type=int,
                   help='time of Tcal_on divided by sampling time')
group.add_argument('-n','--n', '--n_off', type=int,
                    help='time of Tcal_off divided by sampling time')
group.add_argument('--frange', type=float, nargs=2,
                   help='freq range',default= [1000,1500] )
group.add_argument('--noise_mode', default='high', choices=['high','low'],
                    help='noise_mode, high or low')
group.add_argument('--noise_date', default='auto',type=str,
                    help='noise obs date, default auto')
group.add_argument('--smt_sigma', type=float,
                   help='smooth sigma (in MHz) of gaussian smooth along frequency',default=1)

group = parser.add_argument_group(f'*Parameters in MultiBeamCalibration or OnOff\n{sep_line}')
group.add_argument('--t_src', type=int,
                   help='Tracking time of On-Source in second. (It is used in conjunction with Off-Source and change time to determine the On/Off source spectra. \
                   In the `OnOff` mode, the Off-Source time used is identical to the On-Source time, with a change time of 30 seconds. \
                   In the `MultiBeamCalibration` mode, the change time is 40 seconds.)', default=60)
group.add_argument('--n_cir', '--n_repeat', type=int, default=1,
                   help='The number of on-source off-source cycles in `OnOff`. It is 1 in `MultiBeamCalibration`')
# group.add_argument('--T_select', type=bool_fun, choices=[True, False], default='True',
#                    help='Select spectra within `--rlim_src` from calibrater source if `--obsmode MultiBeamCalibration`')
group.add_argument('--src_drange', nargs=2, type=float, default=[0., 20.],
                   help='Further limit the distance of On-Source spectra to the Calibrator coordinates in this range. unit is arcsecond',)
group.add_argument('--ref_drange', nargs=2, type=float, default=[180., 3600.],
                   help='Further limit the distance of Off-Source spectra to the Calibrator coordinates in this range. unit is arcsecond',)
group = parser.add_argument_group(f'*Others\n{sep_line}')
group.add_argument('--saveT', type=bool_fun, choices=[True, False], default='False',
                   help='Save T files or not',)

if __name__ == '__main__':

    args = parser.parse_args()
    fname= args.fname
    obsfile= fname.split('/')[-1]
    print('Flux Calibration of '+obsfile)
    d, m, n = args.d, args.m, args.n
    frange= args.frange
    noise_mode = args.noise_mode
    noise_date = args.noise_date
    sigma= args.smt_sigma
    obsmode= args.obsmode
    modes= ['Drift','MultiBeamCalibration','OnOff','MultiBeamOTF','DriftWithAngle','DecDriftWithAngle']
    if obsmode not in modes:
        print('Please check observation mode.')
        sys.exit(1) 
    CalinBs= args.nBs 
    if CalinBs== None:
        if obsmode in ['MultiBeamCalibration','MultiBeamOTF',]:
            nBs= np.arange(1,20)
        else:
            nBs= [1,]
    else:
        nBs= CalinBs
    print('Calibrating beams:')
    print(nBs)
    outdir= args.outdir
    if outdir==None:
        outdir= '.'
    calname= args.calname
    ccrd = load_crd(calname,args.crd)
    fpparas= args.fluxProfilePara
    t_src= args.t_src/3600/24#change into days
    ncir= args.n_cir
    T_select=args.T_select
    rlim_src = args.rlim_src
    saveT= args.saveT
    oname= f'{outdir}/{calname}-{obsmode}-'

    outdata={}
    for i in nBs:
        global nB
        nB= i
        K_Jy,Tcal_s,tcal_file,mjd,ONp,OFFp,TMAX= load_data(fname)
        ra,dec,ZD = radec[f'ra{i}'], radec[f'dec{i}'],radec[f'ZD{i}']
        outdata[f'M{nB:02d}']= K_Jy
        outdata[f'Tcal{nB}']= Tcal_s
        outdata[f'ra{i}']= ra
        outdata[f'dec{i}']= dec
        outdata[f'ZD{i}']= ZD
        outdata[f'ONp{i}']= ONp
        outdata[f'OFFp{i}']= OFFp
        outdata[f'Tmax{i}']= TMAX
        print(f'Finished calibration of M{nB:02d}, {fit_k*100:.3f}% input data has bean ignored.')
        print('#'*50)
        
    outdata['Tcal_file']= tcal_file
    print('Used Tcal of '+tcal_file)
    outdata['freq']= freq
    outdata['cal_flux']=flux
    outdata['mjd']= mjd
    if len(nBs)==19 :
        outname= outname+'-FluxGain-All.hdf5'
    else:
        outname= outname+'-FluxGain.hdf5'
    save_dict_hdf5(outname,outdata) 
    print('Saved to '+outname)
