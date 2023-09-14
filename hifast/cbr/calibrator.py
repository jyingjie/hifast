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

def load_crd(obj,crds):
    if crds==None:
        from astroquery.simbad import Simbad
        try:
            sbo = Simbad.query_object(obj)
            crd = SkyCoord((sbo['RA'][0]+sbo['DEC'][0]),unit=(u.hourangle,u.deg))
        except:
            print('Calibrator can not find in Simbad.query_object(), please input coordinates! ')
            os._exit(0)
    else:
        crd= SkyCoord(*crds,unit=(u.deg))
    return crd

def load_flux_profile(calname,freq_key,fpparas):
    """
    freq_key should in MHz
    calname should be string
    
    """
    freq_key= freq_key/1000
    if fpparas!= None:
        a0,a1,a2,a3= fpparas
    else:
        import pandas as pd
        fpf= pd.read_csv(os.path.dirname(__file__) + '/data/FluxProfiles.csv')
        try:
            narg= np.where(fpf['name']==calname)[0][0]
            a0,a1,a2,a3= fpf['a0'][narg],fpf['a1'][narg],fpf['a2'][narg],fpf['a3'][narg]
        except:
            print('Can not find calibrator name in FluxProfile.csv. Please check calname or input flux profile parameters!')
            os._exit()
    logS = a0+a1*np.log10(freq_key)+a2*(np.log10(freq_key))**2+a3*(np.log10(freq_key))**3
    flux= 10**(logS)
    return flux

#fitting method of Scan modes
def fit_curve(T,ra,):
    rcon= (ra>-0.025)&(ra<0.025)
    rcoff= ~((ra>-0.18)&(ra<0.18))
    tmax= np.max(T)
    tsys= np.mean(T[rcoff])
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
            mjd_src,mjd_ref= is_ref,is_src,
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
    Tsmt= smooth1d(T,'gaussian_fft',sigma,axis=0)
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
    T_itpsmt= smooth1d(T_itp,'gaussian_fft',sigma*10,axis=0)
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
    para['s_para'] = {'s_sigma':1}
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
    pcal_smt= smooth1d(pcal_smt,'gaussian_fft',sigma,axis=0)
    return pcal_smt

def load_Tsc(spec,T,mjd):
    global radec
    global fit_k
    global rlim_src # max allowed radius[arcsec] to src
    if nB==1:
        radec = get_radec(mjd, guess_str=fname,tol=30,nBs=list(nBs))
        save_dict_hdf5(outname+'-radec.hdf5',radec) 
        print('Saved to '+outname+'-radec.hdf5')
    else:
        radec = get_radec(mjd, guess_str=fname,tol=30,nBs=[nB,])
    ra0,dec0 = radec[f'ra{nB}'], radec[f'dec{nB}']
    if obsmode in ['Drift','DriftWithAngle','DecDriftWithAngle','MultiBeamOTF']:
        T= smooth1d(T,'gaussian_fft',sigma,axis=1)
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
                tmaxXX, onpXX, offpXX  = fit_curve(T_fit[:,0],rause,)
                tmaxYY, onpYY, offpYY  = fit_curve(T_fit[:,1],rause,)
            except:                
                print('Fit failed in '+str(int(freq_key[i]))+'MHz, please check.')
                tmaxXX, onpXX, offpXX= np.nan,[np.nan]*4,[np.nan]*2
                tmaxYY, onpYY, offpYY= np.nan,[np.nan]*4,[np.nan]*2
                fit_k+=1
                if fit_k>=int(0.1*len(freq_key)):
                    print(f'Fit failed at too many Freq. Please check input parameters!')
                    os._exit(0)
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
        if T_select==True:
            Tsrc,_= radec_select(T_src,crds[mjd_src],[0., rlim_src])
            fit_k= 1- Tsrc.shape[0]/T_src.shape[0]
            Tref,_= radec_select(T_ref,crds[mjd_ref],[180,36000])
        else:
            Tsrc,Tref= T_src,T_ref
            fit_k= 1
        Tsrc= np.mean(Tsrc,axis=0)
        Tref= np.mean(Tref,axis=0)
        ONp,OFFp,TMAX= Tsrc,Tref,np.nan
    return Tsrc,Tref,ONp,OFFp,TMAX

def radec_select(T,crds,seprange):
    seps= ccrd.separation(crds).to_value(u.arcsec)
    is_use= (seps>seprange[0])&(seps<seprange[1])
    num_used= T[is_use].shape[0]
    if num_used==0:
        print('No point inside source area, please check input parameters!')
        os._exit(0)
    else:
        print(f'Input {T.shape[0]} points, used {num_used} points')
    return T[is_use],is_use

def load_data(fname,nB):
    dataname = replace_nB(fname, nB)
    spec= sep_spe(dataname,d,m,n)
    global freq
    freq= spec.freq_use
    mjd= spec.get_mjds()
    global obsdate
    obsdate= (Time(np.median(mjd),format='mjd').to_value('iso','date_hm'))
    obsdate= re.sub('\D','',obsdate)
    global outname
    outname= oname+obsdate
    mjd_on = mjd[spec.inds_on]
    mjd_off = mjd[spec.inds_off]
    
    Tcal_s= spec.get_Tcal_s()
    spec.sep_on_off_inds()
    p_on = spec.get_field(spec.inds_on, 'DATA',)
    p_off= spec.get_field(spec.inds_off, 'DATA', close_file=False)
    c_on, c_off, pcal_s,_= spec.get_count_tcal(spec.inds_on,spec.inds_off)
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
        print('Saved to '+outname+f'-M{nB:02d}_T.hdf5')

    Tsrc,Tref,ONp,OFFp,TMAX= load_Tsc(spec,T_off,mjd_off,)
    Tsc= Tsrc-Tref
    Tsc_rfi= rfi_mask(Tsc)
    global flux
    flux= load_flux_profile(calname,freq,fpparas)
    flux= flux[:,np.newaxis]
    K_Jy= Tsc_rfi/flux
    K_Jy= K_Jy[np.newaxis,:]
    return K_Jy,Tcal_s,spec.tcal_file,mjd_off,ONp,OFFp,TMAX

def bool_fun(s):
    """
    used as the type of parser.add_argument
    """
    s = str(s)
    if s == 'True' or s.lower() == 'yes':
        return True
    elif s == 'False' or s.lower() == 'no':
        return False
    else:
        raise(ValueError("input must be 'True'('yes') or 'False'('no')"))

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('fname',
                        help='file name')
    parser.add_argument('-d','--d', '--n_delay', type=int, required=True,
                       help='time of delay divided by sampling time')
    parser.add_argument('-m','--m', '--n_on', type=int, required=True,
                       help='time of Tcal_on divided by sampling time')
    parser.add_argument('-n','--n', '--n_off', type=int, required=True,
                        help='time of Tcal_off divided by sampling time')
    parser.add_argument('--frange', type=float, nargs=2,
                       help='freq range',default= [1000,1500] )
    parser.add_argument('--noise_mode', default='high', choices=['high','low'],
                        help='noise_mode, high or low')
    parser.add_argument('--noise_date', default='auto',type=str,
                        help='noise obs date, default auto')
    parser.add_argument('--smt_sigma', type=float,
                       help='smooth sigma',default= 1)
    parser.add_argument('--obsmode',type=str,
                        help='observation mode',)
    parser.add_argument('--nBs', type=int, nargs='*',
                        help='beam numbers',default=None)
    parser.add_argument('--outdir', type=str,
                       help='output path',default=None)

    parser.add_argument('--calname', type=str, required=True,
                       help='calibrator name',)
    parser.add_argument('--crd', type=float, nargs=2,
                       help='calibrator coordinate in deg',default=None)
    parser.add_argument('--fluxProfilePara', type=float, nargs=4,
                       help='calibrator flux',default= None)
    parser.add_argument('--t_src', type=int,
                       help='tracking time of On Source',default=60)
    parser.add_argument('--n_cir', type=int,default=1,
                       help='tracking time of On Source')
    parser.add_argument('--T_select', type=bool_fun, choices=[True, False], default='True',
                       help='Select T with radec if `--obsmode MultiBeamCalibration`')
    parser.add_argument('--rlim_src', type=float, default=16,
                       help='if `--T_select`, maximum allowed offset (in arcsec) to calibrator. Default is 16.', )
    parser.add_argument('--saveT', type=bool_fun, choices=[True, False], default='True',
                       help='Save T files or not',)

    args = parser.parse_args()
    fname= args.fname
    obsfile= fname.split('/')[-1]
    print('Flux Calibration of '+obsfile)
    d, m, n = args.d, args.m, args.n
    frange= args.frange
    noise_mode = args.noise_mode
    noise_date = args.noise_date
    sigma= int(args.smt_sigma*65536/500)
    obsmode= args.obsmode
    modes= ['Drift','MultiBeamCalibration','OnOff','MultiBeamOTF','DriftWithAngle','DecDriftWithAngle']
    if obsmode not in modes:
        print('Please check observation mode!')
        os._exit(0) 
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
        K_Jy,Tcal_s,tcal_file,mjd,ONp,OFFp,TMAX= load_data(fname,i)
        ra,dec,ZD = radec[f'ra{i}'], radec[f'dec{i}'],radec[f'ZD{i}']
        outdata[f'M{nB:02d}']= K_Jy
        outdata[f'Tcal{nB}']= Tcal_s
        outdata[f'ra{i}']= ra
        outdata[f'dec{i}']= dec
        outdata[f'ZD{i}']= ZD
        outdata[f'ONp{i}']= ONp
        outdata[f'OFFp{i}']= OFFp
        outdata[f'Tmax{i}']= TMAX
        print(f'Finished calibration of M{nB:02d}, {fit_k*100}% input data has bean ignored.')
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
