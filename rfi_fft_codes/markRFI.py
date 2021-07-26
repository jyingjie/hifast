#!/usr/bin/env python
# coding: utf-8

# author: Xu Chen 2021.06

import sys
import os

from astropy import log
from matplotlib import pyplot as plt 

import numpy as np
import h5py
from copy import deepcopy
from tqdm import tqdm

def rms(data,vel,rms_vrange=None):
    if rms_vrange is not None:
        is_use = (vel > rms_vrange[0])&(vel < rms_vrange[1])
        if len(data.shape) == 1:
            data = data[is_use]
        elif len(data.shape) ==2:
            data = data[:,is_use]
        
    ans = np.sqrt(np.nanmean((data)**2,axis = -1))   
    return ans

def real_rms(data,vel,sigma,rms_vrange=None):
    from scipy.ndimage import gaussian_filter1d
    if rms_vrange is not None:
        is_use = (vel > rms_vrange[0])&(vel < rms_vrange[1])
        if len(data.shape) == 1:
            data = data[is_use]
        elif len(data.shape) ==2:
            data = data[:,is_use]
            
    g = gaussian_filter1d(data, sigma, axis = -1)
    ans = np.sqrt(np.nanmean((data - g)**2,axis = -1))        
    
    return ans

def real_std(data,vel,sigma,rms_vrange=None):
    from scipy.ndimage import gaussian_filter1d
    if rms_vrange is not None:
        is_use = (vel > rms_vrange[0])&(vel < rms_vrange[1])
        if len(data.shape) == 1:
            data = data[is_use]
        elif len(data.shape) ==2:
            data = data[:,is_use]
            
    g = gaussian_filter1d(data, sigma, axis = -1)
    ans = np.std(data - g,axis = -1)    
    
    return ans


def find_local_peak(rfi1,freq1,RMS,distance=5,):
    from scipy import signal
    loc_max = signal.find_peaks(rfi1,distance = distance)[0]
    
    loc_shape = loc_max.shape[0]
    if loc_shape == 1:
        loc = loc_max
    elif loc_shape > 1:
        loc = loc_max[np.argmax(rfi1[loc_max])]
    elif loc_shape == 0:
        print("No local peak is found")
        return np.nan,np.nan
        
    freq_max = np.float(freq1[loc])
    flux_peak = np.float(rfi1[loc] - RMS)
    
    return freq_max,flux_peak

def fit_line(rfi,freq,half_factor,**kwargs):
    
    freq_peak,flux_peak = find_local_peak(rfi,freq,**kwargs)
    if np.sum(np.isnan([freq_peak,flux_peak])) > 0:
        return np.nan,np.nan

    limd,limu = flux_peak*.15,flux_peak*.85
    half_peak = flux_peak * half_factor

    left_use = (freq > freq_peak)&(rfi > limd)&(rfi < limu)
    righ_use = (freq < freq_peak)&(rfi > limd)&(rfi < limu)

    freq_left = freq[left_use]
    flux_left = rfi[left_use]
    freq_righ = freq[righ_use]
    flux_righ = rfi[righ_use]
    
    try:
        funcl = np.polyfit(flux_left,freq_left,1)
        f_left = np.poly1d(funcl)
        funcr = np.polyfit(flux_righ,freq_righ,1)
        f_righ = np.poly1d(funcr)
    except TypeError as TE:
        #print(TE)
        return np.nan,np.nan

    half_v_left = f_left(half_peak)
    half_v_righ = f_righ(half_peak)
    
    if (half_v_left <= min(freq))|(half_v_righ >= max(freq)):
        return np.nan,np.nan

    return half_v_left, half_v_righ

def get_startend(is_rfi,rfi_width_lim = None,ext_sec= 0):
    if is_rfi.any() == False:
        raise ValueError("Input array has no Trues.")
    
    if is_rfi[0] == True:
        is_rfi[0] = False
        
    starend = np.hstack((np.diff(is_rfi+0),0))
    start_ = np.where(starend==1)[0]
    end_ = np.where(starend==-1)[0]
    
    if end_[0] < start_[0]:
        end_ = np.delete(end_,0)
    if end_[-1] < start_[-1]:
        start_ = np.delete(start_,-1)
    
    if rfi_width_lim is None:
        start = start_ - ext_sec
        end = end_ + ext_sec
    else:
        if isinstance(rfi_width_lim,int)|isinstance(rfi_width_lim,float):
            starend_use = (end_-start_)>rfi_width_lim
        elif len(rfi_width_lim) == 2:
            starend_use = ((end_-start_)>rfi_width_lim[0])&((end_-start_)<rfi_width_lim[1])

        if starend_use.any() == False:
            raise ValueError("No Trues meet width condition.")

        start = start_[starend_use] - ext_sec
        end = end_[starend_use] + ext_sec

    if min(start) < 0:
        start[start<0] = 0 
    N = len(is_rfi) - 1 
    if max(end) > N:
        end[end>N] =N
    
    return start,end


def find_center(spec,freq,is_rfi,RMS,freq_thr = .5,freq_step = 8.1,check = True,
                rfi_fit_use = 'two groups',**kwargs):

    start,end = get_startend(is_rfi,**kwargs)
    # big rfi position
    fc0 = [];peaks = []
    for s,e in zip(start,end):
        rfi1 = spec[s:e]
        freq1 = freq[s:e]
        freq_peak,flux_peak = find_local_peak(rfi1,freq1,RMS=RMS)
        #half_v_left, half_v_righ = fit_line(rfi1,freq1,half_factor =.5,RMS=RMS)
        #if np.sum(np.isnan([half_v_left, half_v_righ])) == 0:
            #vc = np.float((half_v_left + half_v_righ)/2)
            #if (vc > freq1[0])&(vc < freq1[-1]):  
        vc = freq_peak
        fc0 += [vc,]; peaks += [flux_peak,]
    
    peaks = np.array(peaks)
    # from the biggest one to guess rfi set 1
    fc00 = fc0[np.argmax(peaks)]
    
    def rfi_set1(fc0,fc00,freq,freq_step,freq_thr,peaks):
        fcenter1 = []; peak1 = []
        former = int((fc00 - freq[0])//freq_step + 10)
        fc1 = np.arange(fc00-freq_step*former,1440,freq_step)
        
        fcenter2_loc = []
        ii = -1
        for j in range(len(fc0)):
            vc = fc0[j]
            ii += 1
            if np.min(np.abs(vc - fc1)) < freq_thr:
                # record rfi set 1 center
                fcenter1 += [vc,]; peak1 += [peaks[j],]
            else:   
                fcenter2_loc += [int(ii),]
                
        fcenter1 = np.array(fcenter1)
        if check:
            peak1 = np.array(peak1)
            _,fcenter1 = check_repeat_center(fcenter1,freq_step=8.1,x_start = '1st point',peaks = peak1)        
                   
        return fcenter1,fcenter2_loc,fc1
    
    fcenter1,fcenter2_loc,fc1 = rfi_set1(fc0,fc00,freq,freq_step,freq_thr,peaks)
    
    if len(fcenter1)<2:
        log.warning(f"rfi set 1 only has {len(fcenter1)} items, please check if it is a high-flux source or non-period RFI. Try again...")
        peaks[np.argmax(peaks)] = 0
        # find another one
        fc00 = fc0[np.argmax(peaks)]
        fcenter1,fcenter2_loc,fc1 = rfi_set1(fc0,fc00,freq,freq_step,freq_thr,peaks)
    
    # guess rfi set 2 besides set 1                
    peak_use = np.full(len(peaks),False)
    peak_use[fcenter2_loc] = True
    peaks[~peak_use] = 0
    fc01 = fc0[np.argmax(peaks)]

    def rfi_set2(fc0,fc1,fc01,freq,freq_step,freq_thr,peaks):
        fcenter2 = []; peak2 = []
        former = int((fc01 - freq[0])//freq_step + 10)
        fc2 = np.arange(fc01-freq_step*former,1440,freq_step)

        fcenter3_loc = []
        ii = -1
        for j in range(len(fc0)):
            vc = fc0[j]
            ii += 1
            if np.min(np.abs(vc - fc2)) < freq_thr:
                fcenter2 += [vc,]; peak2 += [peaks[j],]
            elif (np.min(np.abs(vc - fc1)) >= freq_thr)&(np.min(np.abs(vc - fc2)) >= freq_thr):
                fcenter3_loc += [int(ii),]
        
        fcenter2 = np.array(fcenter2)
        if check:
            peak2 = np.array(peak2)
            _,fcenter2 = check_repeat_center(fcenter2,freq_step=8.1,x_start = '1st point',peaks = peak2)         
        
        return fcenter2,fcenter3_loc
    
    fcenter2,fcenter3_loc = rfi_set2(fc0,fc1,fc01,freq,freq_step,freq_thr,peaks)
    
    if len(fcenter2)<2:
        log.warning(f"rfi set 2 only has {len(fcenter2)} items, please check if it is a high-flux source or non-period RFI. Try again...")
        peaks[np.argmax(peaks)] = 0
        fc01 = fc0[np.argmax(peaks)]
        # find another one 
        fcenter2,fcenter3_loc = rfi_set2(fc0,fc1,fc01,freq,freq_step,freq_thr,peaks)
    
    fcenter3 = []
    if rfi_fit_use == 'three groups':  
        # guess rfi set 3 besides set 1 & 2 
        peak_use = np.full(len(peaks),False)
        peak_use[fcenter3_loc] = True
        peaks[~peak_use] = 0

        if (peaks == 0.).all() == False:
            fc02 = fc0[np.argmax(peaks)]

            former = int((fc02 - freq[0])//freq_step + 10)
            fc3 = np.arange(fc02-freq_step*former,1440,freq_step)  
            peak3 = []
            for j in range(len(fc0)):
                vc = fc0[j]
                if np.min(np.abs(vc - fc3)) < freq_thr:
                    fcenter3 += [vc,]; peak3 +=[peaks[j]]
            fcenter3 = np.array(fcenter3)
            if check:
                if fcenter3.size > 0:
                    peak3 = np.array(peak3)       
                    _,fcenter3 = check_repeat_center(fcenter3,freq_step=8.1,x_start = '1st point',peaks = peak3)
    
    fcenter3 = np.array(fcenter3)
    fc0 = np.array(fc0).flatten()
    
    return fcenter1,fcenter2,fcenter3,fc0

def check_repeat_center(item,freq_step=8.1,x_start = '1st point',peaks = None):
    if len(item) <= 0:
        raise ValueError(f"Input item length is {len(item)}!")
    if x_start == 'zero':
        x = np.around(item/freq_step)
    elif x_start == '1st point':
        x = np.around((item-item[0])/freq_step)

    if len(set(x)) < len(x):
        from collections import Counter 
        c = Counter(x)
        repeat_loc = np.array(list(c.values())) > 1
        repeat_key = np.array(list(c.keys()))[repeat_loc]

        if isinstance(peaks,np.ndarray):
            delete_loc = []
            for i in range(len(repeat_key)):
                repeat_loc_in_x = np.where(x == repeat_key[i])[0]
                peak = np.full(len(x),-1)
                peak[repeat_loc_in_x] = peaks[repeat_loc_in_x]
                Max = np.max(peak)
                rep = np.sum(peak == Max)
                if rep > 1:
                    repeat_peak = np.where(peak == np.max(peak))[0]
                    del_loc = np.hstack((repeat_peak[:-1],
                                        np.where((peak >= 0)&(peak < np.max(peak)))[0]))
                else:
                    del_loc = np.where((peak >= 0)&(peak < Max))[0]
                delete_loc.append(del_loc)
            delete_loc = np.hstack(delete_loc)
        else:
            repeat_loc_in_x = np.array([np.where(x == repeat_key[i])[0] for i in range(len(repeat_key))])
            delete_loc = repeat_loc_in_x[:,1:].flatten()
            
        x = np.delete(x,delete_loc,)
        item = np.delete(item,delete_loc,)
        #log.warning("RFIs are too closed. So x has repeated item, delete {len(delete_loc)}.")
        #raise ValueError(f" RFIs are too closed. So x has repeated item, delete {len(delete_loc)}.")
    return x,item


def polyfit1order(item,plot = False,pdf = None,**kwargs):
    x,item = check_repeat_center(item,**kwargs)
        
    if len(x) <= 1:
        raise ValueError(f" x only has {len(x)} item, can't polyfit.")
        
    pfit = np.polyfit(x,item,1)
    pfunc = np.poly1d(pfit)
    if plot:
               
        if pdf is not None:
            plt.switch_backend('agg')
        fig,ax = plt.subplots()
        ax.plot(x,pfunc(x),label = f'k = {pfunc[1]:.3f}')
        ax.plot(x,item,'x',label = f'b = {pfunc[0]:.2f}')
        ax.legend()
        ax.grid()
        if pdf is not None:
            pdf.savefig();plt.close()
    return pfunc

def mask_RFI(freq,is_rfi,theory,mask_all_theory = False,freq_from_theory = None,**kwargs):
    # fix width
    rfi_width_lim = kwargs['rfi_width_lim']
    ext_sec = kwargs['ext_sec']
    
    start,end = get_startend(is_rfi,rfi_width_lim=rfi_width_lim/2,ext_sec= ext_sec)
    bigrfi = np.full(is_rfi.shape[0],False)
    for s,e in zip(start,end):
        fs,fe = freq[s],freq[e]
        mark = (theory > fs) & (theory < fe)
        if mark.any() == True:
            bigrfi[s:e] = True
    
    if mask_all_theory:
        bigrfi2 = np.zeros_like(bigrfi,dtype = 'bool')
        for th in theory:
            bigrfi2[(freq>th - freq_from_theory/2) & (freq<th + freq_from_theory/2)] = True
        bigrfi3 = bigrfi | bigrfi2
    else:
        bigrfi3 = bigrfi
    
    return bigrfi,bigrfi3

def find_edge_2sides(spec,vel,peak_position,step,rms_thresh,Print=False,small_rfi_times=2):
    start = deepcopy(peak_position)
    part_rms = rms(spec,vel,rms_vrange=[start-step/2,start+step/2])
    
    if part_rms < rms_thresh * small_rfi_times:
        #print("theory peak rms < rms_thresh")
        edge = np.array([np.nan,np.nan])
        return edge
    
    while part_rms > rms_thresh:
        start += -step
        if start < min(vel):
            break
        part_rms = rms(spec,vel,rms_vrange=[start-step/2,start+step/2])           
    
    end = deepcopy(peak_position)
    part_rms = rms(spec,vel,rms_vrange=[end-step/2,start+end/2])
    while part_rms > rms_thresh:
        end += step
        if end > max(vel):
            break
        part_rms = rms(spec,vel,rms_vrange=[end-step/2,end+step/2]) 
      
    edge = np.array([start,end]) #+ np.array([-1,1]) * ext_times * step
    edge[edge < min(vel)] = min(vel)
    edge[edge > max(vel)] = max(vel)
    if Print:
        print(start,end)  
        print(peak_position,edge)
    return edge

def mask_RFI_2sides(spec,freq,is_rfi,theory0,mark_rfi_width,small_rfi_times,RMS,
                    ext_times = 0,chan_step=3,**kwargs):
    flim_l = theory0 - mark_rfi_width/2
    flim_r = theory0 + mark_rfi_width/2
    fdelta = freq[1] - freq[0]
    
    rfi_width_lim = kwargs['rfi_width_lim']
    ext_sec = kwargs['ext_sec']
    
    start,end = get_startend(is_rfi,rfi_width_lim=rfi_width_lim/2,ext_sec= ext_sec)
    '''
    bigrfi = np.full(is_rfi.shape[0],False)
    for s,e in zip(start,end):
        fs,fe = freq[s],freq[e]
        mark = (theory0 > fs) & (theory0 < fe)
        if mark.any() == True:
            bigrfi[s:e] = True
    '''       
    bigrfi2 = np.full(is_rfi.shape[0],False)
    for j in range(len(flim_l)):
        rfi_part = (freq > flim_l[j])&(freq < flim_r[j])
        rspec = deepcopy(spec)[rfi_part]
        
        edge = find_edge_2sides(rspec,freq[rfi_part],peak_position = theory0[j],small_rfi_times = small_rfi_times,
                         step=chan_step*fdelta,rms_thresh=RMS,Print=False)
        if edge[0] < edge[1]:
            edge_use = (freq >= edge[0])&(freq <= edge[1])
            bigrfi2[edge_use] = True
    bigrfi = bigrfi2    
        
    return bigrfi,bigrfi2       
            
def center_theory(freq,fcenter1,fcenter2,fcenter3,freq_step,plot,pdf,rfi_fit_use = 'two groups'):
    if fcenter3.size <= 1:
        rfi_fit_use = 'two groups'
    
    if (rfi_fit_use == 'two groups')|(rfi_fit_use == 'three groups'):
        pfunc1 = polyfit1order(fcenter1, freq_step=freq_step, plot = plot, pdf = pdf)
        pfunc2 = polyfit1order(fcenter2, freq_step=freq_step, plot = plot, pdf = pdf)
        
        former = int((pfunc1[0] - freq[0])//pfunc1[1] + 10)
        theory1 = np.arange(pfunc1[0]-pfunc1[1]*former,1460,pfunc1[1])
        theory1 = np.delete(theory1,np.where((theory1 < freq[0])|(theory1 > freq[-1]))[0])
        former = int((pfunc2[0] - freq[0])//pfunc2[1] + 10)
        theory2 = np.arange(pfunc2[0]-pfunc2[1]*former,1460,pfunc2[1])
        theory2 = np.delete(theory2,np.where((theory2 < freq[0])|(theory2 > freq[-1]))[0])
        
        if rfi_fit_use == 'three groups':  
            try:
                pfunc3 = polyfit1order(fcenter3, freq_step=freq_step, plot = plot, pdf = pdf)
                former = int((pfunc3[0] - freq[0])//pfunc3[1] + 10)
                theory3 = np.arange(pfunc3[0]-pfunc3[1]*former,1460,pfunc3[1])
                theory3 = np.delete(theory3,np.where((theory3 < freq[0])|(theory3 > freq[-1]))[0])
            except ValueError:
                theory3 = np.array([])
            theory0 = np.hstack((theory1,theory2,theory3)) 
        elif rfi_fit_use == 'two groups':      
            theory3 = np.array([])
            theory0 = np.hstack((theory1,theory2))
        theory0.sort()  
        
    elif rfi_fit_use == 'all':
        former = int((pfunc0[0] - freq[0])//pfunc0[1] + 10)
        pfunc0 = polyfit1order(fc0, freq_step=freq_step/2)
        theory0 = np.arange(pfunc0[0]-pfunc0[1]*4*former,1460,pfunc0[1])
        theory0 = np.delete(theory0,np.where((theory0 < freq[0])|(theory0 > freq[-1]))[0])
        theory1,theory2,theory3 = [],[],[]
    
    return theory0,theory1,theory2,theory3
        
    
    
def find_RFI(spec,freq,is_rfi,freq_step=8.1,RMS = None,freq_thr = 0.5, ext_edge = 0,
             plot = False,pdf = None,ylim = None,rfi_fit_use = 'two groups',
             mask_RFI_method = 'fixed freq',small_rfi_times = 2,chan_step = 3,
             mask_all_theory = False,freq_from_theory = None,mask_thr = 15,
             ret_type = 'all theory',**kwargs):
         
    if plot:
             
        if pdf is not None:
            plt.switch_backend('agg')

    fdelta = freq[1] - freq[0]

    fcenter1,fcenter2,fcenter3,fc0 = find_center(spec,freq,is_rfi,RMS,freq_step=freq_step,
                                                 freq_thr =freq_thr,rfi_fit_use = rfi_fit_use,**kwargs)
    theory0,theory1,theory2,theory3 = center_theory(freq,fcenter1,fcenter2,fcenter3,freq_step =freq_step,
                            plot = plot,pdf = pdf,rfi_fit_use =rfi_fit_use)

    is_rfi_ = (spec > RMS * mask_thr)
    if mask_RFI_method == '2 sides':
        bigrfi,bigrfi2 = mask_RFI_2sides(spec,freq,is_rfi,theory0,mark_rfi_width = freq_from_theory,
                                 small_rfi_times = small_rfi_times,chan_step = chan_step,RMS = RMS,**kwargs)
    elif mask_RFI_method == 'fixed freq':
        bigrfi,bigrfi2 = mask_RFI(freq,is_rfi_,theory0,mask_all_theory,freq_from_theory,**kwargs)
        
    if ext_edge > 0:
        from hifast.utils.misc import extend_Trues
        bigrfi2 = extend_Trues(bigrfi2,axis = -1,ext_add = ext_edge)
     
    if plot:   
        fig,ax = plt.subplots(figsize=(15,3))
        ax.plot(freq,spec)
        spec_ = np.full(spec.shape[0],np.nan)
        spec_[bigrfi] = spec[bigrfi]
        ax.plot(freq,spec_)
        
        ymin = np.min(spec);ymax = np.max(spec)*.8
        
        ax.vlines(fc0,ymin=ymin,ymax=ymax,linestyles='--',colors='k',label = 'all')
        ax.vlines(fcenter1,ymin=ymin,ymax=ymax,linestyles='--',colors='r',label = '1')
        ax.vlines(fcenter2,ymin=ymin,ymax=ymax,linestyles='--',colors='g',label = '2')
        if (rfi_fit_use == 'three groups') & (fcenter3.size>0):
            ax.vlines(fcenter3,ymin=ymin,ymax=ymax,linestyles='--',colors='b',label = '3')
        if ylim is not None:
            ax.set_ylim(ylim[0],ylim[1])
        ax.set_title('RFI positions coarse measure')
        plt.grid();plt.legend()
        if pdf is not None:
            pdf.savefig();plt.close()
            
        fig,ax = plt.subplots(figsize=(15,3))    
        ax.plot(freq,spec)
        spec_ = np.full(spec.shape[0],np.nan)
        spec_[bigrfi2] = spec[bigrfi2]
        ax.plot(freq,spec_)
        if (rfi_fit_use == 'two groups')|(rfi_fit_use == 'three groups'):
            ax.vlines(theory1,ymin=ymin,ymax=ymax,linestyles='--',colors='r',label = 'theory1')
            ax.vlines(theory2,ymin=ymin,ymax=ymax,linestyles='--',colors='g',label = 'theory2')
            if (rfi_fit_use == 'three groups') & (theory3.size>0):
                ax.vlines(theory3,ymin=ymin,ymax=ymax,linestyles='--',colors='b',label = 'theory3')
        elif rfi_fit_use == 'all':
            ax.vlines(theory0,ymin=ymin,ymax=ymax,linestyles='--',colors='k',label = 'theory0')
        if ylim is not None:
            ax.set_ylim(ylim[0],ylim[1])
        ax.set_title('RFI positions prediction')
        plt.grid();plt.legend()
        if pdf is not None:
            pdf.savefig();plt.close()
            
    if ret_type == 'all theory':
        return bigrfi2,theory0
    elif ret_type == 'theory 123':
        return bigrfi2,(theory1,theory2,theory3)





####################### time rfi ########################

def find_t(data,freq,frange,times = 10,thr = 20,rfi_width_lim = 20,
           ext_add = 0,plot = False,pdf = None,ylim = None):
    t = np.arange(data.shape[0])
    is_timerfi = np.zeros_like(t,dtype = 'bool')
    
    if (freq is not None)&(len(frange) == 2):
        pattern_use = (freq>frange[0])&(freq<frange[1])
        
        pat_data = data[:,pattern_use]
    else:
        pat_data = data

    pat_mean = np.nanmean(pat_data,axis = 1)
    pat_med = np.nanmedian(pat_mean)
    
    pat_mean[np.isnan(pat_mean)] = 0
    
    is_pat = pat_mean > pat_med * times
    if is_pat.any() == False:
        if plot:  
            if pdf is not None:
                plt.switch_backend('agg')
            fig,ax = plt.subplots(figsize = (15,3))
            ax.plot(t,pat_mean)
            ax.axhline(pat_med*times,c = 'k',linestyle = '--',label='median')
            ax.grid()
            ax.set_ylabel('mean along freq axis')
            ax.set_xlabel('spec number (time)')
            ax.set_title(f'freq in {frange} MHz')
            if ylim is not None:
                ax.set_ylim(ylim[0],ylim[1])
            if pdf is not None:
                pdf.savefig();plt.close()
        return  is_timerfi
    
    #from markRFI import get_startend
    try:
        start,end = get_startend(is_pat,rfi_width_lim = rfi_width_lim,ext_sec=0)
    except ValueError :
        import traceback
        traceback.print_exc() 
        return  is_timerfi
       
    pat_diff = np.abs(np.diff(pat_mean,append = 0))

    for s,e in zip(start,end):
        #cond = ((pat_diff[s] - pat_diff[s-1])/pat_med > thr) & ((pat_diff[e] - pat_diff[e+1])/pat_med > thr)
        cond = (pat_diff[s]/pat_med > thr)& (pat_diff[e]/pat_med > thr)
        if cond:
            is_timerfi[s:e] = True
            
    if ext_add > 0:
        from hifast.utils.misc import extend_Trues
        is_timerfi = extend_Trues(is_timerfi,axis = 0,ext_add = ext_add)
    
    if plot:  
        if pdf is not None:
            plt.switch_backend('agg')
        fig,ax = plt.subplots(figsize = (15,3))
        ax.plot(t,pat_mean)
        ax.axhline(pat_med*times,c = 'k',linestyle = '--',label='median')
        ax.plot(t,is_timerfi*np.max(pat_mean),label='is_timeRFI')
        ax.plot(t,pat_diff-np.max(pat_mean),label='abs(diff)')
        ax.grid()
        ax.plot(t[start],pat_mean[start],'r.',label = 'start')
        ax.plot(t[end],pat_mean[end],'g.',label = 'end')
        ax.legend()
        ax.set_ylabel('mean along freq axis')
        ax.set_xlabel('spec number (time)')
        ax.set_title(f'freq in {frange} MHz')
        if ylim is not None:
            ax.set_ylim(ylim[0],ylim[1])
        if pdf is not None:
            pdf.savefig();plt.close()
            
    return  is_timerfi

def mask_sf(data,freq,T_thr_times = None,frange = None,RMS = None,**kwargs):
    ret = np.zeros_like(data,dtype = 'bool')
    
    log.info(f"Looking for short-freq time RFI in {frange} ...")
    f_use = (freq>frange[0])&(freq<frange[1])
    
    is_timerfi = find_t(data,freq,frange =frange, **kwargs)
    
    if is_timerfi.any() == False:
        return ret
    use1 = np.zeros_like(data,dtype = 'bool')
    use2 = np.zeros_like(data,dtype = 'bool')
    use1[is_timerfi,:]  = True
    use2[:,f_use]  = True
    use = use1 & use2

    tmp = deepcopy(data[use])
    T_thr = RMS * T_thr_times
    ret[use] = (tmp > T_thr)
    
    ext_add = kwargs['ext_add']
    if ext_add is not None:
        from hifast.utils.misc import extend_Trues
        ret = extend_Trues(ret,axis = -1,ext_add = ext_add)

    print("Found :D")
    return ret

def mask_time_rfi(data,freq,T_thr_times = None,rtype = 'short-freq',frange = None,file = None,RMS = None,**kwargs):
    ret = np.zeros_like(data,dtype = 'bool')

    if rtype == 'short-freq':
        if isinstance(file,str):
            franges = np.load(file)
            if franges.shape[1] != 2:
                raise ValueError("time RFI freq shape must like (n,2)")

            for nf in range(franges.shape[0]):
                ret = ret | mask_sf(data,freq,T_thr_times,frange = franges[nf],RMS = RMS,**kwargs)

        elif len(frange) == 2:
            ret = mask_sf(data,freq,T_thr_times,frange,RMS = RMS,**kwargs)
        else:
            raise ValueError("frange should like [fmin,fmax] or a string (RFI npy filepath).")

        if ret.any() == False:
            print("No short-freq time rfi is found.")
            
    elif rtype == 'long-freq':
        log.info(f"Looking for long-freq time RFI ...")
        if len(frange) == 2:
            is_timerfi = find_t(data,freq,frange =frange, **kwargs)
            if is_timerfi.any() == False:
                print("No long-freq time rfi is found.")
                return ret
            ret[is_timerfi,:] = True
            print("Found :D")
        else:
            raise ValueError("frange should like [fmin,fmax].")
    
    print("Finish")
    return ret


