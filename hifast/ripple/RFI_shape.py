#!/usr/bin/env python
# coding: utf-8

# author: Xu Chen, 2021.07
# code：Xu Chen

import numpy as np
from copy import deepcopy
from matplotlib import pyplot as plt 

def rfis_in_one_spec(spec,freq,is_rfi,mw_use,freq_step=8.1,RMS = None,freq_thr = 0.5,
                     rfi_fit_use = 'two groups',step=0,ylim = None,plot = False,
                     pdf = None,**kwargs): 

    is_rfi_ = deepcopy(is_rfi)
    is_rfi_[mw_use] = False
    syn_rfi1 = [];syn_rfi2 = [];syn_rfi3 = []
    syn_freq1 = [];syn_freq2 = [];syn_freq3 = []
    
    from markRFI import find_center,find_local_peak,get_startend
    start,end = get_startend(is_rfi_,**kwargs)
    fcenter1,fcenter2,fcenter3,fc0 = find_center(spec,freq,is_rfi_,RMS,freq_step=freq_step,
                                    freq_thr =freq_thr,rfi_fit_use = rfi_fit_use,**kwargs)
    if plot:
        
        if pdf is not None:
            plt.switch_backend('agg')

        fig2,ax = plt.subplots(figsize=(15,3))
        ax.plot(freq,spec)
        ymin = np.min(spec);ymax = np.max(spec)*.8
        ax.vlines(fc0,ymin=ymin,ymax=ymax,linestyles='--',colors='k',label = 'all')
        ax.vlines(fcenter1,ymin=ymin,ymax=ymax,linestyles='--',colors='r',label = '1')
        ax.vlines(fcenter2,ymin=ymin,ymax=ymax,linestyles='--',colors='g',label = '2')
        if (rfi_fit_use == 'three groups') & (fcenter3.size>0):
            ax.vlines(fcenter3,ymin=ymin,ymax=ymax,linestyles='--',colors='b',label = '3')
        ax.grid();ax.legend();
        if ylim is not None:
            ax.set_ylim(ylim[0],ylim[1])
        if pdf is not None:
            pdf.savefig();plt.close()

        fig1 = plt.figure(figsize=(10,4))
        ax1 = fig1.add_subplot(131)
        ax2 = fig1.add_subplot(132)
        ax3 = fig1.add_subplot(133)
        ax1.grid();ax2.grid();ax3.grid()
    
    k1,k2,k3 = 0,0,0
    
    for s,e in zip(start,end):
        rfi1 = spec[s:e]
        freq1 = freq[s:e]
        freq_peak,flux_peak = find_local_peak(rfi1,freq1,RMS=RMS,)

        vc = freq_peak
        if (vc > freq1[0])&(vc < freq1[-1]):
            norm_factor = 1/(np.sum(rfi1))
            nrfi = rfi1 * norm_factor
            nfreq = freq1 - vc

            if vc in fcenter1:
                if plot:
                    ax1.step(nfreq,nrfi+k1,label = f'{vc:.3f}')
                syn_rfi1.append(nrfi)
                syn_freq1.append(nfreq)
                k1 += step
            elif vc in fcenter2:
                if plot:
                    ax2.step(nfreq,nrfi+k2,label = f'{vc:.3f}')
                syn_rfi2.append(nrfi)
                syn_freq2.append(nfreq)
                k2 += step
            elif vc in fcenter3:
                if plot:
                    ax3.step(nfreq,nrfi+k3,label = f'{vc:.3f}')
                syn_rfi3.append(nrfi)
                syn_freq3.append(nfreq)
                k3 += step
    if plot:
        #ax.legend(loc=2,borderaxespad=1,bbox_to_anchor=(1.0,1.0))
        ax1.set_xlabel('MHz');ax2.set_xlabel('MHz');ax3.set_xlabel('MHz');
        ax1.set_ylabel('flux');ax2.set_ylabel('flux');ax3.set_ylabel('flux')
        ax1.set_title('RFI1 in one spec in a beam');
        ax2.set_title('RFI2 in one spec in a beam');
        ax3.set_title('RFI3 in one spec in a beam');
        if pdf is not None:
            pdf.savefig();plt.close()

    return (syn_rfi1,syn_rfi2,syn_rfi3),(syn_freq1,syn_freq2,syn_freq3)


def synthetic_rfi(syn_rfi,syn_freq,fdelta,plot = False,pdf = None,):
    
    neg_num = [np.sum(syn_freq[fn] < 0) for fn in range(len(syn_freq))]
    pos_num = [len(syn_freq[fn]) - neg_num[fn] for fn in range(len(syn_freq))] 
    
    m = np.median(neg_num) 
    n = np.median(pos_num) 
    x_inpd = np.arange(-fdelta*m,fdelta*n,fdelta)
    from scipy.interpolate import interp1d
    syn_rfi_inpd = np.zeros((len(syn_rfi),len(x_inpd)))
    
    if plot:
        if pdf is not None:
            plt.switch_backend('agg')
        fig,ax = plt.subplots(figsize=(5,4))
        
    for j in range(len(syn_rfi)):
        y_interp = interp1d(syn_freq[j],syn_rfi[j],kind='linear',fill_value="extrapolate")
        y_inpd = y_interp(x_inpd)

        syn_rfi_inpd[j,:] = y_inpd
        if plot:
            ax.step(x_inpd,y_inpd)

    syn_rfi_med = np.median(syn_rfi_inpd,axis = 0)
    syn_peak = np.max(syn_rfi_med)
    if plot:
        ax.plot(x_inpd,syn_rfi_med,'k',label = 'interplate median')
        ax.set_title('RFIs after interpolation');
        ax.legend()
        if pdf is not None:
            pdf.savefig();plt.close()

    return x_inpd,syn_rfi_inpd,syn_rfi_med, syn_peak, (m,n)

def common_shape(spec,ori_spec,freq,pd_rfi,theorys,syn_rfis,syn_freqs,m31_use,RMS,
                 freq_thr=.3,plot = False,pdf = None,only_M31 = False,ylim = None,**kwargs):

    fdelta = freq[1]-freq[0]
    #############
    theory1,theory2,theory3 = theorys
    syn_rfi1,syn_rfi2,syn_rfi3 = syn_rfis
    syn_freq1,syn_freq2,syn_freq3 = syn_freqs
    
    x_inpd1,_,syn_rfi_med1, syn_peak1,(m1,n1) = synthetic_rfi(syn_rfi1,syn_freq1,fdelta,plot,pdf)
    x_inpd2,_,syn_rfi_med2, syn_peak2,(m2,n2) = synthetic_rfi(syn_rfi2,syn_freq2,fdelta,plot,pdf)
    if len(syn_rfi3) > 0:
        x_inpd3,_,syn_rfi_med3, syn_peak3,(m3,n3) = synthetic_rfi(syn_rfi3,syn_freq3,fdelta,plot,pdf)

    #################
    pd_rfi_ = deepcopy(pd_rfi)
    pd_rfi_[m31_use] = False
    from markRFI import get_startend,find_local_peak
    start,end = get_startend(pd_rfi_,**kwargs)
    
    def replace_rfi(cfrq,peak,fdelta,m,n,com_rfi,x_inpd,syn_rfi_med,syn_peak):
        fl,fr = 0 - fdelta * m, 0 + fdelta * n
        if fl < freq[0] - cfrq: fl = freq[0] - cfrq 
        if fr > freq[-1] - cfrq: fr = freq[-1] - cfrq 
        rfi_use = (freq >= cfrq - fdelta * m)&(freq < cfrq + fdelta * n)
        syn_use = (x_inpd >= fl)&(x_inpd <= fr)
        syn_rfi_med_ = syn_rfi_med[syn_use]
        
        rfi_len = np.sum(rfi_use)
        syn_len = len(syn_rfi_med_)
        diflen = syn_len - rfi_len
        if diflen > 0:
            syn_rfi_med_ = syn_rfi_med_[:-int(diflen)]
        elif diflen < 0:
            rfi_use[:-int(diflen)] = False
            
        com_rfi[rfi_use] = syn_rfi_med_ * peak / syn_peak
        return com_rfi
    
    com_rfi = np.full_like(freq,np.nan)
    for s,e in zip(start,end):
        rfi1 = spec[s:e]
        freq1 = freq[s:e]
        cfrq,flux_peak = find_local_peak(rfi1,freq1,RMS=RMS,)
        peak = flux_peak + RMS
        if np.min(np.abs(cfrq - theory1)) < freq_thr:
            com_rfi = replace_rfi(cfrq,peak,fdelta,m1,n1,com_rfi,x_inpd1,syn_rfi_med1,syn_peak1)
            
        elif np.min(np.abs(cfrq - theory2)) < freq_thr:
            com_rfi = replace_rfi(cfrq,peak,fdelta,m2,n2,com_rfi,x_inpd2,syn_rfi_med2,syn_peak2)
        else:
            if theory3.size > 0 & len(syn_rfi3) > 0:
                if np.min(np.abs(cfrq - theory3)) < freq_thr:
                    com_rfi = replace_rfi(cfrq,peak,fdelta,m3,n3,com_rfi,x_inpd3,syn_rfi_med3,syn_peak3)
                    
    ret = ori_spec - com_rfi
    ret[np.isnan(ret)] = spec[np.isnan(ret)]
    
    if only_M31:
        m31_unuse = (freq<m31_frange[0]-3)|(freq>m31_frange[1]+3)
        ret[m31_unuse] = spec[m31_unuse]
    
    if plot:
        if pdf is not None:
            plt.switch_backend('agg')
        fig,ax = plt.subplots(figsize=(15,3))
        ax.plot(freq,spec,label='ori')
        ax.grid()
        if ylim is not None:
            plt.ylim(ylim[0],ylim[1])
        ax.plot(freq,com_rfi,label = 'common rfi')
        ax.plot(freq,ret,label = 'diff')
        ax.legend()

        fig,ax = plt.subplots(figsize=(15,3))
        ax.plot(freq,spec,label='ori')
        ax.grid()
        if ylim is not None:
            plt.ylim(ylim[0],ylim[1])
        ax.plot(freq,com_rfi,label = 'common rfi')
        ax.plot(freq,ret,label = 'diff')
        ax.legend();plt.xlim(1400,1440)
        if pdf is not None:
            pdf.savefig();plt.close()
    
    return ret,com_rfi
