#!/usr/bin/env python
# coding: utf-8

# Chuan-Peng Zhang cpzhang@nao.cas.cn
# Yingjie Jing

import numpy as np
import scipy.interpolate as interp
from scipy import ndimage

from ..utils.misc import average_every_n, boxcar_smooth1d, mask_Trues, extend_Trues, smooth1d
from scipy import ndimage

def get_mean_rms(freq, T2p, **kwargs):
    rfi_nan = (freq>1390.0) & (freq<1410.0)
    m = np.mean(T2p[:, rfi_nan], axis=1)[:,None]
    std= np.std(T2p[:, rfi_nan], axis=1)[:,None]
    return m, std

def get_rfi_folder(freq, T2p, folder, times, ext_add=0, ext_frac=0):
    rms_times = times
    nth = np.arange(len(T2p))
    T2p_m = average_every_n(T2p, folder, axis=0, drop=False)
    nth_m = average_every_n(nth, folder, axis=-1,drop=False)

    m, std = get_mean_rms(freq, T2p_m)
    is_rfi = (T2p_m > m+std*rms_times) | (T2p_m < m-std*rms_times)
    
    if ext_add > 0 or ext_frac > 0:
        is_rfi = extend_Trues(is_rfi, axis=1, ext_add=ext_add, ext_frac=ext_frac)
    
    inds = interp.interp1d(nth_m, range(len(nth_m)), kind='nearest',fill_value='extrapolate')(nth).astype('int')

    return is_rfi[inds]


def find_signal(T2p, *, s_method_t='gaussian', s_sigma_t=10, s_method_freq='gaussian', s_sigma_freq=5, 
              times_s=2, times=5, ext_add=0, ext_frac=0):
    """
    T2p: array, shape=(x,N,2); T or flux with two polar and sorted by time
    s_sigma: int; gaussian smooth size along time (axis=0)
    chan_smooth_method: method of smoothing along channel axis
    times_s: # for smoothed channel
    times: # for not smoothed channel
    """
    # smooth along time
    if s_method_t is not None: 
        T2p_s = smooth1d(T2p, s_method_t, s_sigma_t, axis=0)
    else:
        T2p_s = T2p

    # smooth channel
    if s_method_freq is not None:
        T2p_sc = smooth1d(T2p_s, s_method_freq, s_sigma_freq, axis=1)
    else:
        T2p_sc = T2p_s

    # estimated rms
    pers = np.nanpercentile(T2p_s, [10, 50, 90], axis=1)
    rms = np.min(abs(np.diff(pers,axis=0)),axis=0)
    # find rfi
    # compare channel smoothed 
    is_exceeded = abs(T2p_sc) >= times_s*rms[:,None]
    # compare origin 
    is_exceeded = is_exceeded | (abs(T2p_s) >= times*rms[:,None])
    
    is_exceeded = is_exceeded[:,:,0] | is_exceeded[:,:,1]
    if ext_add > 0 or ext_frac > 0:
        is_exceeded = extend_Trues(is_exceeded, axis=1, ext_add=ext_add, ext_frac=ext_frac)
    #return is_exceeded
    # n continue along t axis as rfi
    return is_exceeded

def get_rfi_c(T2p, n_continue=50, return_sig=False, **kwargs):
    is_sig = find_signal(T2p, **kwargs)
    is_rfi = mask_Trues(is_sig, axis=0, leng_lim=n_continue)
    if return_sig:
        return is_rfi, is_sig
    return is_rfi

def mask_rfi_t(freq, T2p, method='smooth', **kwargs):
    if method=='smooth':
        return get_rfi_c(T2p, **kwargs)
    elif method=='folder':
        return get_rfi_folder(freq, T2p, **kwargs)
    return is_rfi

## compare "signal" with rfi
def cross_rfi_1d(arr, is_rfi, frac_match=0.3):
    """
    is_rfi:
    frac_match:
    """
    leng_lim = 1
    # add False in the beginning and ending of arr
    arr = np.hstack([[False],arr,[False]])
    diff = np.diff(arr.astype('int16'))

    ind_neg = np.where(diff==-1)[0]
    ind_posi = np.where(diff==1)[0]

    leng = (ind_neg - ind_posi)
    is_use = leng > leng_lim 

    for i,j,_leng in zip(ind_posi[is_use],ind_neg[is_use], leng[is_use]):
        leng_in = np.sum(is_rfi[i:j])
        if leng_in/_leng < frac_match:
            arr[i+1:j+1] = False
    return arr[1:-1] | is_rfi

def cross_rfi_axis1_d2(arr, is_rfi, inplace=False, **kwarg):
    if not inplace:
        out = np.zeros_like(is_rfi)
    else:
        out = is_rfi
    for i in range(len(arr)):
        out[i] = cross_rfi_1d(arr[i], is_rfi[i], **kwarg)
    return out
