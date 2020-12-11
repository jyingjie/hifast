#!/usr/bin/env python
# coding: utf-8

# Chuan-Peng Zhang cpzhang@nao.cas.cn
# Yingjie Jing

import numpy as np
import scipy.interpolate as interp
from scipy import ndimage

from .util import average_every_n, boxcar_smooth1d, mask_Trues, extend_Trues
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


def get_rfi_c(T2p, n_continue=50, s_sigma=10, chan_smooth_method='gaussian', times_s=2, times=5, ext_add=0, ext_frac=0):
    """
    T2p: array, shape=(x,N,2); T or flux with two polar and sorted by time
    s_sigma: int; gaussian smooth size along time (axis=0)
    chan_smooth_method: method of smoothing along channel axis
    times_s: # for smoothed channel
    times: # for not smoothed channel
    """
    if s_sigma is not None: 
        T2p_s = ndimage.gaussian_filter1d(T2p, s_sigma, axis=0)
    else:
        T2p_s = T2p

    # smooth channel
    chan_smooth_sigma=5
    if chan_smooth_method == 'gaussian':
        T2p_sc = ndimage.gaussian_filter1d(T2p_s, chan_smooth_sigma, axis=1)
    elif chan_smooth_method == 'boxcar':
        T2p_sc = util.boxcar_smooth1d(T2p_s, chan_smooth_sigma, axis=1)

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
    is_rfi = mask_Trues(is_exceeded, axis=0, leng_lim=n_continue)
    return is_rfi

def mask_rfi_t(freq, T2p, method='smooth', **kwargs):
    if method=='smooth':
        is_rfi = get_rfi_c(T2p, **kwargs)
    elif method=='folder':
        is_rfi = get_rfi_folder(freq, T2p, **kwargs)
    return is_rfi
            