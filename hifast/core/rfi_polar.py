#!/usr/bin/env python
# coding: utf-8

import h5py
import numpy as np
import matplotlib.pyplot as plt
#from fast_python.baseline import baseline
from scipy import ndimage
from ..utils.misc import extend_Trues

def mask_rfi_p(T2p, s_sigma=10, chan_smooth_method='gaussian', times_s=1, times=6, ext_add=0, ext_frac=0):
    """
    T2p: array, shape=(x,N,2); T or flux with two polar and sorted by time
    s_sigma: int; gaussian smooth size along time (axis=0)
    chan_smooth_method: method of smoothing along channel axis
    times_s: # for smoothed (XX-YY)
    times: # for origin (XX-YY)
    """
    if s_sigma is not None: 
        T2p_s = ndimage.gaussian_filter1d(T2p, s_sigma, axis=0)
    else:
        T2p_s = T2p

    diff = T2p_s[:,:,0] - T2p_s[:,:,1]
    # smooth channel
    chan_smooth_method = 'gaussian'
    chan_smooth_sigma = 5
    if chan_smooth_method == 'gaussian':
        diff_s = ndimage.gaussian_filter1d(diff, chan_smooth_sigma, axis=1)
    else:
        diff_s = np.zeros_like(diff)
        for i in range(len(diff)):
            diff_s[i]= baseline(x=None, y=diff[i], method = 'arPLS', lam=400, deg=3, offset=2)[1]
    # residual
    resi = diff - diff_s
    vdown = np.nanpercentile(resi, 10, axis=1)
    vmean = np.nanpercentile(resi, 50, axis=1)
    vupper = np.nanpercentile(resi, 90, axis=1)
    # find rfi
    # compare smoothed diff
    is_rfi = diff_s < (vmean - times_s*(vmean - vdown))[:,None]
    is_rfi = is_rfi | (diff_s > (vmean + times_s*(vupper - vmean))[:,None])
    # compare origin diff
    is_rfi = is_rfi | (diff < (vmean - times*(vmean - vdown))[:,None])
    is_rfi = is_rfi | (diff > (vmean + times*(vupper - vmean))[:,None])
    # relatively diff
    diff_frac = diff_s/np.nanmean(T2p_s,axis=2)
    is_rfi_and= abs(diff_frac) > 0.25
    diff_frac_2 = diff/np.nanmean(T2p_s,axis=2)
    is_rfi_and= is_rfi_and | (abs(diff_frac) > 0.25)
    
    is_rfi = is_rfi & is_rfi_and
    
    if ext_add > 0 or ext_frac > 0:
        is_rfi = extend_Trues(is_rfi, axis=1, ext_add=ext_add, ext_frac=ext_frac)
    return is_rfi
