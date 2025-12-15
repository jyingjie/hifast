

__all__ = ['PolarMjdChan_to_MjdChanPolar', 'MjdChanPolar_to_PolarMjdChan', 'Test', 'phrase_ylim', 'main']


import numpy as np
from scipy import ndimage
import matplotlib.pyplot as plt
import h5py

import re
import os
from glob import glob
import json

from ..core.baseline import get_baseline, sub_baseline
from ..utils import io
from .widgets import *

import warnings
warnings.filterwarnings("ignore", r'overflow encountered in exp')

import matplotlib as mpl

PolarMjdChan_to_MjdChanPolar = lambda x:io.PolarMjdChan_to_MjdChanPolar(x, check_polar=False)
MjdChanPolar_to_PolarMjdChan = lambda x:io.MjdChanPolar_to_PolarMjdChan(x, check_polar=False)


class Test(object):
    def __init__(self, T2p, freq, frange=None, is_excluded=None, trans=False, save=False):
        """
        T2p: shape (Polar,Mjd,Chan)

        """
        self.trans = trans
        if T2p.ndim != 3:
            raise(ValueError('T2p should be 3 dim'))
        self.T2p = T2p

        self.is_excluded = is_excluded

        freq = freq[:]
        if self.trans:
            freq = np.arange(T2p.shape[1])
        if frange is not None:
            self.is_use = (freq > frange[0]) & (freq < frange[1])
            self.freq = freq[self.is_use]
        else:
            self.is_use = None
            self.freq = freq
        self.select()

        if save:
            # increase memory usage
            # (Mjd,Chan,Polar)
            self.bld_out = np.full((self.T2p.shape[1], len(self.freq), self.T2p.shape[0]), np.nan)
        else:
            self.bld_out = None
        self.save = save

    def select(self, start=0, length=20, polar=0):
        """
        select one polar of spectra from `start` to `start + length`
        """
        if self.trans:
            self.T2p_t = self.T2p[polar:polar+1,:, start:start+length]
            self.T2p_t = PolarMjdChan_to_MjdChanPolar(self.T2p_t)
            self.T2p_t = self.T2p_t.transpose((1, 0, 2))
        else:
            self.T2p_t = PolarMjdChan_to_MjdChanPolar(self.T2p[polar:polar+1, start:start+length])
        # also change self.is_excluded
        if self.is_excluded is not None:
            if self.is_excluded.shape == self.T2p.shape:
                if self.trans:
                    self.is_excluded_t = self.is_excluded[polar:polar+1, :, start:start+length]
                    self.is_excluded_t = PolarMjdChan_to_MjdChanPolar(self.is_excluded_t).transpose((1, 0, 2))
                else:
                    self.is_excluded_t = PolarMjdChan_to_MjdChanPolar(self.is_excluded[polar:polar+1, start:start+length])
            elif self.is_excluded.ndim == 2 and self.is_excluded.shape == T2p.shape[1:]:
                if self.trans:
                    self.is_excluded_t = self.is_excluded[:, start:start+length][None,:]
                    self.is_excluded_t = PolarMjdChan_to_MjdChanPolar(self.is_excluded_t).transpose((1, 0, 2))
                else:
                    self.is_excluded_t = self.is_excluded[start:start+length][..., None]
            else:
                raise(ValueError('shape of ``is_excluded``'))
        else:
            self.is_excluded_t = None

        if self.is_use is not None:
            self.T2p_t = self.T2p_t[:, self.is_use]
            if self.is_excluded_t is not None:
                self.is_excluded_t = self.is_excluded_t[:, self.is_use]

    def sub(self, x, start=0, length=20, polar=0, frange_excluded=None,
            exclude_add='none', verbose=False, **kwargs):

        #store
        self.start = start
        self.length = length
        self.polar = polar

        self.select(start, length, polar)

        if frange_excluded is not None and frange_excluded[0] != frange_excluded[1]:
            is_ = (x > frange_excluded[0]) & (x < frange_excluded[1])
            # full
            is_ = np.full(self.T2p_t.shape, is_[None,...,None])
            if self.is_excluded_t is not None:
                self.is_excluded_t &= is_
            else:
                self.is_excluded_t = is_

        self.frange_excluded = frange_excluded
        self.sub_baseline_para = kwargs
        self.sub_baseline_para['is_excluded'] = self.is_excluded_t
        self.sub_baseline_para['exclude_add'] = exclude_add
        self.sub_baseline_para['knots'] = json.loads(self.sub_baseline_para['knots'])

        self.bld = sub_baseline(self.freq, self.T2p_t, verbose=verbose, **self.sub_baseline_para)
        if self.bld_out is not None:
            self.bld_out[start:start+length, :, polar:polar+1] = self.bld

        return np.full(x.shape, np.nan)

    def get_ori(self, x, i, **kwargs):
        res = self.T2p_t[i, :, 0]
        return np.vstack([res, ndimage.gaussian_filter1d(res, 3)]).T

    def get_bld(self, x, i, **kwargs):
        res = self.bld[i, :, 0]
        return np.vstack([res, ndimage.gaussian_filter1d(res, 3)]).T

    def get_bl(self, x, i,  **kwargs):
        return self.T2p_t[i, :, 0] - self.bld[i, :, 0]

    def get_bld_mean(self, x, start_stop,  **kwargs):
        start, stop = start_stop
        res = np.mean(self.bld[int(start):int(stop)+1, :, 0], axis=0)
        return np.vstack([res, ndimage.gaussian_filter1d(res, 3)]).T


def phrase_ylim(ylim, vals):
    is_ = np.isfinite(vals)
    ylim = list(map(lambda x: float(x) if 'per' != str(x)[:3] else
                    np.nanpercentile(vals[is_], float(str(x)[3:]), interpolation='nearest'), ylim))
    return ylim



def main():
    global ylim

    tes = Test(T2p, freq, frange, is_excluded, trans, save)

    rew_type = ['asym1', 'asym2', 'asym3', 'sym1',]
    method = ['PLS-'+r for r in rew_type]
    method += ['poly-'+r for r in rew_type]
    method += ['Gauss-'+r for r in rew_type]
    method += ['knspline-'+r for r in rew_type]
    method += ['knpoly-'+r for r in rew_type]
    method += ['spline-'+r for r in rew_type]
    method += ['masPLS-'+r for r in rew_type]
    method += ['asPLS',]

    sliders = {}
    if trans:
        start_e = tes.T2p.shape[2]-length
    else:
        start_e = tes.T2p.shape[1]-length
    if start_init is None:
        start_i = start_e//2
    else:
        start_i = min(start_init, start_e)
    sliders.update(_BoundedIntText(
                       start=(start_i, 0, start_e, 1), polar=(0, 0, tes.T2p.shape[0], 1)))
    sliders.update(_IntSlider(njoin=(0, 1, length, 1),
                              average_every_freq=(0, 1, 40, 1)))
    sliders.update(_Dropdown(method=method,
                             s_method_freq=('none', 'gaussian', 'boxcar',),
                             s_method_t=('none', 'gaussian', 'boxcar'),
                             exclude_add=('none', 'auto1', 'auto2'),
                             ))
    sliders.update(_IntSlider(s_sigma_freq=(3, 1, 20, 1), s_sigma_t=(3, 1, 20, 1)))
    sliders.update(_FloatLog10Slider(lam=(1e8, 0, 13, 0.2), readout_format='.2e'))
    sliders.update(_IntSlider(deg=(2, 1, 10, 1),
                              niter=(100, 1, 200, 1)))
    sliders.update(_FloatSlider(offset=(2, 0.1, 4, 0.1)))
    sliders.update(_FloatSlider(ratio=(0.01, 0.001, 0.015,0.001),
                                readout_format='.3f'))

    sliders['knots'] = widgets.Text(
                    value=json.dumps(tes.freq[[1,-2]].tolist()),
                    description='knots',
                    layout=widgets.Layout(width=f"{mpl.rcParams['figure.dpi']*figsize[0]*1.2}px",),
                    disabled=False)


    bak = w_conf.pop('layout')
    #w_conf['style'] = {'description_width': 'initial'}
    w_conf['layout'] = widgets.Layout(width=f"{mpl.rcParams['figure.dpi']*figsize[0]*1.2}px",)
    sliders.update(_FloatRangeSlider(frange_excluded=([tes.freq[0], tes.freq[0]], tes.freq[0], tes.freq[-1], tes.freq[2]-tes.freq[0])))
    w_conf['layout'] = bak

    plt.ioff()
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=figsize)
    fig.canvas.header_visible = False

    # ylim xlim in plot
    if isinstance(ylim, str):
        if ylim == 'auto':
            ylim = phrase_ylim(['per0.1', 'per99.5'], PolarMjdChan_to_MjdChanPolar(T2p[:,:5]))
        elif ylim == 'full':
            ylim = 'fixed'
    else:
        ylim = phrase_ylim(ylim, tes.T2p_t)


    if isinstance(ylim, str):
        ylim2 = ylim
    else:
        ylim2 = [(ylim[0]-ylim[1])/2, (ylim[1]-ylim[0])/2]
    xlim = [tes.freq.min(), tes.freq.max()]

    widget_i = _IntSlider(i=(0, 0, length-1, 1))['i'] if length ==1 else range(length)


    # init
    controls = iplt.plot(tes.freq, tes.sub, nproc=nproc, length=length,
                         **sliders,
                         xlim=xlim, ylim=ylim,
                         color='r', ax=ax1, play_buttons=True, display_controls=False)
    controls2 = iplt.plot(tes.freq, tes.get_ori, i=widget_i,
                          controls=controls, xlim=xlim, ylim=ylim, ax=ax1, display_controls=False, play_buttons=True)
    controls3 = iplt.plot(tes.freq, tes.get_bl,
                          controls=controls2, xlim=xlim, ylim=ylim, ax=ax1)

    iplt.plot(tes.freq, tes.get_bld,
              controls=controls2, xlim=xlim, ylim=ylim2, ax=ax2)
    plt.axhline(y=0,)

    controls4 = iplt.plot(tes.freq, tes.get_bld_mean, start_stop=('r', 0, length, 11),
                          controls=controls, xlim=xlim, ylim=ylim2, ax=ax3)

    plt.axhline(y=0, )
    ax1.set_title('orange line: gaussian smoothed for better view. \n origial spectra', fontsize=8)
    ax2.set_title('spectra after baseline removed', fontsize=8)

    def return_v(*args, **kwargs): return tes.frange_excluded[0]
    [iplt.axvline(return_v, controls=controls, ax=ax, color='r') for ax in [ax1, ax2, ax3]]
    def return_v(*args, **kwargs): return tes.frange_excluded[1]
    [iplt.axvline(return_v, controls=controls, ax=ax, color='r') for ax in [ax1, ax2, ax3]]

    fig.tight_layout()

    w = controls.controls
    hbs = [
        widgets.HBox([widgets.Label(value=rf"$Select$ spectra from $start$(max:{controls.controls['start'].max}) to $start+{length}$ to fit:"),
                      w['start'], w['polar'],
                      widgets.Label('.'),
                      w['ratio']]),
        widgets.HBox(
            [widgets.Label(value=r"$One$ spectrum in top & middle:"),
             w['i'],
             widgets.Label(value=r"$Stacked$ spectra in bottom:"),
             w['start_stop']]),
        widgets.HBox([w['njoin'],
                      w['s_method_t'],
                      w['s_sigma_t'],
                      w['niter'],
                     ]),
        widgets.HBox([w['average_every_freq'],
                      w['s_method_freq'],
                      w['s_sigma_freq'],
                      w['offset']]),
        widgets.HBox([w['method'], w['lam'], w['deg'],
                      w['exclude_add']]),
        widgets.HBox([w['knots'], ]),
        widgets.HBox([w['frange_excluded'], ]),
    ]

    BOX = widgets.VBox(hbs)
    display(BOX)
    plt.show()
    display(BOX)

    return tes, w, BOX
