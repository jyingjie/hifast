

__all__ = ['PolarMjdChan_to_MjdChanPolar', 'MjdChanPolar_to_PolarMjdChan', 'Test', 'phrase_ylim', 'main']


import numpy as np
from scipy import ndimage
import matplotlib.pyplot as plt
import h5py

import re
import os
from glob import glob

from ..core.baseline import get_baseline, sub_baseline, get_exclude_fun
from ..utils import io
from .widgets import *

import warnings
warnings.filterwarnings("ignore", r'overflow encountered in exp')

import matplotlib as mpl

PolarMjdChan_to_MjdChanPolar = lambda x:io.PolarMjdChan_to_MjdChanPolar(x, check_polar=False)
MjdChanPolar_to_PolarMjdChan = lambda x:io.MjdChanPolar_to_PolarMjdChan(x, check_polar=False)


T2p = None
freq = None
frange = None
nproc = 10

length = 20
figsize = (10, 7)
ylim = 'auto'

trans = False
start_init = None

save = True



from .bld_i import Test as Test_bld

class Test(Test_bld):


    def sub(self, x, start=0, length=20, polar=0, frange_excluded=None,
            exclude_add='none', verbose=False,
            bound_f=None, exclude_m=0, **kwargs):
        self.select(start, length, polar)
        if frange_excluded is not None and frange_excluded[0] != frange_excluded[1]:
            is_ = (x > frange_excluded[0]) & (x < frange_excluded[1])
            self.is_excluded_t = np.full(self.T2p_t.shape, is_[None,...,None])
        else:
            self.is_excluded_t = None
        self.frange_excluded = frange_excluded

        self.sub_baseline_para = kwargs
        self.sub_baseline_para['is_excluded'] = self.is_excluded_t
        self.sub_baseline_para['exclude_add'] = exclude_add
        self.sub_baseline_para['method'] = 'sin_poly'
        self.sub_baseline_para['niter'] = 1
        if exclude_m >= 0:
            self.sub_baseline_para['exclude_fun'] = get_exclude_fun(exclude_m)
        else:
            self.sub_baseline_para['exclude_fun'] = None

        bounds = [(0., 1.), bound_f, (0, 2*np.pi), (-1, 1)] + [(-np.inf, np.inf), ]*kwargs['deg']  # optimize.minimize
        self.sub_baseline_para['opt_para'] = {'bounds': bounds, }
        self.bld = sub_baseline(self.freq, self.T2p_t, verbose=verbose, **self.sub_baseline_para)
        if self.bld_out is not None:
            self.bld_out[start:start+length, :, polar:polar+1] = self.bld
        return np.full(x.shape, np.nan)


def phrase_ylim(ylim, vals):

    ylim = list(map(lambda x: float(x) if 'per' != str(x)[:3] else
                    np.nanpercentile(vals, float(str(x)[3:]), interpolation='nearest'), ylim))
    return ylim


def main():
    global ylim

    tes = Test(T2p, freq, frange, save=save)

    sliders = {}

    start_e = tes.T2p.shape[1]-length
    if start_init is None:
        start_i = start_e//2
    else:
        start_i = min(start_init, start_e)
    sliders.update(_BoundedIntText(
                       start=(start_i, 0, start_e, 1), polar=(0, 0, tes.T2p.shape[0], 1)))
    sliders.update(_IntSlider(njoin=(0, 1, length, 1),
                              average_every_freq=(0, 1, 40, 1)))
    sliders.update(_Dropdown(s_method_freq=('none', 'gaussian', 'boxcar'),
                             s_method_t=('none', 'gaussian', 'boxcar'),
                             ))
    sliders.update(_IntSlider(s_sigma_freq=(3, 1, 20, 1), s_sigma_t=(3, 1, 20, 1)))
    sliders.update(_FloatSlider(sin_f=(0.929, 0.910, 0.940, 0.002), readout_format='.3f'))
    sliders.update(_IntSlider(deg=(0, 0, 5, 1),
                              exclude_m=(0, -1, 1, 1)))
    sliders.update(_FloatRangeSlider(bound_f=([.90, .95], 0.8, 1.15, 0.02), readout_format='.3f',))

    bak = w_conf.pop('layout')
    #w_conf['style'] = {'description_width': 'initial'}
    w_conf['layout'] = widgets.Layout(width=f"{mpl.rcParams['figure.dpi']*figsize[0]*1.2}px",)
    sliders.update(_FloatRangeSlider(frange_excluded=([tes.freq[0], tes.freq[0]], tes.freq[0], tes.freq[-1], tes.freq[2]-tes.freq[0])))
    w_conf['layout'] = bak


    plt.ioff()
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=figsize)

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

    # bld parameter
    # slider_formats = {'lam': "{:.2e}",
    #                   'njoin': "{:d}",
    #                   }

    controls = iplt.plot(tes.freq, tes.sub, nproc=nproc, length=length,
                         **sliders,
                         color='r', linewidth=2, xlim=xlim, ylim=ylim, ax=ax1, play_buttons=True, display_controls=False)

    widget_i = _IntSlider(i=(0, 0, length-1, 1))['i'] if length == 1 else range(length)

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
    ax1.set_title('origial spectra', fontsize=8)
    ax2.set_title('spectra after baseline removed', fontsize=8)

    def return_v(*args, **kwargs): return tes.frange_excluded[0]
    [iplt.axvline(return_v, controls=controls, ax=ax, color='r') for ax in [ax1, ax2, ax3]]
    def return_v(*args, **kwargs): return tes.frange_excluded[1]
    [iplt.axvline(return_v, controls=controls, ax=ax, color='r') for ax in [ax1, ax2, ax3]]

    fig.tight_layout()

    w = controls.controls
    hbs = [
        widgets.HBox([widgets.Label(value=f"1. Select spectra from $start$(max:{controls.controls['start'].max}) to $start+{length}$ to fit:"),
                      w['start'], w['polar']]),
        widgets.HBox([widgets.Label(value=f"2. Fitting:"), w['njoin'], ]),
        widgets.HBox([w['s_method_t'], w['s_sigma_t'],
                     w['s_method_freq'], w['s_sigma_freq']]),
        widgets.HBox([w['exclude_m'], w['deg'], w['sin_f'], w['bound_f']]),
        widgets.HBox(
            [widgets.Label(value=f"3. Show one spectra in top and middle panels:"), w['i']]),
        widgets.HBox([widgets.Label(
            value=f"4. Show the stack spectra of in the range in the bottom panel:"), w['start_stop']]),
        widgets.HBox([w['frange_excluded'], ]),
    ]

    BOX = widgets.VBox(hbs)
    display(BOX)
    plt.show()
    display(BOX)
    return tes, w, BOX
