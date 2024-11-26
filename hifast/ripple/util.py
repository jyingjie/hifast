#!/usr/bin/env python
# coding: utf-8

import numpy as np
from copy import deepcopy
from hifast.utils.io import BaseIO

def vopt2vrad(vopt):
    """
    velocity: optical to radio
    """
    c = 299792.458
    vrad = c -c**2/(c + vopt)
    return vrad

def vrad2vopt(vrad):
    """
    velocity: radio to optical
    """
    c = 299792.458
    vopt = c**2/(c - vrad) - c
    return vopt

def redshift(v , relative = False):
    c = 299792.458
    beta = v / c
    
    if relative:
        g = 1 / np.sqrt(1 - beta ** 2)
        z = (1 + v / c) * g - 1
        return z
    else:
        return betas

def percent_vminmax(data,percent = None):
    if percent != None:
        vmin = np.nanpercentile(data,q = (1-percent)/2* 100,interpolation='nearest')
        vmax = np.nanpercentile(data,q = (1+percent)/2* 100,interpolation='nearest')
        return vmin,vmax
    else:
        return np.nanmin(data),np.nanmax(data)

def _round_up_to_odd_integer(value):
    c = np.ceil(value)
    if np.isscalar(c):
        if c % 2 == 0:
            r = c + 1
        else:
            r = c
    else:
        r = np.zeros(c.shape[0])
        for i in range(c.shape[0]):
            if c[i] % 2 == 0:
                r[i] = c[i] + 1
            else:
                r[i] = c[i]
    return r.astype('int')

def date2mjd(date):
    """
    date: str; UTC+8
    """
    from astropy.time import Time; import astropy.units as u
    t = Time(date, format='iso', scale='utc') - 8*u.hour
    
    return t.mjd

def do_smooth(s1p, s_method_t = 'none', s_sigma_t = None, 
              s_method_freq = 'none', s_sigma_freq = None, is_rfi = None):
    T = deepcopy(s1p)
    if is_rfi is None: is_rfi = np.full(T.shape[:2], False, dtype=bool)
    is_excluded = np.all(is_rfi,axis = 1)
    T[is_rfi] = 0

    from ..utils.misc import smooth1d
    if s_method_t in ['gaussian', 'boxcar', 'median']:
        print('Smoothing t ...')
        T[~is_excluded] = smooth1d(T[~is_excluded], axis=0, sigma=s_sigma_t, method=s_method_t)
    elif s_method_t == 'iter_median':
        print('Smoothing t ...')
        from .sw_fft import med_fit_ripple
        T[~is_excluded] = med_fit_ripple(T[~is_excluded], nspec = s_sigma_t, func='iter')

    if s_method_freq in ['gaussian', 'boxcar', 'median']:
        print('Smoothing freq ...')
        T[~is_excluded] = smooth1d(T[~is_excluded], axis=1, sigma=s_sigma_freq, method=s_method_freq)
    return T

def do_smooth_onoff(s1p, is_on = None,is_rfi = None, **kwargs):
    T = deepcopy(s1p)
    if is_on is None: 
        T = do_smooth(T,is_rfi = is_rfi, **kwargs)
    else:
        if is_rfi is None: is_rfi = np.full(T.shape[:2], False, dtype=bool)
        T[is_on] = do_smooth(T[is_on],is_rfi = is_rfi[is_on], **kwargs)
        T[~is_on] = do_smooth(T[~is_on],is_rfi = is_rfi[~is_on], **kwargs)
    return T

def line_set(ax,xlabel,ylabel,direction='in',xlim=None,ylim=None,legend=True,
             title=None,loc='best', size = None,frameon=True):
    """
    plot settings
    """
    if size is None:
        size = {'mz': 1,   # Set thickness of the tick marks
                'lz': 3,   # Set length of the tick marks
                'lbz': 14,  # Set label size
                'tkz': 12,  # Set tick size
                }
    mz = size['mz']; lz = size['lz']
    lbz = size['lbz']; tkz = size['tkz']
    # Make tick lines thicker
    for l in ax.get_xticklines():
        l.set_markersize(lz)
        l.set_markeredgewidth(mz)
    for l in ax.get_yticklines():
        l.set_markersize(lz)
        l.set_markeredgewidth(mz)

    # Make figure box thicker
    for s in ax.spines.values():
        s.set_linewidth(mz)
    ax.minorticks_on()
    ax.tick_params("both",which = 'both',direction=direction,labelsize=tkz,
                  bottom=True, top=True,left=True,right=True)
    if xlim is not None: ax.set_xlim(xlim)
    if ylim is not None: ax.set_ylim(ylim)
    ax.set_xlabel(xlabel,fontsize=lbz)
    ax.set_ylabel(ylabel,fontsize=lbz)
    if legend: ax.legend(loc = loc,fontsize=tkz, frameon=frameon)
    if title is not None: ax.set_title(title,fontsize=lbz)

    
class Args(object):
    def __init__(self, fpath, frange = None, outdir = None,):
        self.fpath = fpath
        self.outdir = outdir
        self.frange = frange

class Read_hdf5(BaseIO):
    ver = 'old'
    
    def __init__(self, args, dict_in=None, inplace_args=False, HistoryAdd=None):
        """
        args: class
              including attributes: fpath, outdir, frange
        dict_in: if None, load data from args.fpath, if set, omit data in args.fpath
        """
        self.args = args if inplace_args else deepcopy(args)
        self.dict_in = dict_in
        self.HistoryAdd = HistoryAdd
        self._gen_fpath_out()
        if self.dict_in is None:
            self._check_fout()
        self.nB = self.get_nB(self.args.fpath)
        self._import_m()
        self.open_fpath()
        self.load_specs()
        self.load_radec()
        try:
            self.load_and_add_Header()
        except TypeError:
            pass
    
    def get_data(self, polar = 'none'):
        data = deepcopy(self.s2p)
        if data.shape[0] == 2 or data.shape[0] == 1:
            from ..utils.io import PolarMjdChan_to_MjdChanPolar
            data = PolarMjdChan_to_MjdChanPolar(data)
        if len(data.shape) == 3:
            if polar == 'xx':
                data = data[:,:,0]
            elif polar == 'yy':
                data = data[:,:,1]
            elif polar == 'average':
                data = np.mean(data,axis = 2)
        return data
    
    def plot_waterfall(self,data = None,xtype = 'freq', cmap = 'rainbow',per_vmin_max = None, polar = None, 
                       vmin_max = None,xylim = None,outdir = './',xrange = None,interp_method = 'nearest',
                       plot = True,pdf = None,time_label = False, figsize=(15,4),title = None,**kwargs):
        from matplotlib import pyplot as plt
        #import matplotlib
        #matplotlib.rcParams['image.interpolation'] = 'none'

        if pdf is not None:
            plt.switch_backend('agg')

        restfreq = 1420.405751#7667
        c = 299792.458
        freq = self.freq
        if xtype == 'freq':
            x = self.freq
        elif xtype == 'vrad':
            if 'vel' in self.fs.keys():
                x = self.fs['vel'][()]
            else:
                x = c*(restfreq-freq)/restfreq
        elif xtype == 'vopt':
            x = c*(restfreq-freq)/freq
        
        if not isinstance(data,np.ndarray):
            data = self.get_data(polar)
        
        if len(data.shape) != 2:
            raise ValueError(f"Check your input data shape {data.shape}. Are they 2D? ")

        if xrange != None:
            x1,x2 = np.min(xrange),np.max(xrange)
            is_use = (x>=x1)&(x<=x2)
            x = x[is_use]
#             if data.shape[1] != x.shape[0]:
#                 data = data[:,is_use]
        else:
            xrange = [x[0],x[-1]]

        if plot:
            if time_label:
                mjds = self.mjd
                extent = (xrange[0],xrange[1],mjds[0],mjds[-1]) 
            else:
                extent = (xrange[0],xrange[1],0,data.shape[0])            

            fig,ax = plt.subplots(figsize=figsize)

            if vmin_max != None:
                im=ax.imshow(data,vmin=vmin_max[0],vmax=vmin_max[1],origin='lower', cmap = cmap,
                             aspect='auto',extent = extent,interpolation=interp_method,)
            else:
                if per_vmin_max == None:   
                    im=ax.imshow(data,origin='lower', cmap = cmap,aspect='auto',extent = extent,
                                 interpolation=interp_method,)
                else:
                    vmin,vmax = percent_vminmax(data,percent = per_vmin_max)
                    im=ax.imshow(data,vmin=vmin,vmax=vmax,origin='lower', cmap = cmap,
                                aspect='auto',extent = extent,interpolation=interp_method,)

            ax.set_xlabel(xtype)
            ax.set_ylabel("specs")
            #ax.tick_params(labelsize=14)    
            ax.set_title(title)
            plt.colorbar(im,pad=.01)
            plt.tight_layout()
            if xylim is not None:
                xlim1,xlim2,ylim1,ylim2 = xylim
                ax.set_xlim(xlim1,xlim2)
                ax.set_ylim(ylim1,ylim2)
            ax.minorticks_on()
            if time_label:
                from astropy.time import Time; import astropy.units as u
                labels_new = ((Time(ax.axes.get_yticks(), format='mjd')) + 8*u.hour).strftime('%H:%M')
                ax.axes.set_yticklabels(labels_new)

            if pdf is not None:
                pdf.savefig();plt.close()
            else:
                plt.show()

