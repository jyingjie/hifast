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

def do_smooth(s1p, s_method_t, s_sigma_t, s_method_freq, s_sigma_freq, is_rfi = None):
    T = deepcopy(s1p)
    if is_rfi is None: is_rfi = np.full(T.shape, False, dtype=bool)
    is_excluded = np.all(is_rfi,axis = 1)
    T[is_rfi] = 0

    from hifast.utils.misc import smooth1d
    if s_method_t in ['gaussian', 'boxcar', 'median']:
        print('Smooth ing ...')
        T[~is_excluded] = smooth1d(T[~is_excluded], axis=0, sigma=s_sigma_t, method=s_method_t)

    if s_method_freq in ['gaussian', 'boxcar', 'median']:
        print('Smooth ing ...')
        T[~is_excluded] = smooth1d(T[~is_excluded], axis=1, sigma=s_sigma_freq, method=s_method_freq)
    return T
    
class Args(object):
    def __init__(self, fpath, frange = None, outdir = None,):
        self.fpath = fpath
        self.outdir = outdir
        self.frange = frange

class Read_hdf5(BaseIO):
    ver = 'old'
    
    def load_and_add_Header(self,):
        pass
    
    def get_data(self, polar = 'none'):
        data = deepcopy(self.s2p)
        if data.shape[0] == 2 or data.shape[0] == 1:
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
        if data.shape[1] != len(freq):
            data = data.T

        if xrange != None:
            x1,x2 = np.min(xrange),np.max(xrange)
            is_use = (x>=x1)&(x<=x2)
            x = x[is_use]
            if data.shape[1] != x.shape[0]:
                data = data[:,is_use]
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

