#!/usr/bin/env python
# coding: utf-8

import numpy as np
from copy import deepcopy

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

def plot_waterfall(f,corvel=None,data = None,xtype = 'freq',polar = 'xx',cmap = 'rainbow',per_vmin_max = None,
                   vmin_max = None,xylim = None,outdir = './',xrange = None,ynsmall = None,
                   plot = True,pdf = None,time_label = False, figsize=(18,9),title = None,**kwargs):
    from matplotlib import pyplot as plt
    
    if pdf is not None:
        plt.switch_backend('agg')
    
    ind_sort = np.argsort(f['mjd'])
    restfreq = 1420.405751#7667
    if xtype == 'freq':
        x = f['freq'][()]
    elif xtype == 'vrad':
        if 'vel' in f.keys():
            x = f['vel'][()]
        else:
            x = 299792.458*(restfreq-freq)/restfreq
    elif xtype == 'vopt':
        freq = f['freq'][()]
        x = 299792.458*(restfreq-freq)/freq
    
    if not isinstance(data,np.ndarray):
        if 'flux' in f.keys():
            T = f['flux'][()][ind_sort]
        elif 'Ta' in f.keys():
            T = f['Ta'][()][ind_sort]
        if 'T' in f.keys():
            T = f['T'][()][ind_sort]

        if len(T.shape) == 2:
            data = T
        else:
            if polar == 'xx':                  
                data = T[:,:,0]
            elif polar == 'yy':                  
                data = T[:,:,1]
    

    
    if xrange != None:
        x1,x2 = np.min(xrange),np.max(xrange)
        is_use = (x>=x1)&(x<=x2)
        x = x[is_use]
        if data.shape[1] != x.shape[0]:
            data = data[:,is_use]
    else:
        xrange = [x[0],x[-1]]
        
    if ynsmall != None:
        if 'ra' in f.keys():
            nspec_use = which_specs(f,nsmall = ynsmall,**kwargs)
        else:
            nspec_use = which_specs(corvel,nsmall = ynsmall,**kwargs)
        mjd = f['mjd'][nspec_use]
        ind_sort2 = np.argsort(mjd)
        data = data[nspec_use,:][ind_sort2]
    
    if plot:
        if time_label:
            mjds = f['mjd'][()]
            extent = (xrange[0],xrange[1],mjds[0],mjds[-1]) 
        else:
            extent = (xrange[0],xrange[1],0,data.shape[0])            

        fig,ax = plt.subplots(figsize=figsize)

        if vmin_max != None:
            im=ax.imshow(data,vmin=vmin_max[0],vmax=vmin_max[1],origin='lower', cmap = cmap,
                         aspect='auto',extent = extent)
        else:
            if per_vmin_max == None:   
                im=ax.imshow(data,origin='lower', cmap = cmap,aspect='auto',extent = extent)
            else:
                #from Fits_Inspect_Overlap.utils import percent_vminmax

                vmin,vmax = percent_vminmax(data,percent = per_vmin_max)
                im=ax.imshow(data,vmin=vmin,vmax=vmax,origin='lower', cmap = cmap,
                            aspect='auto',extent = extent)

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
    
    #return data
    
def get_data(f,polar='average',xrange = None):
    ind_sort = np.argsort(f['mjd'])
    x = f['freq'][()]
    if 'flux' in f.keys():
        T = f['flux'][()][ind_sort]
    elif 'Ta' in f.keys():
        T = f['Ta'][()][ind_sort]
    if 'T' in f.keys():
        T = f['T'][()][ind_sort]
        
    if polar=='average':
        data = np.mean(T,axis = 2)
    elif polar=='xx':
        data =T[:,:,0]
    elif polar=='yy':
        data =T[:,:,1]
    elif polar=='merged':
        data =T
    if xrange != None:
        x1,x2 = np.min(xrange),np.max(xrange)
        is_use = (x>=x1)&(x<=x2)
        x = x[is_use]
        if data.shape[1] != x.shape[0]:
            data = data[:,is_use]
    else:
        xrange = [x[0],x[-1]]   
    
    return np.float64(data)

