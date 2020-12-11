#!/usr/bin/env python
# coding: utf-8

import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u

def _get_start_stop(arr,arr_in):
    """
    arr_in: array
       sorted
    """
    if not np.all(arr_in[:-1] <= arr_in[1:]):
        raise('arr_in should be sorted')
    inds_left= np.searchsorted(arr_in, arr, side='left',)
    inds_right=np.searchsorted(arr_in, arr, side='right',)
    return inds_left,inds_right

def pixel_spec(spec, dis, method='reweight', sigma=1.275088, statistic='median'):
    """
    -----------------
    spec: flux
    dis: degree
    sigma: arcmin
    """
    if method=="reweight":
        #Barnes el. al. 2001, MNRAS 322, 486 https://ui.adsabs.harvard.edu/abs/2001MNRAS.322..486B/abstract
        sigma= sigma/60*u.deg
        if statistic=='median':
            wei_m= np.nanmedian(np.exp(- (dis/sigma)**2/2))
            return np.nanmedian(spec,axis=0)/wei_m
        elif statistic=='mean':
            wei_m= np.nanmean(np.exp(- (dis/sigma)**2/2))
            return np.nanmean(spec,axis=0)/wei_m
    elif method=='mean':
        return np.nanmean(spec,axis=0)
    elif method=='median':
        return np.nanmedian(spec,axis=0)
    elif method=='gaussian':
        wei= np.exp(- (dis/sigma)**2/2)
        return np.nansum(spec*wei.reshape((-1,)+(1,)*(spec.ndim-1)),axis=0)/np.nansum(wei)
    else:
        raise(ValueError('method'))

def gridding(ra, dec, spectra, ra_grid, dec_grid, r=1.5, **kwargs):
    """
    ra, dec: array, shape (m,); degree
        The ra dec of observed spectra; degree
    spectra: array, shape (m,n)
        flux
    ra_grid, dec_grid: array, shape (x,y) or (x,); degree
        The center of the grid 
    r: scalar; arcmin
       The spectra separated from the center of the grid less than r will be considered. 
    ------------
    other parameters
    sigma: float; unit is arcmin; default is 1.5
    statistic: str; median or mean, default is median
    """
    if np.isscalar(ra) or np.isscalar(dec) or np.isscalar(ra_grid) or np.isscalar(dec_grid):
        raise ValueError('One of the inputs is a scalar.')
    cata = SkyCoord(ra,dec,unit=(u.degree, u.degree))
    grid= SkyCoord(ra_grid,dec_grid, unit=(u.degree, u.degree))
    grid_ori_ndim=grid.ndim
    if grid_ori_ndim==1:
        grid=grid[:,None]
    grid_f= grid.flatten()
    #find the spec in r arcmin
    ind_g, ind_cata, d2d, d3d=cata.search_around_sky(grid_f,r*u.arcmin)
    ind_g_uni=np.unique(ind_g)# index in grid flatten
    start, stop= _get_start_stop(ind_g_uni,ind_g)
    
    out= np.full(spectra.shape[1:2] + grid.shape, np.nan)
    nums= np.zeros(grid.shape)
    for i,start_, stop_ in zip(ind_g_uni,start, stop):
        dis= d2d[start_ : stop_ ] #distance of spec from the center of grid
        ind_use_= ind_cata[start_ : stop_ ] # index in cata, ra, dec, spectra
        spec_= spectra[ind_use_]
        m,n=i//grid.shape[1], i%grid.shape[1] #index in grid before flatten
        nums[m,n]=len(dis)
        out[:,m,n]= pixel_spec(spec_, dis, **kwargs)
    
    if grid_ori_ndim==1:
        out= out[:,:,0]
        nums= nums[:,0]
    return out, nums