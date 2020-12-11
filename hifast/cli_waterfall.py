#!/usr/bin/env python
# coding: utf-8

import h5py
from glob import glob
import matplotlib.pyplot as plt
import numpy as np
import re
import os
from scipy import interpolate as inter

#import mpl_rcParams_paper
plt.switch_backend('agg')

def _tight_ra(ra):
    ra_s= np.sort(ra)
    diff_s= np.diff(ra_s)
    ind_max= np.argmax(diff_s)
    if diff_s[ind_max] < ra_s[0]+360-ra_s[-1]:
        return ra
    else:
        ra= np.copy(ra)
        is_c= ra>=ra_s[ind_max+1]
        ra[is_c]= ra[is_c]-360
        return ra

def plot(fname, ax=None, imshow_kwargs={}, vlines=[-80, -600], colorbar=True, sec_ytick=False, xrange=None):
    f = h5py.File(fname,'r')
    if 'T' in f.keys():
        vals = f['T']
    elif 'Ta' in f.keys():
        vals = f['Ta']
    elif 'flux' in f.keys():
        vals = f['flux']
    else:
        raise()
    if 'vel' in f.keys():
        f_axis = f['vel'][()]
    elif 'freq' in f.keys():
        f_axis = f['freq'][()]
    else:
        raise()
    try:
        t_axis = f['ra'][()]
        t_axis = _tight_ra(t_axis)
        t_axis_2 = f['dec'][()]
    except:
        t_axis = t_axis_2 = np.arange(len(vals))

    if vals.ndim ==3:
        vals = vals[...,0]
    
    if xrange is not None:
        is_ = (f_axis >= xrange[0]) & (f_axis <= xrange[1])
        f_axis = f_axis[is_]
        vals = vals[:,is_]
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(15,15))
#     if imshow_kwargs == {}:
#         imshow_kwargs = {}
    else:
        fig = None
    im = ax.imshow(vals[:], aspect='auto', extent=(f_axis[0], f_axis[-1], t_axis[0], t_axis[-1]), origin='lower', **imshow_kwargs)
    #rasterized=False
    if vlines is not None:
        for vline in vlines:
            ax.axvline(x=vline, color='r', linestyle=':')
    if sec_ytick:
        try:
            def forward(x):
                return inter.interp1d(t_axis, np.linspace(np.min(t_axis_2),np.max(t_axis_2),len(t_axis_2[:])), fill_value='extrapolate')(x)
            def inverse(x):
                return inter.interp1d(np.linspace(np.min(t_axis_2),np.max(t_axis_2),len(t_axis_2[:])), t_axis, fill_value='extrapolate')(x)
            ax.secondary_yaxis('right', functions=(forward, inverse))
        except:
            pass
    ax.minorticks_on()
    f.close()
    return im, fig

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(allow_abbrev=False)
    parser.add_argument('fnames', nargs='+',
                        help='fnames')
    parser.add_argument('--outdir', default='./',
                       help='default is ./')
    parser.add_argument('-s', '--single', action='store_true',
                       help="19 beams in 19 figures")
    parser.add_argument('--xrange', type=float, nargs=2,
                       help='x axis range')
    parser.add_argument('--vmin', type=float,
                       help='')
    parser.add_argument('--vmax', type=float,
                       help='')
    parser.add_argument('--vlines', type=float, nargs='+',
                       help='')
    
    args = parser.parse_args()
    fnames = args.fnames
    vlines = args.vlines
    outdir = args.outdir
    xrange = args.xrange
    
    fnames.sort()
    single = False
    if len(fnames)==1:
        single = True
    single = args.single
    imshow_kwargs = {}
    imshow_kwargs['vmax'] = args.vmax
    imshow_kwargs['vmin'] = args.vmin
    imshow_kwargs['cmap'] = 'jet'
    
    from tqdm import tqdm
    if not single:
        import pandas as pd
        #keys = list(map(lambda x: re.sub('-M[0-1][0-9]','-M00', os.path.basename(x).split('-specs_T')[0]), fnames))
        keys = list(map(lambda x: re.sub('-M[0-1][0-9]','-M00', os.path.basename(x)), fnames))
        files = pd.DataFrame({'key':keys, 'fname':fnames})
        for key, _files in tqdm(files.groupby('key')):
            nrows=4
            ncols=5
            fig, axs = plt.subplots(nrows, ncols, figsize=(160/3,90/3), sharex=True, sharey=True)
            axs = axs.flatten()
            for i, (fname, ax) in enumerate(zip(_files['fname'], axs[:19])):
                sec_ytick = True if (i+1)%ncols==0 else False
                im, _ = plot(fname, ax, imshow_kwargs, vlines=vlines, sec_ytick=sec_ytick, xrange=xrange)
            fig.colorbar(im, ax= axs[-1])
            ax = axs[-1]
            ax.plot([], [], label=key)
            ax.legend()
            
            fig.tight_layout()
            fig.savefig(f'{outdir}/{key}.19.pdf')
            fig.clear()
    else:
        for fname in tqdm(fnames):
            fbasename = os.path.basename(fname)
            nrows=1
            ncols=1
            fig, ax = plt.subplots(nrows, ncols, figsize=(15,12), sharex=True, sharey=True)
            sec_ytick = True
            im, _ = plot(fname, ax, imshow_kwargs, vlines=vlines, sec_ytick=sec_ytick, xrange=xrange)
            fig.colorbar(im, ax= ax)
            ax.set_title(fbasename)
            
            fig.tight_layout()
            fig.savefig(f'{outdir}/{fbasename}.pdf')
            fig.clear()
