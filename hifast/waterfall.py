

__all__ = ['plot_im', 'plot']


import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as interp
import h5py


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


def plot_im(vals, x=None, y=None, x_plot_range=None, y_plot_range=None, ax=None, vlines=None, colorbar=True, y2=None,
         vmin='per0.01', vmax='per95', **kwargs):
    """
    vals: shape (m,n)
    x: None or shape (n,)
    y: None or shape (m,)
    """

    if x is None:
        x = np.arange(vals.shape[1])
    if y is None:
        y = np.arange(vals.shape[0])

    if x_plot_range is not None:
        is_ = (x >= x_plot_range[0]) & (x <= x_plot_range[1])
        x = x[is_]
        vals = vals[:,is_]
    if y_plot_range is not None:
        is_ = (y >= y_plot_range[0]) & (y <= y_plot_range[1])
        y = y[is_]
        vals = vals[is_,:]
        if y2 is not None:
            y2 = y2[is_]
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(10,10))
#     if imshow_kwargs == {}:
#         imshow_kwargs = {}
    else:
        fig = None
    # convert to value
    if vmax is not None:
        vmax = float(vmax) if 'per' != str(vmax)[:3] else np.nanpercentile(vals, float(str(vmax)[3:]), interpolation='nearest')
    if vmin is not None:
        vmin = float(vmin) if 'per' != str(vmin)[:3] else np.nanpercentile(vals, float(str(vmin)[3:]), interpolation='nearest')

    im = ax.imshow(vals[:], aspect='auto', extent=(x[0], x[-1], y[0], y[-1]), origin='lower', vmin=vmin, vmax=vmax, **kwargs)
    #rasterized=False
    if vlines is not None:
        for vline in vlines:
            ax.axvline(x=vline, color='r', linestyle=':')
    if y2 is not None:
        try:
            def forward(x):
                return interp.interp1d(y, y2, fill_value='extrapolate')(x)
            def inverse(x):
                return interp.interp1d(y2, y, fill_value='extrapolate')(x)
            ax.secondary_yaxis('right', functions=(forward, inverse))
        except:
            pass
    ax.minorticks_on()
    if colorbar:
        plt.colorbar(im, ax=ax)
    return im, fig


def plot(fname, ax, polar=0, ytick1='index', ytick2=None, replace_rfi = False, **kwargs):
    """
    fname:
    ax:
    polar: polarization, 0 or 1
    ytick1: str; index, ra, dec, mjd or time
    ytick2: str; index, ra, dec, mjd or None
    """
    f = h5py.File(fname,'r')
    S = f['S']
    if 'T' in S.keys():
        vals = S['T']
    elif 'Ta' in S.keys():
        vals = S['Ta']
    elif 'flux' in S.keys():
        vals = S['flux']
    else:
        raise(ValueError('can not find spec'))

    if 'vel' in S.keys():
        x = S['vel'][()]
    elif 'freq' in S.keys():
        x = S['freq'][()]
    else:
        raise()

    if vals.ndim == 3:
        if polar == -1:
            vals = np.mean(vals, axis = 0)
        else:
            vals = vals[polar]

    if replace_rfi:
        if 'is_rfi' in S.keys():
            is_rfi = S['is_rfi'][:]
            vals[is_rfi] = np.nan

    if ytick1 in ['ra', 'dec', 'mjd', 'time']:
        key = 'mjd' if ytick1=='time' else ytick1
        try:
            y = S[key][()]
        except Exception as err:
            print(err)
    elif ytick1 == 'index':
        y = np.arange(vals.shape[1])
    y2 = None
    if ytick2 in ['ra', 'dec', 'mjd', 'time']:
        key = 'mjd' if ytick2=='time' else ytick2
        try:
            y2 = S[key][()]
        except Exception as err:
            print(err)
    elif ytick2 == 'index':
        y2 = np.arange(vals.shape[1])

    im, _ = plot_im(vals, x, y=y, y2=y2, ax=ax, **kwargs)
    if ytick1=='time':
        from astropy.time import Time, TimezoneInfo
        import astropy.units as u
        labels_new = ((Time(ax.axes.get_yticks(), format='mjd')) + 8*u.hour).strftime('%H:%M:%S')
        ax.axes.set_yticklabels(labels_new)
    return im

def _check_fout(args, fpath_out):
    """
    check whether fout exists
    """
    ret = None
    if os.path.exists(fpath_out):
        if args.force:
            print(f"will overwrite the existing output file {fpath_out}")
        else:
            print(f"File exists {fpath_out}")
            print("exit... Use ' -f ' to overwrite it.")
            ret = "exit"
    return ret


if __name__ == '__main__':
    import argparse
    import os
    import re

    parser = argparse.ArgumentParser(allow_abbrev=False)
    parser.add_argument('fnames', nargs='+',
                        help='fnames')
    parser.add_argument('--outdir', default='./',
                       help='default is ./')
    parser.add_argument('-s', '--single', action='store_true',
                       help="")
    parser.add_argument('-f', action='store_true', dest='force',
                        help='if set, overwriting file if output file exists')
    parser.add_argument('--xrange', type=float, nargs=2,
                       help='x axis range')
    parser.add_argument('--vmin', type=str, default='per0.1',
                       help='')
    parser.add_argument('--vmax', type=str, default='per95',
                       help='')
    parser.add_argument('--interpolation', type=str, default='nearest',
                       help="'none', 'antialiased', 'nearest', 'bilinear', 'bicubic', 'spline16', 'spline36',")
    parser.add_argument('--ytick1', choices=['index', 'ra', 'dec', 'mjd', 'time'], default='index',
                       help='')
    parser.add_argument('--ytick2', choices=['index', 'ra', 'dec', 'mjd'],
                       help='')
    parser.add_argument('--vlines', type=float, nargs='+',
                       help='')
    parser.add_argument('--show', action='store_true',
                       help='')
    parser.add_argument('--replace_rfi', action='store_true',
                       help='rfi')
    parser.add_argument('--polar', type=int, default=0,
                       help='polar')

    args = parser.parse_args()
    fnames = args.fnames
    vlines = args.vlines
    outdir = args.outdir
    xrange = args.xrange

    if not args.show:
        plt.switch_backend('agg')
    fnames.sort()
    single = args.single
    if len(fnames)==1:
        single = True
    rasterized = True if args.interpolation == 'none' else False
    imshow_kwargs = {}
    imshow_kwargs['interpolation'] = args.interpolation

    from tqdm import tqdm
    if not single:
        from itertools import groupby
        
        #keys = list(map(lambda x: re.sub('-M[0-1][0-9]','-M00', os.path.basename(x).split('-specs_T')[0]), fnames))
        keys = list(map(lambda x: re.sub('-M[0-1][0-9]','-M00', os.path.basename(x)), fnames))
        
        # Sort by keys for groupby
        combined = sorted(zip(keys, fnames), key=lambda x: x[0])
        
        for key, group in tqdm(groupby(combined, lambda x: x[0])):
            # group is iterator of (key, fname)
            _files_fname = [x[1] for x in group]
            fpath_out = f'{outdir}/{key}.19.pdf'
            ret = _check_fout(args, fpath_out)
            if ret == "exit": continue

            nrows=4
            ncols=5
            fig, axs = plt.subplots(nrows, ncols, figsize=(160/3,90/3), sharex=True, sharey=True)
            axs = axs.flatten()
            for i, (fname, ax) in enumerate(zip(_files_fname, axs[:19])):
                im = plot(fname, ax, polar=args.polar, ytick1=args.ytick1, ytick2=args.ytick2, vmin=args.vmin, vmax=args.vmax,
                          x_plot_range=xrange, vlines=vlines, colorbar=False,replace_rfi = args.replace_rfi,
                          **imshow_kwargs)
            fig.colorbar(im, ax= axs[-1])
            ax = axs[-1]
            ax.plot([], [], label=key)
            ax.legend()

            fig.tight_layout()
            if args.show:
                fig.show()
                input()
            fig.savefig(fpath_out)
            fig.clear()
    else:
        for fname in tqdm(fnames):
            fbasename = os.path.basename(fname)
            fpath_out = f'{outdir}/{fbasename}.pdf'
            ret = _check_fout(args, fpath_out)
            if ret == "exit": continue

            nrows=1
            ncols=1
            fig, ax = plt.subplots(nrows, ncols, figsize=(15,12), sharex=True, sharey=True)
            im = plot(fname, ax, polar=args.polar, ytick1=args.ytick1, ytick2=args.ytick2, vmin=args.vmin, vmax=args.vmax,
                      x_plot_range=xrange, vlines=vlines, colorbar=False,replace_rfi = args.replace_rfi,
                      **imshow_kwargs)
            fig.colorbar(im, ax= ax)
            ax.set_title(fbasename)

            fig.tight_layout()
            if args.show:
                fig.show()
                input()
            fig.savefig(fpath_out)
            fig.clear()
