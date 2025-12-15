

__all__ = ['plot', 'colors', 'linestyles', 'linewidths']


import h5py
import os
import re
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt


# for 19 beams
colors = ['k'] + ['C0', 'C1', 'C2', 'C3', 'C4', 'C5']*3
linestyles = ['-'] +  ['-']*6 + ['-.']*6 + ['--']*6
linewidths = [3] + [1]*6 + [2]*6 + [2]*6

def plot(info, frange=[1050, 1450], suptitle='', K_Jy_ratio=False, field_start=None):

    fig, axs = plt.subplots(2, 1, figsize=(10, 14), sharex=True, sharey=True)
    f = info
    freq = f['freq'][:]
    is_use = (freq > frange[0]) & (freq < frange[1])
    freq = freq[is_use]
    keys = []
    if field_start is None:
        for key in f.keys():
            if key.startswith('M'):
                keys.append(key)
        if len(keys) == 0:
            for key in f.keys():
                if key.startswith('K_Jy'):
                    keys.append(key)
    else:
        for key in f.keys():
            if key.startswith(field_start):
                keys.append(key)
    if len(key) > 19:
        raise(ValueError('find more than 19beams?'))
    # plot
    ax = axs[0]
    polar = 0

    if K_Jy_ratio:
        refer_d = {}

    for i, key in enumerate(keys):
        y = f[key][0, :, polar][is_use]
        label = key
        if K_Jy_ratio:
            nB = int(re.sub(r'\D', '', key))
            print(nB)
            ra = f[f'ra{nB}'][:1]
            dec = f[f'dec{nB}'][:1]
            mjd = f[f'mjd'][:1]

            refer_y = gain.Get_gain(ra, dec, mjd, nB, freq)[0][0]
            refer_d[nB] = refer_y
            y = y/refer_y
            label += '_r'
        ax.plot(freq, y, c=colors[i], ls=linestyles[i], lw=linewidths[i], label=label)
    ax.grid()
    ax.set_xlabel('freq')
    ax.set_title(f'Polar:{polar}')
    ## plot
    if f[key].shape[2] == 2:
        polar = 1
        ax = axs[1]

        for i, key in enumerate(keys):
            y = f[key][0, :, polar][is_use]
            label = key
            if K_Jy_ratio:
                nB = int(re.sub(r'\D', '', key))
                refer_y = refer_d[nB]
                y = y/refer_y
                label += '_r'
            ax.plot(freq, y, c=colors[i], ls=linestyles[i], lw=linewidths[i], label=label)
    ax.grid()
    ax.set_title(f'polar:{polar}')
    axs[1].legend(loc=(1.01, 0))
    axs[0].legend(loc=(1.01, 0))

    fig.suptitle(suptitle)
    fig.tight_layout()
    # for ax in axs:
    #     ax.set_ylim(13, 17.9)

    return fig


if __name__ == '__main__':
    import argparse
    import os
    import re

    parser = argparse.ArgumentParser(allow_abbrev=False)
    parser.add_argument('fpaths', nargs='+',
                        help='paths of FluxGain files')
    parser.add_argument('--outdir', default='./',
                       help='default is ./')
    args = parser.parse_args()
    plot_files(args.fpaths, args.outdir)
