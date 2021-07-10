#!/usr/bin/env python
# coding: utf-8

import os
import re
import sys
from glob import glob

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('fname',
                        help='file name; hdf5 file with "mjd" filed or KY file(.xlsx)')
    parser.add_argument('-f', '--force', action='store_true',
                        help='overwriting file if out file exists')
    parser.add_argument('--ky_files', nargs='*',
                        help='KY files, if not given, guessing from fname')
    parser.add_argument('--tol', type=float, default=1,
                       help='max allowed extrapolate time; unit: second')
    parser.add_argument('--ky_fixed', action='store_true',
                       help='')
    parser.add_argument('--no_cache', action='store_true',
                       help='whether use the cached')
    parser.add_argument('-n', '--nproc', type=int, default=1,
                       help='parallel process number')
    parser.add_argument('--plot', action='store_true',
                       help='plot the ra dec in pdf image')
    parser.add_argument('--outdir',
                       help='output file directory; Default is same with input file if input hdf5 file, "./" if input .xlsx file.')
    

    args = parser.parse_args()

    fname = args.fname
    ky_files = args.ky_files
    tol = args.tol
    ky_fixed = args.ky_fixed
    use_cache = not args.no_cache
    plot = args.plot
    outdir = args.outdir
    nproc = args.nproc
    
    outpart = '-radec'
    if outdir is None: outdir = os.path.dirname(fname)
    fileout = os.path.join(outdir, '.'.join(os.path.basename(fname).split('.')[:-1]) + f'{outpart}.hdf5')
    if os.path.exists(fileout):
        if args.force:
            print(f"will overwrite the existing out file {fileout}")
        else:
            print(f"File exists {fileout}")
            print('exit... Using -f to overwrite it.')
            sys.exit()
    from .core.radec import plot_radec, get_radec
    radec = get_radec(fname, ky_files=ky_files, tol=tol, ky_fixed=ky_fixed, use_cache=use_cache, nproc=nproc)
    #saving
    ##record history
    from .utils.io import *
    header = rec_his(args=args)
    print('Saving...')
    save_dict_hdf5(fileout, radec, header=header)
    print(f"Saved to {fileout}")
    if plot:
        plot_radec(radec, fileout + '.pdf')
