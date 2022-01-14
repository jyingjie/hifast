#!/usr/bin/env python
# coding: utf-8

import os
import re
import sys
from glob import glob

if __name__ == '__main__':
    from .utils.io import *

    sep_line = '##'+'#'*70+'##'
    parser = ArgumentParser(prog=f"python -m hifast.{os.path.basename(sys.argv[0])[:-3]}",
                        formatter_class=formatter_class, allow_abbrev=False,
                        description='Calculate RA DEC', )
    parser.add_argument('fname',
                        help='file name; hdf5 file with "S/mjd" filed or KY file(.xlsx)')
    parser.add_argument('-f', dest='force', action='store_true',
                        help='overwriting file if out file exists')
    parser.add_argument('--ky_files', nargs='*',
                        help='KY files, if not given, guessing from fname')
    parser.add_argument('--backend', choices=['erfa', 'astropy'], default='astropy',
                       help='using erfa or astropy. The astropy consider dUT1, xp, yp')
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
    
    group = parser.add_argument_group(f'environment parameters')
    group.add_argument('--phpa', type=float, default=925.,
                   help='atmospheric pressure in hPa')
    group.add_argument('--temperature', type=float, default=15.,
                   help='The ground-level temperature in deg C.')
    group.add_argument('--humidity', type=float, default=0.8,
                   help='The relative humidity as a dimensionless quantity between 0 to 1')
    parser.add_argument('--dUT1', default=0.1,
                   help='UT1-UTC')
    
    
    

    args = parser.parse_args()
#     print('#'*35+'Args'+'#'*35)
#     print(parser.format_values())  # useful for logging where different settings came from
#     print('#'*35+'####'+'#'*35)

    fname = args.fname
    ky_files = args.ky_files
    tol = args.tol
    ky_fixed = args.ky_fixed
    use_cache = not args.no_cache
    plot = args.plot
    outdir = args.outdir
    nproc = args.nproc
    
    outpart = '-radec'
    if outdir is None:
        if fname[-5:] == '.xlsx':
            outdir = './'
        else:
            outdir = os.path.dirname(fname)
    fileout = os.path.join(outdir, '.'.join(os.path.basename(fname).split('.')[:-1]) + f'{outpart}.hdf5')
    if os.path.exists(fileout):
        if args.force:
            print(f"will overwrite the existing out file {fileout}")
        else:
            print(f"File exists {fileout}")
            print('exit... Using -f to overwrite it.')
            sys.exit()
            
    from .core.radec import plot_radec, get_radec
    env_para = {}
    env_para['phpa'] = args.phpa
    env_para['temperature'] = args.temperature
    env_para['humidity'] = args.humidity
    try:
        dUT1 = float(args.dUT1)
    except:
        dUT1 = args.dUT1
    radec = get_radec(fname, ky_files=ky_files, tol=tol, ky_fixed=ky_fixed, use_cache=use_cache, nproc=nproc, 
                      backend=args.backend, env_para=env_para, dUT1=dUT1)
    #saving
    ##record history
    from .utils.io import *
    import json
    header = rec_his(args=json.dumps(args.__dict__))
    print('Saving...')
    radec['Header'] = header
    save_specs_hdf5(fileout, radec)
    print(f"Saved to {fileout}")
    if plot:
        plot_radec(radec, fileout + '.pdf')
