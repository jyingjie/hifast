"""
hifast.radec
============

Calculate celestial coordinates (RA, DEC) for FAST drift scan data.

This module matches spectral data with feed cabin position data (KY files)
to compute the Right Ascension (RA) and Declination (DEC) for each spectrum.
It supports interpolation, coordinate transformation using `erfa` or `astropy`,
and environmental corrections.

It can also process a KY file directly to calculate and plot the feed cabin's
trajectory independent of spectral data.
"""

import os
import re
import sys
from glob import glob
import argparse
from .utils.io import *

def create_parser():
    sep_line = '##'+'#'*70+'##'
    parser = ArgumentParser(prog=f"python -m hifast.{os.path.basename(sys.argv[0])[:-3]}",
                        formatter_class=argparse.ArgumentDefaultsHelpFormatter, allow_abbrev=False,
                        description='Calculate RA/DEC coordinates from feed cabin data.')
    
    # --- Input/Output ---
    group = parser.add_argument_group('Input/Output')
    group.add_argument('fname', metavar='FILE',
                        help='Input file path. Can be an HDF5 spectral data file (must contain "S/mjd") or a KY feed cabin log file (.xlsx).')
    group.add_argument('--outdir', metavar='DIR',
                        help='Output directory. Defaults to the input file directory (for HDF5) or current directory (for XLSX).')
    group.add_argument('-f', dest='force', action='store_true',
                        help='Force overwrite of the output file if it already exists.')
    group.add_argument('--plot', action='store_true',
                        help='Create a PDF plot showing the calculated RA/DEC trajectory.')

    # --- Configuration ---
    group = parser.add_argument_group('Configuration')
    group.add_argument('--ky_files', nargs='*', metavar='FILE',
                        help='Explicitly specify one or more KY (feed cabin log) file paths. Separate multiple files with spaces. If omitted, searches in `--ky_dir` based on observation time.')
    group.add_argument('--ky_dir', metavar='DIR',
                        help="Base directory to search for KY files. Defaults to standard locations, with `~/KY` checked first, followed by others like `/data31/KY`.")
    group.add_argument('--backend', choices=['erfa', 'astropy'], default='astropy',
                       help='Coordinate calculation backend. "astropy" (default) is recommended (auto-handles Earth orientation). "erfa" is a lower-level alternative.')
    group.add_argument('--tol', type=float, default=1, metavar='SEC',
                       help='Maximum allowed extrapolation time [seconds]. Coordinates are extrapolated for gaps smaller than this; larger gaps raise an error.')
    group.add_argument('--ky_fixed', action='store_true',
                       help='Assume fixed feed position. Useful for early drift scans with incomplete trajectory logs.')

    # --- Environment ---
    group = parser.add_argument_group('Environment')
    group.add_argument('--phpa', type=float, default=925., metavar='HPA',
                   help='Atmospheric pressure [hPa]. Used for refraction correction.')
    group.add_argument('--temperature', type=float, default=15., metavar='DEG_C',
                   help='Ground-level temperature [Celsius]. Used for refraction correction.')
    group.add_argument('--humidity', type=float, default=0.8, metavar='0-1',
                   help='Relative humidity [0.0 - 1.0]. Used for refraction correction.')
    group.add_argument('--dUT1', metavar='SEC',
                   help='UT1-UTC time difference [seconds]. Only used for "erfa" backend (ignored for "astropy" as it handles it automatically).')

    # --- Performance ---
    group = parser.add_argument_group('Performance')
    group.add_argument('-n', '--nproc', type=int, default=1, metavar='INT',
                       help='Number of parallel processes for multi-beam calculation.')
    group.add_argument('--no_cache', action='store_true',
                       help='Disable using cached intermediate results (forces fresh KY file processing).')
    
    return parser

parser = create_parser()

if __name__ == '__main__':
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
    if args.ky_dir is not None:
        from .core import radec as radec_py
        radec_py.ky_dir_default = [args.ky_dir]
    
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
    if radec is None:
        raise(ValueError('The attempt to process RA-DEC was unsuccessful. Exiting the procedure...'))
        sys.exit(1)
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
