"""
hifast.cube
===========

This module handles the regridding of spectra into data cubes.
It supports various convolution kernels (Gaussian, Bessel-Gaussian, Sinc-Gaussian)
and allows flexible configuration of the output WCS grid.
"""

__all__ = ['sep_line', 'parser']


from .utils.io import *


sep_line = '##'+'#'*70+'##'
parser = ArgumentParser(prog=f"python -m hifast.{os.path.basename(sys.argv[0])[:-3]}",
                        formatter_class=formatter_class, allow_abbrev=False,
                        description='Regrid spectra into a data cube.',
                        epilog='For more details, see the full documentation at https://hifast.readthedocs.io')


# --- Input/Output Arguments ---
group = parser.add_argument_group('Input/Output')
group.add_argument('fpatterns', nargs='+',
                    help='Input file names or pattern (e.g., "data/*.hdf5").')
group.add_argument('--outname', required=True,
                    help='Output FITS file name (full path).')
group.add_argument('-f', action='store_true', dest='force',
                    help='Overwrite the output file if it already exists.')
group.add_argument('-k', '--key', choices=['flux', 'Ta', 'Power'],
                       help='Data field to use from the input files: "flux", "Ta", or "Power". Default: Auto-detected.')
group.add_argument('--polar', choices=['XX','YY','M'], default='M',
                       help='Polarization to process: "XX", "YY", or "M" (Merge/Average).')
group.add_argument('--wcs_from', metavar='FILE',
                   help='Path to a FITS file to copy spatial WCS parameters (RA/DEC) from. Note: The spectral axis is NOT copied.')

# --- Grid & WCS Arguments ---
group = parser.add_argument_group('Grid & WCS')
group.add_argument('--ra_range', type=float, nargs=2, metavar=('MIN', 'MAX'),
                   help='RA range for the output grid [deg]. Followed by two values (min, max). Default: Derived from input files.')
group.add_argument('--dec_range', type=float, nargs=2, metavar=('MIN', 'MAX'),
                   help='DEC range for the output grid [deg]. Followed by two values (min, max). Default: Derived from input files.')
group.add_argument('--wcs_ra_center', type=float, metavar='DEG',
                   help='Manually set the CRVAL1 (RA center) in the WCS header [deg].')
group.add_argument('--wcs_dec_center', type=float, metavar='DEG',
                   help='Manually set the CRVAL2 (DEC center) in the WCS header [deg].')
group.add_argument('--range3', type=float, nargs=2, metavar=('MIN', 'MAX'),
                   help='Range for the third axis (velocity or frequency). Followed by two values (min, max).')
group.add_argument('--type3', choices=['vopt', 'vrad', 'freq'], default='vrad',
                   help='Type of the third axis: "vopt" (Optical Velocity), "vrad" (Radio Velocity), or "freq" (Frequency). '
                        'Note: Velocity is always calculated from the frequency information in the input data; '
                        'existing velocity columns in the input are ignored.')
group.add_argument('--bwidth', type=float, default=[60.], nargs='+', metavar='ARCSEC',
                   help='Grid pixel size [arcsec]. Followed by one value (used for both RA and DEC) or two values (RA_pix, DEC_pix).')
group.add_argument('-p', '--proj', default='AIT',
                   help='WCS projection type. Common options: "SIN", "TAN", "CAR", etc. See arXiv:astro-ph/0207413. '
                        'Special mode: "rCAR" forces a Plate Carree projection '
                        'resulting in straight and orthogonal RA/DEC lines.')

# --- Convolution Arguments ---
group = parser.add_argument_group('Convolution')
group.add_argument('-m', '--method', default='gaussian', choices=['gaussian'],
                   help='Convolution kernel type. Currently only "gaussian" is fully supported.')
group.add_argument('--beam_fwhw', type=float, default=2.9, metavar='ARCMIN',
                   help='Beam Full Width at Half Maximum (FWHM) [arcmin].')
group.add_argument('--gaussian_fwhw', type=float, metavar='ARCMIN',
                   help='FWHM of the Gaussian convolution kernel [arcmin]. Default: beam_fwhw / 2.')
group.add_argument('--apply_beam_correction', type=bool_fun, choices=[True, False], required=True, metavar='BOOL',
                   help=('Apply beam correction. This updates BMAJ/BMIN in the header to reflect the new effective resolution '
                         'and applies a flux correction factor for Jy/beam units.'))
group.add_argument('--r_cut', type=float, metavar='ARCSEC',
                   help="Cutoff radius for the convolution kernel [arcsec]. " +
                        "Spectra within this radius from a grid point are included. " +
                        "Default (Gaussian): 3 * sigma, where sigma is derived from --gaussian_fwhw (sigma ≈ gaussian_fwhw / 2.355).")
group.add_argument('--frac_finite_min', type=float, default=1, metavar='FLOAT',
                   help='Minimum fraction of finite values required in a channel to produce a valid output. ' +
                        'If the fraction is lower, the output is NaN. (requires all spectra to be valid).')

# --- Performance Arguments ---
group = parser.add_argument_group('Performance')
group.add_argument('--nproc', '-n', type=int, default=1, metavar='INT',
                   help='Number of processes to use for parallel processing.')
group.add_argument('--step', type=int, default=5, metavar='INT',
                   help='Number of files to process in each batch to control memory usage.')
group.add_argument('--share_mem', type=bool_fun, choices=[True, False], default='False', metavar='BOOL',
                   help='Enable shared memory for multiprocessing (Linux only). When enabled, uses /dev/shm by default.')
group.add_argument('--temp_dir', type=str, default='/dev/shm', metavar='DIR',
                   help='Temporary directory for multiprocessing when --share_mem is True. ' +
                        '(fast, in-memory). ' +
                        'If /dev/shm is too small, specify a disk directory with more space. ' +
                        'Note: Remember to clean up residual temp files after processing.')
group.add_argument('--use_optimized_search', type=bool_fun, choices=[True, False], default='False', metavar='BOOL',
                   help='Use optimized search_around_sky implementation for better performance. ' +
                        'Set to False to use original astropy implementation for comparison.')

# --- Testing/Hidden Arguments ---
parser.add_argument('--scale_beams_file', help=argparse.SUPPRESS)


if __name__ == '__main__':
    import os
    from .core.image import Imaging
    args_ = parser.parse_args()

    img = Imaging(args_)
    img()
