
"""
hifast.bld
==========

Baseline fitting and subtraction for spectral data.

This module provides tools to fit and subtract baselines from spectra. It supports various fitting methods (PLS, polynomial, spline) and preprocessing steps (smoothing, binning) along both time and frequency axes.
"""


from .utils.io import *


__all__ = ['IO', 'IO_i', 'create_parser', 'parser']

sep_line = '##'+'#'*70+'##'

def create_parser():
    parser = ArgumentParser(prog=f"python -m hifast.{os.path.basename(sys.argv[0])[:-3]}",
                            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                            allow_abbrev=False,
                            description='Fit and subtract the baseline from spectral data.')
    add_common_argument(parser)
    
    # Input/Output
    parser.add_argument('fpath', metavar='FILE',
                        help='Path to the input HDF5 file containing spectra (e.g., "data.hdf5").')
    parser.add_argument('--frange', type=float, nargs=2, default=[0, float('inf')], metavar=('START', 'END'),
                        help='Frequency range to process [MHz]. Default: process all.')
    parser.add_argument('--no_radec', action='store_true', not_in_write_out_config_file=True,
                        help="Skip verification/addition of RA/Dec coordinates.")
    parser.add_argument('--show_prog', type=bool_fun, choices=[True, False], default='True', env_var='HIFAST_SHOW_PROG',
                        help='Show progress bar during processing.')
    parser.add_argument('--nproc', '-n', type=int, default=1, metavar='N',
                        help='Number of parallel processes to use.')

    # Baseline Fitting Method
    group = parser.add_argument_group('Baseline Fitting Method')
    
    rew_type = ['asym1', 'asym2', 'asym3', 'sym1',]
    method = ['none']
    method += ['PLS-'+r for r in rew_type]
    method += ['poly-'+r for r in rew_type]
    method += ['knpoly-'+r for r in rew_type]
    method += ['Gauss-'+r for r in rew_type]
    method += ['knspline-'+r for r in rew_type]
    method += ['masPLS-'+r for r in rew_type]
    method += ['asPLS', 'arPLS']
    method += ['original']
    
    group.add_argument('--method', default='arPLS', choices=method,
                       help='Baseline fitting algorithm. "arPLS" is robust for most cases. "poly-*" uses polynomial fitting.')
    group.add_argument('--lam', type=float, default=1.0e8, metavar='LAMBDA',
                       help='Smoothing parameter (lambda) for PLS, Spline, and Gauss methods. Larger values = smoother (stiffer) baseline. Typical range: 1e4 - 1e9.')
    group.add_argument('--deg', type=int, default=2, metavar='DEG',
                       help='Degree for polynomial or spline methods.')
    group.add_argument('--knots', type=str, metavar='JSON_FILE',
                       help='Path to JSON file with knots for spline methods.')
    group.add_argument('--offset', type=float, default=2,
                       help='Offset parameter.')
    group.add_argument('--ratio', type=float, default=0.01,
                       help='Ratio parameter.')
    group.add_argument('--niter', type=int, default=100, metavar='N',
                       help='Maximum number of iterations.')
    group.add_argument('--exclude_add', '--exclude_type', default='none', choices=['none', 'auto1', 'auto2'],
                       help='Auto-exclusion method. "auto1": exclude low weights; "auto2": exclude outliers (>3 sigma).')

    # Preprocessing
    group = parser.add_argument_group('Preprocessing (Before Fitting)')
    group.add_argument('-T', '--trans', type=bool_fun, choices=[True, False], default='False',
                       help='Fit baseline along TIME axis instead of FREQUENCY axis. WARNING: This swaps the roles of time and frequency parameters.')
    group.add_argument('--njoin', '--njoin_t', '--average_every_t', type=int, default=0, metavar='N',
                       help='Average N adjacent time samples to increase SNR before fitting. Reduces time resolution of the baseline model.')
    group.add_argument('--s_method_t', default='none', choices=['none', 'gaussian', 'boxcar', 'median'],
                       help='Smoothing method along TIME axis to suppress noise. Options: gaussian, boxcar, median.')
    group.add_argument('--s_sigma_t', type=int, default=5, metavar='SIGMA',
                       help='Width parameter for time-axis smoothing (sigma for Gaussian, window size for others).')
    group.add_argument('--s_method_freq', default='none', choices=['none', 'gaussian', 'boxcar', 'median', 'PLS', 'fft'],
                       help='Smoothing method along FREQUENCY axis.')
    group.add_argument('--s_sigma_freq', type=float, default=5, metavar='SIGMA',
                       help='Width parameter for frequency-axis smoothing.')
    group.add_argument('--average_every_freq', '--njoin_freq', type=int, default=0, metavar='N',
                       help='Binning factor along FREQUENCY axis.')

    # Exclusion
    group = parser.add_argument_group('Exclusion Settings')
    group.add_argument('--use_pre_is_excluded', type=bool_fun, choices=[True, False], default='False',
                       help='Use existing "is_excluded" mask from input file (e.g., from previous RFI flagging).')
    group.add_argument('--src_file', metavar='CATALOG',
                       help='Path to source catalog file for masking known sources (columns: RA, Dec, Radius, MinFreq, MaxFreq).')
    group.add_argument('--frame', choices=['BARYCENT', 'HELIOCEN', 'LSRK', 'LSRD'], default='LSRK',
                       help='Velocity frame for source masking.')

    # Post-processing
    group = parser.add_argument_group('Post-processing (Optional)')
    
    rew_type = ['asym1', 'asym2', 'asym3', 'sym1',]
    method = ['none']
    method += ['poly-'+r for r in rew_type]
    
    group.add_argument('--post_method', default='none', choices=method,
                       help='Method for secondary baseline fitting on individual spectra (e.g., to remove residuals after time-averaged fitting).')
    group.add_argument('--post_s_method_freq', default='none', choices=['none', 'gaussian', 'boxcar', 'median', 'PLS', 'fft'],
                       help='Smoothing method for post-processing.')
    group.add_argument('--post_s_sigma_freq', type=float, default=5, metavar='SIGMA',
                       help='Smoothing sigma for post-processing.')
    group.add_argument('--post_average_every_freq', type=int, default=0, metavar='N',
                       help='Binning factor for post-processing.')
    group.add_argument('--post_deg', type=int, default=2, metavar='DEG',
                       help='Polynomial degree for post-processing.')
    group.add_argument('--post_offset', type=float, default=2,
                       help='Offset for post-processing.')
    group.add_argument('--post_ratio', type=float, default=0.01,
                       help='Ratio for post-processing.')
    group.add_argument('--post_niter', type=int, default=100, metavar='N',
                       help='Iterations for post-processing.')
    group.add_argument('--post_exclude_add', '--post_exclude_type', default='none', choices=['none', 'auto1', 'auto2'],
                       help='Exclusion method for post-processing.')

    # Interaction
    group = parser.add_argument_group('Interaction')
    group.add_argument('-i', '--interact', action='store_true', not_in_write_out_config_file=True,
                       help='Enable interactive mode.')
    group.add_argument('--ylim', nargs='+', default=['auto'], not_in_write_out_config_file=True,
                       help='Y-axis limits for plots (e.g., "auto" or "-1 5").')
    group.add_argument('--figsize', type=float, nargs=2, default=(10, 7), not_in_write_out_config_file=True, metavar=('W', 'H'),
                       help='Figure size in inches.')
    group.add_argument('--length', type=int, default=20, not_in_write_out_config_file=True, metavar='N',
                       help='Number of spectra to test in interactive mode.')
    group.add_argument('--start_init', type=int, default=0, not_in_write_out_config_file=True, metavar='INDEX',
                       help='Starting spectrum index.')
    
    return parser

parser = create_parser()


class IO(BaseIO):
    ver = 'old'

    def _get_fpart(self,):
        fpart = '-bld'
        if self.args.post_method is not None and self.args.post_method != 'none':
            fpart += '_p'
        return fpart

    def _import_m(self,):
        """
        Lazy import of heavy dependencies (numpy, h5py, etc.).
        """
        super()._import_m()
        global h5py, OrderedDict, sub_baseline, np

        import h5py
        import numpy as np

        from collections import OrderedDict
        from .core.baseline import sub_baseline

    @staticmethod
    def fit_baseline(s2p, freq, mjd, args, is_excluded=None):
        fit_kwargs = {}
        keys = ['nproc',
                'method',
                'njoin',
                's_method_t',
                's_sigma_t',
                's_method_freq',
                's_sigma_freq',
                'average_every_freq',
                'lam',
                'deg',
                'offset',
                'ratio',
                'niter',
                'exclude_add',
               ]
        for key in keys:
            fit_kwargs[key] = getattr(args, key)
        fit_kwargs['is_excluded'] = is_excluded

        if fit_kwargs['method'].startswith('knspline') or fit_kwargs['method'].startswith('knpoly'):
            # load knots from knots file
            import json
            with open(args.knots, 'r') as f:
                fit_kwargs['knots'] = json.load(f)

        trans = getattr(args, 'trans', False)
        if trans:
            if fit_kwargs['is_excluded'] is not None:
                 fit_kwargs['is_excluded'] = fit_kwargs['is_excluded'].transpose((1, 0, 2))
            return sub_baseline(np.arange(len(mjd)), s2p.transpose((1, 0, 2)), subtract=True, **fit_kwargs).transpose((1, 0, 2))
        else:
            return sub_baseline(freq, s2p, subtract=True, **fit_kwargs)

    @staticmethod
    def post_fit_baseline(s2p, freq, mjd, args, is_excluded=None):
        fit_kwargs = {}
        keys = ['method',
                's_method_freq',
                's_sigma_freq',
                'average_every_freq',
                'deg',
                'offset',
                'ratio',
                'niter',
                'exclude_add',
               ]
        for key in keys:
            fit_kwargs[key] = getattr(args, 'post_'+key)
        fit_kwargs['nproc'] = args.nproc
        fit_kwargs['njoin'] = 0
        fit_kwargs['s_method_t'] = 'none'
        fit_kwargs['s_sigma_t'] = 0
        fit_kwargs['lam'] = 0
        fit_kwargs['is_excluded'] = is_excluded

        trans = getattr(args, 'trans', False)
        if trans:
            if fit_kwargs['is_excluded'] is not None:
                 fit_kwargs['is_excluded'] = fit_kwargs['is_excluded'].transpose((1, 0, 2))
            return sub_baseline(np.arange(len(mjd)), s2p.transpose((1, 0, 2)), subtract=True, **fit_kwargs).transpose((1, 0, 2))

        return sub_baseline(freq, s2p, subtract=True, **fit_kwargs)

    def _load_is_rfi(self,):
        fs = self.fs
        if 'is_rfi' in fs.keys():
            is_rfi = fs['is_rfi'][:]
            if self.is_use_freq is not None:
                is_rfi = is_rfi[:, self.is_use_freq]

                # M06 bandpass is not flat in [1410, 1416]
            if self.nB == 6:
                whole_rfi = np.all(is_rfi, axis=1)
                freq = self.freq
                is_use = (freq > 1406) & (freq < 1416)
                is_rfi[:, is_use] = False
                is_rfi[whole_rfi] = True

            self.is_rfi = is_rfi

    def _load_is_excluded(self):
        if not self.args.use_pre_is_excluded:
            return
        fs = self.fs
        if 'is_excluded' in fs.keys():
            is_excluded = fs['is_excluded'][:]
            if self.is_use_freq is not None:
                is_excluded = is_excluded[:, self.is_use_freq]
            # M06 bandpass is not flat in [1410, 1416]
            if self.nB == 6:
                freq = self.freq
                is_use = (freq > 1406) & (freq < 1416)
                is_excluded[:, is_use] = False
            self.is_excluded = is_excluded

    def _load_sources(self):
        args = self.args
        fpath = args.src_file
        if fpath is not None:
            from .core.regions import mask_srcs
            is_excluded = self.is_excluded if hasattr(self,'is_excluded') else np.zeros(self.s2p.shape[:2], dtype = bool)
            print(f"loading excluded files from {fpath}")
            self.is_excluded = mask_srcs(fpath, is_excluded, self.ra, self.dec, self.freq,
                                         self.mjd, inplace=True, rest_frame=args.frame)

    def gen_s2p_out(self,):
        args = self.args

        # gen self.s2p_out
        s2p = self.s2p[:]
        # rfi not in is_excluded
        self._load_is_rfi()
        # is_excluded
        self._load_is_excluded()
        # src add to is_excluded
        self._load_sources()

        # fit baseline:
        if args.method is not None and args.method != 'none':
            if hasattr(self,'is_excluded'):
                is_excluded = self.is_excluded
                # suppress nan warning in baseline fitting in the future
                # whole_rfi = np.all(is_excluded, axis=1)
                # is_excluded[whole_rfi] = False # I don't like some warnings when whole spec is nan ...
                # check is_excluded shape
                if is_excluded.ndim == 2:
                    is_excluded = np.full(is_excluded.shape + s2p.shape[2:], is_excluded[..., None])
                else:
                    raise(ValueError('Now `is_excluded` only supports 2-dim.'))
            else:
                is_excluded = None

            print('fit and substract baseline')
            s2p = self.fit_baseline(s2p, self.freq, self.mjd, args, is_excluded)
            if args.post_method is not None and args.post_method !='none':
                print('applying low-order polynomial fitting on each spectrum')
                s2p = self.post_fit_baseline(s2p, self.freq, self.mjd, args, is_excluded)
        else:
            print('no baseline method assigned, skip. Make sure the input spectra have been baselined.')
        self.s2p_out = s2p

    def __call__(self, save=True):
        self.gen_s2p_out()
        self.gen_dict_out()

        if hasattr(self, 'is_excluded'):
            # update
            self.dict_out['is_excluded'] = self.is_excluded
        if hasattr(self, 'is_rfi'):
            # update
            self.dict_out['is_rfi'] = self.is_rfi

        # save to hdf5 file
        if save:
            self.save()


def check_backend():
    import matplotlib as mpl
    if 'ipympl' not in mpl.get_backend():
        print('Please use interaction mode in Jupyert and run \'%matplotlib ipympl\' in the notebook cell first')
        sys.exit()

interact_save_lim = 20
def interact(args):
    check_backend()
    import h5py
    from .interaction import bld_i as interact
    f = h5py.File(args.fpath, 'r')
    fs = f['S']
    if 'T' in fs.keys():
        T = fs['T']
    elif 'Ta' in fs.keys():
        T = fs['Ta']
    elif 'flux' in fs.keys():
        T = fs['flux']
    if args.use_pre_is_excluded and 'is_excluded' in fs.keys():
        is_excluded = fs['is_excluded']
    else:
        is_excluded = None
    interact.T2p = T
    interact.is_excluded = is_excluded
    interact.freq = fs['freq'][:]
    interact.frange = args.frange
    interact.nproc = args.nproc
    interact.length = min(args.length, interact.T2p.shape[1])
    interact.figsize = args.figsize
    interact.ylim = args.ylim[0] if len(args.ylim) == 1 else args.ylim
    interact.trans = args.trans
    interact.start_init = args.start_init
    if T.shape[1] <= interact_save_lim:
        interact.save = True
    else:
        interact.save = False
    return interact.main()
    # sys.exit()


class IO_i(IO):

    def set_interact_spec(self, interact_spec):
        self.interact_spec = interact_spec

    def gen_s2p_out(self,):
        args = self.args
        if args.interact:
            print('Please ensure you have adjusted the parameters for each `polar`, otherwise you will get nan vaules in the file')
            self.s2p_out = self.interact_spec.bld_out
            return


if __name__ == '__main__':
    args_ = parser.parse_args()
    if args_.interact:
        interact_spec, widgets = interact(args_)[0:2]
        save = IO_i(args_, HistoryAdd={'interact':str(widgets)})
        save.set_interact_spec(interact_spec)
        if interact_spec.save:
            print('Please run \'save()\' in the notebook cell to save your results')
    else:
        print('#'*35+'Args'+'#'*35)
        args_from = parser.format_values()
        print(args_from)
        print('#'*35+'####'+'#'*35)

        HistoryAdd = {'args_from': args_from} if args_.my_config is not None else None
        io = IO(args_, HistoryAdd=HistoryAdd)
        io()
