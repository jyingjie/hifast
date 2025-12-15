

__all__ = ['IO']


import numpy as np
from .utils.io import *
from .bld import IO as bld_IO


sep_line = '##'+'#'*70+'##'
parser = ArgumentParser(prog=f"python -m hifast.{os.path.basename(sys.argv[0])[:-3]}",
                        formatter_class=formatter_class, allow_abbrev=False,
                        description='Fit and subtract baseline', )
add_common_argument(parser)
parser.add_argument('fpath',
                    help='input spectra temperature or flux file path.')
parser.add_argument('--frange', type=float, nargs=2, default=[0, float('inf')],
                    help='Limit frequence range')
parser.add_argument('--no_radec', action='store_true', not_in_write_out_config_file=True,
                    help="don't check or add ra dec")
parser.add_argument('--show_prog', type=bool_fun, choices=[True, False], default='True', env_var='HIFAST_SHOW_PROG',
                    help='')

# baseline fitting


group = parser.add_argument_group(f'*BaseLine (set --method as none to skip this) \n{sep_line}')
group.add_argument('--nproc', '-n', type=int, default=1,
                   help='number of process used in fitting baseline')

method = ['none']
method += ['MedMed', 'MinMed']

group.add_argument('--method', default='MedMed', choices=method,
                   help='method used to fit baseline, if set as none, skip this')


# MinMed or MedMed
group = parser.add_argument_group(f'*MinMed or MedMed\n{sep_line}')
group.add_argument('--npart', type=int,
                   help='divide data into n parts')
group.add_argument('--nsection', type=int,
                   help='divide each part into n sections. Same with nspec')
group.add_argument('--nspec', type=int,
                   help='divide each part into n sections. One section has n specs.')


group = parser.add_argument_group(f'*Use pre-determined is_excluded array in input file\n{sep_line}')
group.add_argument('--use_pre_is_excluded', type=bool_fun, choices=[True, False], default='True',
                   help='If True, use pre-determined is_excluded array if existed to mask channels')

# Exclude files
group = parser.add_argument_group(f'*Exclude known source catalogs\n{sep_line}')
group.add_argument('--continuum_file',
                   help='txt catalog with #ra[deg], dec[deg], R[arcmin], freq_min[MHz], freq_max[MHz]')
group.add_argument('--src_file',
                   help='txt catalog with #ra[deg], dec[deg], R[arcmin], freq_min[MHz], freq_max[MHz]')
group.add_argument('--frame', choices=['BARYCENT', 'HELIOCEN', 'LSRK', 'LSRD'], default='LSRK',
                   help='Velocity Rest Frames, HELIOCEN or LSRK')

# low-order polynomial on individual spectrum
group = parser.add_argument_group(f'*Post-applying low-order polynomial fitting on individual spectrum if `--s_method_t` or `--njoin` used [optional]\n{sep_line}')

rew_type = ['asym1', 'asym2', 'asym3', 'sym1',]
method = ['none']
method += ['poly-'+r for r in rew_type]

group.add_argument('--post_method', default='none', choices=method,
                   help='method used to fit baseline')
group.add_argument('--post_s_method_freq', default='none', choices=['none', 'median', 'gaussian', 'boxcar', 'PLS', 'fft'],
                   help='method used to smooth each spectrum along frequency axis')
group.add_argument('--post_s_sigma_freq', type=float, default=5,
                   help='used with --post_s_method_freq')
group.add_argument('--post_average_every_freq', type=int, default=0,
                   help='bin the channels by this factor')
group.add_argument('--post_deg', type=int, default=2,
                   help='polynomial degree for `*-poly` methods')
group.add_argument('--post_offset', type=float, default=2,
                   help='baseline fit parameters')
group.add_argument('--post_ratio', type=float, default=0.01,
                   help='baseline fit parameters')
group.add_argument('--post_niter', type=int, default=100,
                   help='baseline fit parameters')
group.add_argument('--post_exclude_add', '--post_exclude_type', default='none', choices=['none', 'auto1', 'auto2'],
                   help='baseline fit parameters')


class IO(bld_IO):

    def _get_fpart(self,):
        fpart = '-ref'
        if self.args.post_method is not None and self.args.post_method != 'none':
            fpart += '_p'
        return fpart

    def _import_m(self,):
        """
        need modify this function
        """
        super()._import_m()
        global h5py, OrderedDict, sub_baseline

        import h5py
        from collections import OrderedDict
        from .core.baseline import sub_baseline


    def _load_continuum(self):
        args = self.args
        fpath = args.continuum_file
        if fpath is not None:
            from .core.regions import mask_srcs

            is_excluded = np.zeros(self.s2p.shape[:2], dtype = bool)
            print(f"loading excluded files from {fpath}")
            self.is_continum = mask_srcs(fpath, is_excluded, self.ra, self.dec, self.freq,
                                         self.mjd, inplace=True, rest_frame=args.frame)

    def med_fit_baseline(self,):
        args = self.args
        # gen self.s2p_out
        s2p = self.s2p[:]
        from copy import deepcopy
        s2p_ori = deepcopy(s2p)

        is_excluded = deepcopy(self.is_excluded) if hasattr(self,'is_excluded') else np.zeros(self.s2p.shape[:2], dtype = bool)
        for key in ['is_rfi', 'is_continum']:
            if hasattr(self, key):
                is_excluded |= getattr(self, key)
        s2p[is_excluded] = np.nan

        from .ripple.running_median import minmed
        kwargs = {}
        keys = ['nsection', 'nspec', 'npart', 'method']
        for key in keys:
            kwargs[key] = getattr(args, key)
        s2p = s2p_ori - minmed(s2p, **kwargs)
        return s2p

    def gen_s2p_out(self,):
        args = self.args

        self._load_is_rfi()
        # is_excluded
        s2p = self.s2p[:]
        # is_excluded
        self._load_is_excluded()
        self._load_sources()

        # fit baseline:
        if 'Med' in args.method:
            self._load_continuum()

            print(f'Use {args.method} to substract baseline. Remember another linear substraction.')
            s2p = self.med_fit_baseline()
        else:
            print('no baseline method assigned, skip. Make sure the input spectra have been baselined.')

        if args.post_method is not None and args.post_method !='none':
            if hasattr(self,'is_excluded'):
                is_excluded = self.is_excluded
                # check is_excluded shape
                if is_excluded.ndim == 2:
                    is_excluded = np.full(is_excluded.shape + s2p.shape[2:], is_excluded[..., None])
                else:
                    raise(ValueError('Now `is_excluded` only supports 2-dim.'))
            else:
                is_excluded = None
            # not use is_continum on baseline fitting
            print('applying low-order polynomial fitting on each spectrum')
            s2p = self.post_fit_baseline(s2p, self.freq, self.mjd, args, is_excluded)
        self.s2p_out = s2p

    def __call__(self, save=True):
        self.gen_s2p_out()
        if hasattr(self,'is_continum'):
            print("save is_continum :D")
            self.gen_dict_out(is_continum = self.is_continum)
        else:
            self.gen_dict_out()
#             if 'is_continum' in self.dict_out.keys():
#                 args = self.args
#                 self.dict_out['is_continum'] = h5py.ExternalLink(os.path.relpath(
#                     args.fpath, os.path.dirname(self.fpath_out)), f'/S/is_continum')
        if hasattr(self,'is_excluded'):
            # update
            self.dict_out['is_excluded'] = self.is_excluded
        if hasattr(self, 'is_rfi'):
            # update
            self.dict_out['is_rfi'] = self.is_rfi
        # save to hdf5 file
        if save:
            self.save()



if __name__ == '__main__':
    args_ = parser.parse_args()

    print('#'*35+'Args'+'#'*35)
    args_from = parser.format_values()
    print(args_from)
    print('#'*35+'####'+'#'*35)

    HistoryAdd = {'args_from': args_from} if args_.my_config is not None else None
    io = IO(args_, HistoryAdd=HistoryAdd)
    io()
