

__all__ = ['IO', 'IO_i']


from .utils.io import *
from copy import deepcopy


sep_line = '##'+'#'*70+'##'
parser = ArgumentParser(prog=f"python -m hifast.{os.path.basename(sys.argv[0])[:-3]}",
                        formatter_class=formatter_class, allow_abbrev=False,
                        description='Fit and subtract standing wave',)
add_common_argument(parser)
parser.add_argument('fpath',
                    help='input baselined spectra file path')
parser.add_argument('--frange', type=float, nargs=2, default=[0, float('inf')],
                    help='Limit frequence range')
parser.add_argument('--no_radec', action='store_true', not_in_write_out_config_file=True,
                    help="if set, don't check or add ra dec")
parser.add_argument('--nproc', '-n', type=int, default=1,
                    help='number of process used in fitting baseline')
parser.add_argument('--show_prog', type=bool_fun, choices=[True, False], default='True', env_var='HIFAST_SHOW_PROG',
                    help='show progress bar')

# smooth
group = parser.add_argument_group('For FFT, smooth to find where should be replaced; for sin-fitting, preprocess.')
group.add_argument('--s_method_t', default='gaussian', choices=['none', 'gaussian', 'boxcar', 'median'],
                   help='smooth method along time axis')
group.add_argument('--s_sigma_t', type=int, default=5,
                   help='smooth sigma along time axis, size = (2*sigma+1) for boxcar and median')
group.add_argument('--s_method_freq', default='gaussian', choices=['none', 'median', 'gaussian', 'boxcar'],
                   help='smooth method along freq axis')
group.add_argument('--s_sigma_freq', type=float, default=3,
                   help='smooth sigma along freq axis, size = (2*sigma+1) for boxcar and median')

group = parser.add_argument_group('For FFT, smooth for finding Troughs of standing waves')
group.add_argument('--s_method_t_T', default='gaussian', choices=['none', 'gaussian', 'boxcar', 'median'],
                   help='smooth method along time axis')
group.add_argument('--s_sigma_t_T', type=int, default=10,
                   help='smooth sigma along time axis, size = (2*sigma+1) for boxcar and median')
group.add_argument('--s_method_freq_T', default='gaussian', choices=['none', 'median', 'gaussian', 'boxcar'],
                   help='smooth method along freq axis')
group.add_argument('--s_sigma_freq_T', type=float, default=5,
                   help='smooth sigma along freq axis, size = (2*sigma+1) for boxcar and median')

# method
group = parser.add_argument_group(f'*select method used to fit standing wave\n{sep_line}')
group.add_argument('--method', default='fft', choices=['sin_poly', 'fft', 'running_median', 'running_mean'],
                   help='method to fit standing wave')
group.add_argument('--nobld', type=bool_fun, choices=[True, False], default='False',
                   help="if True, use the spectra before bld to subtract standing wave and output file name will add 'nobld'")
group.add_argument('--fpattern_nobld',
                   help='if not specify, guess from the History recored in the fpath')

# least square
group = parser.add_argument_group(f'*parameters for --method sin_poly\n{sep_line}:' +
                                  '\n'+'addition options for preprocessing before fitting')
group.add_argument('--njoin', type=int, default=0,
                   help='join nspec along t')
group = parser.add_argument_group('sinusoidal standing wave (sw) fitting')
group.add_argument('--exclude_m', type=int, default=0,
                   help='')
group.add_argument('--sin_f', type=float, nargs='+', default=0.929,
                   help='Initial guess for sw freq')
group.add_argument('--bound_f', type=float, nargs=2, default=[.90, .95],
                   help='Lower and upper bounds for sw freq')
group.add_argument('--deg', type=int, default=1,
                   help='Degree of the added polynomial')
# fft
group = parser.add_argument_group(f'*parameters for --method fft\n{sep_line}:' +
                                  '\n'+'replace big RFI or sources firstly')
group.add_argument('--iter_twice', type=bool_fun, choices=[True, False], default='True',
                   help='iterate twice? first remove the ripples and then back to find replace area again.')
# replace big RFI or sources firstly
group.add_argument('--rfi_method', default='near_ripple',
                   choices=['near_ripple','zero_ripple',],
                   help='method to replace big RFI, recommend the first one')
group.add_argument('--mw_frange', type=float, nargs=2, default=[0, 0],
                   help='milky way freq range. If do not use, keep it as default.')
group.add_argument('--rms_sigma', type=float, default=6,
                   help='gauss filter sigma to compute real rms')
group.add_argument('--rms_frange', type=float, nargs=2, default=[0, 0],
                   help='freq range to compute rms')
group.add_argument('--rms_step', type=float, default=5,
                   help='a step (MHz) to find where to compute real rms')

group.add_argument('--times_thr', type=float, default=4,
                   help='sparks above ~ times of rms will be set noise (unsmooth)')
group.add_argument('--times_s_thr', type=float, default=3,
                   help='above ~ times of rms will be replaced (smoothed once)')
group.add_argument('--times_s_thr2', type=float, default=2,
                   help='above ~ times of rms will be replaced (smoothed twice)')
group.add_argument('--ext_freq', type=float, default=1.3,
                   help='extend freq range to replace (mhz)')
group.add_argument('--rfi_width_lim', type=float, default=20,
                   help='rfi should contain more channels than limit')
group.add_argument('--ext_sec', type=int, default=20,
                   help='extend channel number of start and end of each section')
group.add_argument('--restrict_bound', type=bool_fun, choices=[True, False], default='False',
                   help='restrict replace area bound')
group.add_argument('--save_is_excluded', type=bool_fun, choices=[True, False], default='False',
                   help='save the replaced area as is_excluded')

group = parser.add_argument_group(f'phase space')
# remove which component in fft? (only 4 types now)
group.add_argument('--sw_base', type=bool_fun, choices=[True, False], default='True',
                   help='if True, remove constant components')
group.add_argument('--sw_periods', nargs='+', choices=['1mhz', '2mhz', '0_04mhz', 'none'], default=['1mhz',  '0_04mhz'],
                   help='remove ripple (1mhz: 1.08mhz, 2mhz:1.92mhz, 0_04mhz: 0.039 mhz)')
group.add_argument('--check_2mhz', type=bool_fun, choices=[True, False], default='False',
                   help='if True, remove 2mhz from sw_periods except for Beam 6')

group.add_argument('--amp_thr_mean_factor', type=float, nargs=2, default=[1.05, 1.3],
                   help='above mean amptitude threshold will be chosed, noise off and on')
group.add_argument('--amp_thr_solo_factor', type=float, nargs=2, default=[1.4, 1.7],
                   help='above amptitude threshold will be chosed in every spec, noise off and on')
group.add_argument('--chan_wide', type=int, default=5,
                   help='channel numbers near 1mhz to be chosed (wide)')
group.add_argument('--chan_narr', type=int, default=2,
                   help='channel numbers near 1mhz to be chosed (narrow)')
group.add_argument('--choose_method', default='all', choices=['all', 'interpolate'],
                   help='method to choose components in fft')

# running median or mean
group = parser.add_argument_group(f'*parameters for --method running_median or running_mmean\n{sep_line}')
group.add_argument('--nspec',type=int, default=200,
                    help='average how many specs to fit baseline.')
group.add_argument('--func', default='iter', choices=['iter', 'smooth'],
                       help='function to mean or median')

# Exclude files
group = parser.add_argument_group(f'Exclude known source catalogs\n{sep_line}')
group.add_argument('--src_file',
                   help='Path to a text catalog with columns for RA [deg], Dec [deg], radius [arcmin], minimum frequency [MHz], and maximum frequency [MHz]. Identified sources will be added to the is_excluded array.')
group.add_argument('--frame', choices=['BARYCENT', 'HELIOCEN', 'LSRK', 'LSRD'], default='LSRK',
                   help='Specify the velocity rest frame for excluding source frequency ranges.')

# interaction
group = parser.add_argument_group(f'*Interaction\n{sep_line}')
group.add_argument('-i', '--interact', action='store_true', not_in_write_out_config_file=True,
                   help='interaction')
group.add_argument('--ylim', nargs='+', default=['auto'], not_in_write_out_config_file=True,
                   help='ylim')
group.add_argument('--figsize', type=float, nargs=2, default=(10, 7), not_in_write_out_config_file=True,
                   help='figsize')
group.add_argument('--length', type=int, default=40, not_in_write_out_config_file=True,
                   help='spetra numbers used to test')
group.add_argument('--start_init', type=int, default=0, not_in_write_out_config_file=True,
                   help='the index of the spectrum shown at start')


class IO(BaseIO):
    ver = 'old'

    def __init__(self, *args, dict_bef_bld=None, **kwargs):
        self.dict_bef_bld = dict_bef_bld
        super().__init__(*args, **kwargs)

    def _get_fpart(self,):
        if self.args.nobld:
            return '-sw_nobld'
        return '-sw'

    def _import_m(self,):
        """
        need modify this function
        """
        super()._import_m()
        global h5py, np, OrderedDict, sub_baseline, get_exclude_fun
        import h5py
        import numpy as np
        from collections import OrderedDict
        from .core.baseline import sub_baseline, get_exclude_fun

    @staticmethod
    def find_fname_his(Header, fname):
        """
        from historys in Header to find file:
        return: [path, path related to the 'outdir' in history]
        """
        import json
        fpath_in_list = []
        outdir_list = []
        cwd_list = []
        for key in Header.keys():
            if key.startswith('HISTORY'):
                his = json.loads(Header[key])
                try:
                    args_h = json.loads(his['args'])
                    fpath_in_list += [args_h['fpath']]
                    outdir_list += [args_h['outdir']]
                    cwd_list += [his['cwd']]
                except:
                    try:
                        import re
                        fpath_in_list += [re.findall(r'fpath=.*.hdf5', his['args'])[0].split('=')[1][1:]]
                        outdir_list += [re.findall(r'outdir=[^,]*', his['args'])[0].split('=')[1][1:-1]]
                        cwd_list += [his['cwd']]
                    except:
                        pass
        for i, val in enumerate(fpath_in_list):
            if os.path.basename(val) == fname:
                path_ = os.path.abspath(os.path.join(cwd_list[i], val))
                out = [path_, os.path.relpath(path_, os.path.join(cwd_list[i], outdir_list[i]))]
                break
        return out

    def _gen_s2p_ori_fpath(self,):
        """
        find file path of spectra before baseline subtracted
        """
        args = self.args
        if args.fpattern_nobld is not None:
            print('using the path from --fpattern_nobld as the spectra file before baseline subtracted')
            ori_fpath = self.replace_nB(args.fpattern_nobld, self.nB)
            if not os.path.exists(ori_fpath):
                raise(FileNotFoundError(f'file {ori_fpath} not exists'))

        else:
            print('guess the path of the spectra file before baseline subtracted from the History')
            fname_find = os.path.basename(args.fpath).rsplit('bld')[0][:-1] + '.hdf5'
            ori_fpaths = self.find_fname_his(self.Header, fname_find)
            if os.path.exists(ori_fpaths[0]):
                ori_fpath = ori_fpaths[0]
            elif os.path.exists(ori_fpaths[1]):
                ori_fpath = ori_fpaths[1]
            else:
                raise(FileNotFoundError("can't find the spectra file before baseline subtracted, \
                                        please specify --fpattern_nobld"))
        self.s2p_ori_fpath = ori_fpath
        self.Header['ori_fpath'] = self.s2p_ori_fpath

    def _load_s2p_ori(self):
        """
        load spectra before baseline subtracted
        """
        if self.dict_in is None:
            self._gen_s2p_ori_fpath()
            print(f'load the spectra before baseline subtracted from: \n{self.s2p_ori_fpath}')
            f = h5py.File(self.s2p_ori_fpath, 'r')['S']
        else:
            if self.dict_bef_bld is None:
                raise(ValueError('need input dict_bef_bld'))
            f = self.dict_bef_bld
        freq = f['freq'][:]
        s2p_ori = f[self.outfield][:]
        is_use_freq = (freq >= self.freq.min()) & (freq <= self.freq.max())
        if not np.all(is_use_freq):
            inds = np.where(is_use_freq)[0]
            s2p_ori = s2p_ori[..., inds[0]:inds[-1]+1]  # freq axis at end
        s2p_ori = PolarMjdChan_to_MjdChanPolar(s2p_ori)
        return s2p_ori


    def _load_sources(self):
        args = self.args
        fpath = args.src_file
        if fpath is not None:
            from .core.regions import mask_srcs
            is_excluded = self.is_excluded if hasattr(self,'is_excluded') else np.zeros(self.s2p.shape[:2], dtype = bool)
            print(f"loading excluded files from {fpath}")
            self.is_excluded = mask_srcs(fpath, is_excluded, self.ra, self.dec, self.freq,
                                         self.mjd, inplace=True, rest_frame=args.frame)


    def _load_is_rfi(self):
        # load is_rfi, is_on
        dict_ = self.fs if self.dict_in is None else self.dict_in
        is_rfi = dict_['is_rfi'][:] if 'is_rfi' in dict_.keys() else None
        if is_rfi is not None and self.is_use_freq is not None:
            is_rfi = is_rfi[:, self.is_use_freq]
        self.is_rfi = is_rfi
        self.is_on = dict_['is_on'][:]


    def fit_sw(self, s2p, freq, args, subtract):
        """
        sin poly
        """
        fit_kwargs = {}
        keys = ['method',
                'njoin', 's_method_t', 's_sigma_t', 's_method_freq', 's_sigma_freq',
                'deg', 'sin_f']
        for key in keys:
            fit_kwargs[key] = getattr(args, key)
        fit_kwargs['nproc'] = args.nproc
        fit_kwargs['rew'] = False
        fit_kwargs['niter'] = 1
        if args.method == 'sin_poly':
            # sin_poly: optimize.minimize
            bounds = [(0., 1.), args.bound_f, (0, 2*np.pi), (-1, 1)] + [(-np.inf, np.inf), ]*args.deg
        opt_para = {'bounds': bounds, }
        fit_kwargs['exclude_fun'] = get_exclude_fun(args.exclude_m)
        m, f = s2p.shape[0], s2p.shape[1]
        fit_kwargs['is_excluded'] = np.full((2, m, f), self.is_excluded).transpose(1, 2, 0) if hasattr(self, 'is_excluded') else None
        return sub_baseline(freq, s2p, subtract=subtract, opt_para=opt_para, **fit_kwargs)

    def check_rms_range(self, is_rfi, rms_step=5):
        args = self.args
        rms_frange = args.rms_frange
        if (rms_frange is None) or (rms_frange[1] - rms_frange[0] <= 0):
            from .ripple.markRFI import get_rms_frange
            if is_rfi is not None:
                is_lf = np.all(is_rfi, axis = 1)
                is_excluded = np.any(is_rfi[~is_lf], axis = 0)
                data = self.s2p[~is_lf]
            else:
                data = self.s2p
                is_excluded = None
            spec = np.nanmean(np.nanmean(data, axis = 0),axis = -1)
            self.args.rms_frange = get_rms_frange(spec,self.freq,rms_step=5,is_excluded = is_excluded)
        else:
            freq = self.freq
            if (rms_frange[0] < freq[0]) or (rms_frange[1] > freq[-1]):
                raise ValueError(f"rms frange {rms_frange} is not in freq range {[freq[0],freq[-1]]}.")

    def fft_fit_sw(self, s1p, is_rfi=None, is_on=None,iter_twice = False,
                  sm_find = None, sm_trough = None, sm_res = None, ):
        """
        fft
        """
        from .ripple import sw_fft
        args = self.args

        # replace rfi and mw
        # replace args
        rep_args = {}
        self.check_rms_range(is_rfi, rms_step = args.rms_step)
        if args.rfi_method in ['near_ripple','zero_ripple',]:
            keys = ['rms_sigma', 'rms_frange', 'times_s_thr','times_s_thr2',
                    'times_thr', 'rfi_width_lim', 'ext_sec', 'ext_freq', 'mw_frange',]
        else:
            raise(ValueError('not supported rfi_method'))
        for key in keys:
            rep_args[key] = getattr(args, key)

        if not iter_twice: rep_args['times_s_thr2'] = None
        if args.rfi_method == 'zero_ripple':
            rep_args['fill'] = 'zero'
            args.rfi_method = 'near_ripple'
        rep_args['verbose'] = args.show_prog

        save_is_excluded = args.save_is_excluded & iter_twice
#         print("############ save_is_excluded", save_is_excluded)
        s1p, is_excluded = sw_fft.replace_rfi(s1p, self.freq, time_rfi=is_rfi, method=args.rfi_method,
                                  data_find = sm_find, data_trough = sm_trough, data_restrict = sm_res,
                                 save_is_excluded = save_is_excluded, **rep_args)
        # fft args
        fft_args = {}
        keys = ['chan_wide', 'chan_narr',  'choose_method',
                'choose_method', 'sw_periods', 'sw_base']
        for key in keys:
            fft_args[key] = getattr(args, key)
        if args.check_2mhz:
            print('checking 2mhz')
            try:
                fft_args['sw_periods'].remove('2mhz')
            except:
                pass
            if self.nB in [6, ]:
                fft_args['sw_periods'].insert(-1, '2mhz')
            else:
                print('beam number !=6, remove 2mhz in sw_periods if it exists')

        # use the first factor for noise off or not defined
        for key in ['amp_thr_mean_factor', 'amp_thr_solo_factor']:
            fft_args[key] = getattr(args, key)[0]

        if is_on is not None:
            ret = np.zeros_like(s1p)
            print("# noise off")
            fft_args['is_excluded_mean'] = np.all(is_rfi[~is_on], axis=1) if is_rfi is not None else None
            ret[~is_on] = sw_fft.fit_sw_fft(s1p[~is_on], self.freq, args.nproc, **fft_args)

            print("# noise on")
            for key in ['amp_thr_mean_factor', 'amp_thr_solo_factor']:
                fft_args[key] = getattr(args, key)[1]
            fft_args['is_excluded_mean'] = np.all(is_rfi[is_on], axis=1) if is_rfi is not None else None
            ret[is_on] = sw_fft.fit_sw_fft(s1p[is_on], self.freq, args.nproc, **fft_args)
        else:
            fft_args['is_excluded_mean'] = np.all(is_rfi, axis=1) if is_rfi is not None else None
            ret = sw_fft.fit_sw_fft(s1p, self.freq, args.nproc, **fft_args)
        return ret, is_excluded

    def med_fit_sw(self, s1p,):
        """
        running median or mean
        """
        from .ripple.running_median import running_median
        args = self.args

        fit_args = {}
        keys = ['nspec','func', 'nproc']
        for key in keys:
            fit_args[key] = getattr(args, key)
        fit_args['method'] = args.method[8:]

        if hasattr(self, 'is_excluded'):
            s1p[self.is_excluded | self.is_rfi] = np.nan

#         sw_on = running_median(s1p[self.is_on], **fit_args)
#         sw_off = running_median(s1p[~self.is_on], **fit_args)

#         sw = np.zeros_like(s1p)
#         sw[self.is_on] = sw_on
#         sw[~self.is_on] = sw_off
        sw = running_median(s1p, **fit_args)

        return sw


    def gen_s2p_out(self,):
        args = self.args
        # gen self.s2p_out
        s2p = self.s2p[:]

        # src add to is_excluded
        self._load_sources()
        self._load_is_rfi()

        # fit baseline:
        print(f'standing wave fitting and substract by {args.method}')
        if args.method in ['sin_poly', ]:
            if args.nobld:
                self.s2p_out = self._load_s2p_ori()[:] - self.fit_sw(s2p, self.freq, args, subtract=False)
            else:
                self.s2p_out = self.fit_sw(s2p, self.freq, args, subtract=True)
        else:
            is_rfi = self.is_rfi
            is_on = self.is_on

            s2p_out = deepcopy(s2p)

            if args.method == 'fft':
                if args.iter_twice:
                    s2p_in = s2p
                else:
                    s2p_in = self._load_s2p_ori()[:] if args.nobld else s2p

                # smooth to find areas that need replacement
                find_kwargs = {}
                keys = ['s_method_t', 's_sigma_t', 's_method_freq', 's_sigma_freq',]
                for key in keys:
                    find_kwargs[key] = getattr(args, key)
                find_kwargs['is_rfi'] = is_rfi
                find_kwargs['is_on'] = is_on
                from .ripple.util import do_smooth_onoff
                sm_find = do_smooth_onoff(s2p, **find_kwargs)
                print("------------------")

                # smooth to find troughes
                if args.rfi_method != 'zero_ripple':
                    trough_kwargs = {}
                    keys = ['s_method_t_T', 's_sigma_t_T', 's_method_freq_T', 's_sigma_freq_T',]
                    for key in keys:
                        trough_kwargs[key[:-2]] = getattr(args, key)
                    trough_kwargs['is_rfi'] = is_rfi
                    trough_kwargs['is_on'] = is_on
                    if trough_kwargs == find_kwargs:
                        sm_trough = sm_find
                    else:
                        sm_trough = do_smooth_onoff(s2p, **trough_kwargs)
                        print("------------------")
                else:
                    sm_trough = sm_find
                    trough_kwargs = {}
                    trough_kwargs['s_method_t'] = 'none'

                # smooth to restrict bounds
                if trough_kwargs['s_method_t'] != 'none' and args.restrict_bound:
                    res_kwargs = {}
                    res_kwargs['s_method_freq'] = 'gaussian'
                    res_kwargs['s_sigma_freq'] = args.rms_sigma
                    sm_res = do_smooth_onoff(s2p, **res_kwargs)
                    print("------------------")
                else:
                    sm_res = np.array([None, None])[None,None,:]

                for i in range(s2p.shape[2]):
                    s2p_out[..., i] = s2p_in[..., i] - self.fft_fit_sw(s2p[..., i], is_rfi, is_on,
                          sm_find = sm_find[..., i], sm_trough = sm_trough[..., i], sm_res = sm_res[..., i])[0]

                if args.iter_twice:
                    print("Second iter ...")
                    # smooth
                    s2p_in = self._load_s2p_ori()[:] if args.nobld else s2p
                    sm_find = do_smooth_onoff(s2p_out, **find_kwargs)
                    print("------------------")
                    if trough_kwargs['s_method_t'] != 'none' and args.restrict_bound:
                        sm_res = do_smooth_onoff(s2p_out, **res_kwargs)
                        print("------------------")
                    else:
                        sm_res = np.array([None, None])[None,None,:]

                    is_excluded = self.is_excluded if hasattr(self,'is_excluded') else np.zeros_like(s2p,dtype=bool)
                    sw2 = deepcopy(s2p)
                    for i in range(s2p.shape[2]):
                        sw2[..., i], is_excluded[..., i] = self.fft_fit_sw(s2p[..., i], is_rfi, is_on, iter_twice = True,
                               sm_find = sm_find[..., i], sm_trough = sm_trough[..., i], sm_res = sm_res[..., i])
                        s2p_out[..., i] = s2p_in[..., i] - sw2[..., i]

                    if args.save_is_excluded:
                        self.is_excluded = np.any(is_excluded,axis = -1)
                        if is_rfi is not None:
                            self.is_excluded |= is_rfi


            elif (args.method == 'running_median') or (args.method == 'running_mean'):
                s2p_out = s2p - self.med_fit_sw(s2p_out, )

            self.s2p_out = s2p_out


    def __call__(self, save=True):
        self.gen_s2p_out()
        if hasattr(self,'is_excluded'):
            print("save is_excluded :D")
            self.gen_dict_out(is_excluded = self.is_excluded)
        else:
            self.gen_dict_out()
        # save to hdf5 file
        if save:
            self.save()


def check_backend():
    import matplotlib as mpl
    if 'ipympl' not in mpl.get_backend():
        print('Please run interaction in Jupyert and \'%matplotlib ipympl\' in the notebook cell ')
        sys.exit()

interact_save_lim = 20
def interact(args):
    check_backend()
    import h5py
    from .interaction import sw_i as interact
    f = h5py.File(args.fpath, 'r')
    fs = f['S']
    if 'T' in fs.keys():
        T = fs['T']
    elif 'Ta' in fs.keys():
        T = fs['Ta']
    elif 'flux' in fs.keys():
        T = fs['flux']
    interact.T2p = T
    interact.freq = fs['freq'][:]
    interact.frange = args.frange
    interact.nproc = args.nproc
    interact.length = min(args.length, interact.T2p.shape[1])
    interact.figsize = args.figsize
    interact.ylim = args.ylim[0] if len(args.ylim) == 1 else args.ylim
    interact.start_init = args.start_init
    if T.shape[1] <= interact_save_lim:
        interact.save = True
    else:
        interact.save = False
    return interact.main()


class IO_i(IO):
    def gen_s2p_out(self,):
        args = self.args
        if args.interact:
            print('Please ensure you have adjusted the parameters for each `polar`, otherwise you will get nan vaules in the file')
            self.s2p_out = interact_spec.bld_out
            return


if __name__ == '__main__':
    args_ = parser.parse_args()
    # print(parser.format_help())
    # print("----------")
    # print(parser.format_values())  # useful for logging where different settings came from

    if args_.interact:
        interact_spec, widgets = interact(args_)[0:2]
        save = IO_i(args_, HistoryAdd={'interact':str(widgets)})
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
