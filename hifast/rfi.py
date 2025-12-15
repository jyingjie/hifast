

__all__ = ['IO']


from .utils.io import *
from copy import deepcopy


sep_line = '##'+'#'*70+'##'
parser = ArgumentParser(prog=f"python -m hifast.{os.path.basename(sys.argv[0])[:-3]}",
                        formatter_class=formatter_class, allow_abbrev=False,
                        description='rfi', )
add_common_argument(parser)
parser.add_argument('fpath',
                    help='input baselined spectra file path.')
# not support --frange
parser.add_argument('--no_radec', action='store_true', not_in_write_out_config_file=True,
                    help="don't check or add ra dec")
parser.add_argument('--show_prog', type=bool_fun, choices=[True, False], default='True', env_var='HIFAST_SHOW_PROG',
                    help='')
parser.add_argument('--all_beams', type=bool_fun, choices=[True, False], default='False',
                   help='find time rfi after averaging all 19 beams')

parser.add_argument('--replace_rfi', type=bool_fun, choices=[True, False], default='False',
                   help='If True, replace find rfi as np.nan')
                    #, otherwise store ``is_rfi`` array in output file')

parser.add_argument('--nobld', type=bool_fun, choices=[True, False], default='False',
                   help="if True, output the spectra before bld and output file name will add 'nobld'")

group = parser.add_argument_group(f'*mark rfi from a DS9 regions file\n{sep_line}')
group.add_argument('--reg_from', default='none',
                   help='If not set as ``none``, mark rfi through regions from a DS9 format regions file. ' + \
                        "If set as ``shared``, beams will shared one region file" + \
                        "If set as ``default``, will try to find the file:  the input spectra  fpath + '.reg'. If set as other string, will be treated as a file path")
group.add_argument('--reg_shared_beams', default='all',
                   help="If ``--reg_from`` set as ``shared``, using the region on these beams.")

#################### Time domain continuous RFI ######################
group = parser.add_argument_group(f'*Time domain continuous RFI\n{sep_line}')
group.add_argument('--tr', '--time_rfi', type=bool_fun, choices=[True, False], default='False',
                   help='')
group.add_argument('--tr_s_method_t', default='gaussian',
                   help='smoothing method along time axis')
group.add_argument('--tr_s_sigma_t', type=int, default=10,
                   help='')
group.add_argument('--tr_s_method_freq', default='gaussian',
                   help='smoothing method along time freq')
group.add_argument('--tr_s_sigma_freq', type=float, default=5,
                   help='')
group.add_argument('--tr_n_continue', type=int, default=100,
                   help='')
group.add_argument('--tr_times', type=float, default=6.,
                   help='')
group.add_argument('--tr_times_s', type=float, default=1.5,
                   help='')
group.add_argument('--ext_add', type=int, default=0,
                   help='extend rfi range')
group.add_argument('--ext_frac', type=float, default=0.,
                   help='between 0 and 1, extend rfi range')
# parser.add_argument('--cross_frac', type=float, default=0.,
#                    help='between 0 and 1')

####################### Narrowband single channel RFI #########################
group = parser.add_argument_group(f'*Narrowband single channel RFI\n{sep_line}')
group.add_argument('--nr', '--narr_rfi', type=bool_fun, choices=[True, False], default='False',
                   help='find long time narrow single channel rfi')
group.add_argument('--nr_mean_times', type=float, default=100,
                   help="first threhold, rfi is this times of median value after mean along time axis; \
                   ")
group.add_argument('--nr_diff_times', type=float, default=30,
                   help="second threhold, sharp edge on time axis. diff above this times of median will be recognized; \
                   ")
group.add_argument('--nr_rfi_width_lim',type=int, nargs=2, default= [0, 2],
                   help='rfi channel width limit')
group.add_argument('--nr_mask_rms_times',type = float, default=0,
                   help='if == 0, mask whole channel; if > 0, mask RMS above ~ times of RMS. ')

####################### Polarized RFI #######################################
group = parser.add_argument_group(f'*Polarized RFI\n{sep_line}')
group.add_argument('--pr', '--polar_rfi', type=bool_fun, choices=[True, False], default='False',
                   help='')
group.add_argument('--pr_s_sigma', type=float, default=5,
                   help='gaussian smooth size for spectra smoothing along time axis')
group.add_argument('--pr_times', type=float, default=6,
                   help='')
group.add_argument('--pr_times_s', type=float, default=1,
                   help='')

####################### 19Beam RFI #######################################
group = parser.add_argument_group(f'*19Beam RFI\n{sep_line}')
group.add_argument('--b19', '--beam19_rfi', type=bool_fun, choices=[True, False], default='False',
                   help='Enable or disable 19-beam RFI detection')
group.add_argument('--b19_times', type=float, default=3,
                   help='Threshold multiplier for 19-beam RFI detection')
group.add_argument('--b19_min_percent', type=float, default=0.23,
                   help='Minimum percent of beams that must exceed the threshold to classify as RFI, default >0.23 (4/19 beam)')

###################### Time domain uncontinuous RFI #########################
parser.add_argument('--rms_frange', type=float, nargs=2,
                   help='freq range to compute rms, NEED TO DEFINE when sf, pdr')
parser.add_argument('--mw_frange', type=float, nargs=2,
                   help='milky way freq range')
## long freq time rfi
group = parser.add_argument_group(f'*Long freq \n{sep_line}')
group.add_argument('--lf', '--long_freq', type=bool_fun, choices=[True, False], default='False',
                   help='find time rfi')
group.add_argument('--lsn_thr_type', default='input_absmed_times',
                   choices=['input_med_times','input_absmed_times','input_posimed_times'],
                   help='input times of median value, its absolute value, or add an offset to make it positive.\
                   used for lf, sf, nr')
group.add_argument('--lf_frange', type=float, nargs=2, default = [1400, 1500],
                   help='freq range exists long-freq time rfi')
group.add_argument('--lf_mean_times', type=float, default=2,
                   help="first threhold, rfi is this times of median value after mean along time axis; \
                   ")
group.add_argument('--lf_diff_times', type=float, default=0,
                   help='set 0 and do not change this parameter')
group.add_argument('--lf_rfi_last',type=float, nargs=2, default= [50, float("INF")],
                   help='rfi lasts at least ~ spec numbers')
group.add_argument('--lf_ext_add',type = int,default=50,
                   help='extend edge')
group.add_argument('--lf_mask_rms_times',type = float, default=-1,
                   help='if == -1, mask whole spec; if == 0, only mask region in frange; \
                   if > 0, mask RMS above ~ times of RMS. ')
## short freq time rfi
group = parser.add_argument_group(f'*Short freq \n{sep_line}')
group.add_argument('--sf', '--short_freq', type=bool_fun, choices=[True, False], default='False',
                   help='find time rfi')
group.add_argument('--sf_use_time_only', type=bool_fun, choices=[True, False], default='False',
                   help="if true, only use 'is_timerfi' of 19 beams")
group.add_argument('--sf_frange', type=float, nargs=2, default = [1378, 1385],
                   help='freq range exists short-freq time rfi')
group.add_argument('--sf_frange_step',type = int,
                   help='if sf_frange is None and sf_file is None, cycle in whole freq band.')
group.add_argument('--sf_file',
                   help='freq range exists short-freq time rfi npy filename')
group.add_argument('--sf_mean_times', type=float, default=3,
                   help="first threhold, rfi is this times of median value after mean along time axis; \
                   ")
group.add_argument('--sf_diff_times', type=float, default=.08,
                   help="second threhold, sharp edge on time axis. diff above this times of median will be recognized; \
                   ")
group.add_argument('--sf_rfi_last',type=float, nargs=2, default=[20,float("INF")],
                   help='rfi lasts at least 10 spec numbers')
group.add_argument('--sf_ext_add',type = int,default=3,
                   help='extend edge')
group.add_argument('--sf_mask_rms_times',type = float, default=2,
                   help='mask from peak to 2 sides, until RMS drops to 2 times of RMS')

################## invoking aoflagger #####################
group = parser.add_argument_group(f'*Invoke aoflagger \n{sep_line}')
group.add_argument('--aoflagger', '--af', type=bool_fun, choices=[True, False], default='False',
                   help='Using aoflagger to mask rfi')
group.add_argument('--af_strategy_file', '--afsf', type=str,
                   help='path of strategy file')

################## Period 8 MHZ RFI #######################
parser.add_argument('--rms_sigma', type=float, default =6,
                   help='gauss filter sigma to compute real rms')
## smooth
group = parser.add_argument_group('preprocessing before finding pdr')
group.add_argument('--s_method_t', default='none', choices=['none', 'gaussian', 'boxcar', 'median'],
                   help='')
group.add_argument('--s_sigma_t', type=int, default=5,
                   help='')
group.add_argument('--s_method_freq', default='none', choices=['none', 'median', 'gaussian', 'boxcar'],
                   help='')
group.add_argument('--s_sigma_freq', type=float, default=5,
                   help='')
## find rfi
group = parser.add_argument_group(f'*find Period 8MHz \n{sep_line}')
group.add_argument('--pdr', '--period_rfi', type=bool_fun, choices=[True, False], default='False',
                   help='find period 8 MHz rfi')
group.add_argument('--rfi_thr', type=float, default=3,
                   help='default 3 times of rms threshold')
group.add_argument('--rfi_width_lim', type=float, default=10,
                   help='rfi should contain more channels than limit')
group.add_argument('--ext_sec', type=int, default=30,
                   help='extend channel number of start and end of each section')
group.add_argument('--freq_thr', type=float, default=.5,
                   help='estimate error when looking for peaks to polyfit')
group.add_argument('--freq_step', type=float, default = 8.1,
                   help='rfi period MHz')
group.add_argument('--rfi_groups', default = 'two_groups', choices=['two_groups','three_groups','all'],
                    help='divide rfi into 2 or 3 groups')

## mask rfi
group = parser.add_argument_group('mask 8 Mhz RFI')
group.add_argument('--mask_RFI_method',default='fixed_freq',choices=['2_sides','fixed_freq'],
                     help='from center to two sides, or use a fixed freq width')
group.add_argument('--freq_from_theory', type=float, default=.5,
                    help='mask width from theory center')
group.add_argument('--mask_thr', type=float, default=.7,
                   help=' above ~ times of rms threshold will be masked')
group.add_argument('--ext_edge', type=int, default= 0,
                   help=' extend result channel on freq axis')

### 2 sides
group.add_argument('--small_rfi_times',type=float, default=2,
                    help='small rfi below, eg.2*RMS will not be masked')
group.add_argument('--chan_step' ,type=int, default=5,
                    help='channel step when walk from center to two sides')
### fixed freq
group.add_argument('--mask_all_theory', action= 'store_true',
                   help='mask_all_theory')
### time coherent
group.add_argument('--time_coherent_per', type=float, default = 0.7,
                   help='rfi in one freq appears more than emmm, maybe 0.7, mask them all on time axis.')


class IO(BaseIO):
    ver = 'old'
    def _get_fpart(self,):
        if self.args.nobld:
            return '-rfi_nobld'
        return '-rfi'

    def _import_m(self,):
        """
        need modify this function
        """
        super()._import_m()
        global h5py, OrderedDict, np

        import h5py
        from collections import OrderedDict
        import numpy as np


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
        if hasattr(args, 'fpattern_nobld') and args.fpattern_nobld is not None:
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

    def get_from_regions(self,):

        from .core.regions import read_regions, replace_region
        args = self.args

        if args.reg_from is None or args.reg_from == 'none':
            return None

        nB = get_nB(args.fpath)
        if args.reg_from == 'default':
            print('try to find the default regions file')
            fpath_reg = args.fpath + '.reg'
            if not os.path.exists(fpath_reg):
                print(f'can not find the default regions file {fpath_reg}, skipping')
                return None
        elif args.reg_from == 'shared':
            from glob import glob
            shared_beams = args.reg_shared_beams
            if shared_beams == 'all':
                fpattern = '/*-19rfi.hdf5.reg'
            else:
                fpattern = f'/*{args.fpath[-10:]}.reg'

            fpath_regs = glob(os.path.dirname(args.fpath) + fpattern)
            length = len(fpath_regs)
            if length == 1:
                fpath_reg = fpath_regs[0]
            elif length == 0:
                print('can not find any specified regions files, skipping')
                return None
            else:
                raise FileError('Too many regions file! Beams should share only one.')

            if shared_beams != 'all': # string like '1,2,3'
                beams = [int(i) for i in shared_beams.split(',')]
                if nB not in beams:
                    return None

        else:
            project = get_project(args.fpath)
            date = get_date_from_path(args.fpath)
            fpath_reg = sub_patten(args.reg_from, date=date, nB=f'{nB:02d}', project=project)
            if not os.path.exists(fpath_reg):
                print(f'can not find the specified regions file {fpath_reg}, skipping')
                return None


        print(f'read regions from {fpath_reg}')
        regions = read_regions(fpath_reg)
        if regions is None:
            return None
        is_rfi = np.full(self.s2p.shape[:2], False)
        replace_region(is_rfi, regions, fill_value=True)
        return is_rfi

    def get_tr(self,):
        """
        Time domain continuous RFI
        """
        args = self.args
        from .core.rfi_t import mask_rfi_t
        fit_kwargs = {}
        keys = ['tr_n_continue',
                'tr_s_method_freq',
                'tr_s_method_t',
                'tr_s_sigma_freq',
                'tr_s_sigma_t',
                'tr_times',
                'tr_times_s',
                'ext_add', 'ext_frac']
        for key in keys[:-2]:
            fit_kwargs[key[3:]] = getattr(args, key)
        for key in keys[-2:]:
            fit_kwargs[key] = getattr(args, key)

        return mask_rfi_t(self.freq, self.s2p, method='smooth', **fit_kwargs)

    def get_pr(self,):
        """
        Polarized RFI
        """
        args = self.args
        from .core.rfi_polar import mask_rfi_p
        fit_kwargs = {}
        keys = ['pr_s_sigma',
                'pr_times',
                'pr_times_s',
                'ext_add', 'ext_frac']
        for key in keys[:-2]:
            fit_kwargs[key[3:]] = getattr(args, key)
        for key in keys[-2:]:
            fit_kwargs[key] = getattr(args, key)
        p_rfi = mask_rfi_p(self.s2p_mask, **fit_kwargs)
        p_rfi[:,self.protect_use] = False
        return p_rfi

    def get_nr(self):
        """
        Narrowband single channel RFI
        """
        args = self.args
        T = np.nanmean(self.s2p_mask, axis = 2)
        from .ripple.mark_timeRFI import mask_freq_rfi

        narr_args = {}
        keys = ['nr_mean_times',
                'nr_diff_times',
                'nr_mask_rms_times',
                'nr_rfi_width_lim']
        for key in keys:
            narr_args[key[3:]] = getattr(args, key)
        narr_args['frange'] = None
        narr_args['rms_frange'] = args.rms_frange
        narr_args['ext_add'] = 1
        narr_args['thr_type'] = args.lsn_thr_type

        return mask_freq_rfi(T, self.freq, rtype = 'long-time',plot = False,
                              **narr_args)

    def get_lf(self,T):
        """
        Time domain uncontinuous RFI: long freq time RFI
        """
        args = self.args
        from .ripple.mark_timeRFI import mask_time_rfi

        longf_args = {}
        keys = ['lf_frange',
                'lf_mean_times',
                'lf_diff_times',
                'lf_ext_add',
                'lf_mask_rms_times',]
        for key in keys:
            longf_args[key[3:]] = getattr(args, key)
        longf_args['rfi_width_lim'] = args.lf_rfi_last
        longf_args['rms_frange'] = args.rms_frange
        longf_args['thr_type'] = args.lsn_thr_type

        return mask_time_rfi(T, self.freq, rtype = 'long-freq',plot = False,
                              **longf_args)

    def get_sf(self,T,is_rfi = None):
        """
        Time domain uncontinuous RFI: short freq time RFI
        """
        args = self.args
        from .ripple.mark_timeRFI import mask_time_rfi

        shortf_args = {}
        keys = ['sf_file',
                'sf_frange',
                'sf_frange_step',
                'sf_mean_times',
                'sf_diff_times',
                'sf_ext_add',
                'sf_mask_rms_times',]
        for key in keys:
            shortf_args[key[3:]] = getattr(args, key)
        shortf_args['rfi_width_lim'] = args.sf_rfi_last
        shortf_args['rms_frange'] = args.rms_frange
        shortf_args['thr_type'] = args.lsn_thr_type

        if args.sf_use_time_only and args.all_beams and is_rfi is not None:
            whole_rfi = np.all(is_rfi, axis = 1)
            is_rfi[whole_rfi] = False
            is_timerfi = np.any(is_rfi, axis = 1)
            shortf_args['is_timerfi'] = is_timerfi

        return mask_time_rfi(T, self.freq, rtype = 'short-freq',plot = False,
                               **shortf_args)

    def get_time_rfi(self,is_rfi = None, is_rfi_tmp = None):
        """
        lf & sf
        """
        args = self.args
        T = self.s2p_mean

        t_rfi = np.isnan(T)
        Tt = deepcopy(T)
        Tt[:,self.protect_use] = 0
        if is_rfi_tmp is not None: Tt[is_rfi_tmp] = 0

        if args.lf:
            print('finding lf')
            t_rfi |= self.get_lf(Tt)

        Tt[t_rfi] = 0
        if args.sf:
            print('finding sf')
            t_rfi |= self.get_sf(Tt,is_rfi)

        return t_rfi

    def get_pdr(self):
        """
        Period 8 MHZ RFI
        """
        args = self.args
        T = np.nanmean(self.s2p_mask, axis = 2)
        sm_kwargs = {}
        keys = ['s_method_t', 's_sigma_t', 's_method_freq', 's_sigma_freq',]
        for key in keys:
            sm_kwargs[key] = getattr(args, key)
        sm_kwargs['is_rfi'] = np.isnan(T)
        from .ripple.util import do_smooth
        T = do_smooth(T, **sm_kwargs)

        find_args = {}
        keys = ['rfi_width_lim', 'ext_sec', 'freq_thr', 'freq_step', 'rfi_groups',
                'mask_RFI_method', 'freq_from_theory', 'mask_thr', 'ext_edge',
                'small_rfi_times', 'chan_step',
                'mask_all_theory']
        for key in keys:
            find_args[key] = getattr(args, key)

        from .ripple.markRFI import find_RFI, real_rms
        from tqdm import tqdm
        rfi_thr = args.rfi_thr
        pd_rfi = np.full(T.shape, False, dtype=bool)

        print("Looking for period RFI...")
        for tn in tqdm(range(T.shape[0])):
            if tn in self.not_rfi_num:
                spec = deepcopy(T[tn,:])
                RMS = real_rms(spec, self.freq, args.rms_sigma, args.rms_frange)
#                 is_rfi_mw = (spec > RMS * rfi_thr)
                is_rfi = (spec > RMS * rfi_thr) & (~self.protect_use)
                try:
                    pd_rfi[tn,:],_ = find_RFI(spec,self.freq,is_rfi,RMS = RMS,
                                     plot = False,**find_args)
                except ValueError:
                    pd_rfi[tn,:] = True
                    print(f"tn={tn} has a ValueError !")
#                     import traceback
#                     traceback.print_exc()
        print("Finish finding period RFI...")

        time_coherent_per = args.time_coherent_per
        if (time_coherent_per > 0)&(time_coherent_per <1):
            per_use = (np.sum(pd_rfi,axis = 0)/pd_rfi.shape[0] > time_coherent_per)
            pd_rfi[:,per_use] = True

        return pd_rfi

    def get_aoflagger(self):
        args = self.args
        from .core.aoflagger import flag_data
        is_rfi = flag_data(self.s2p, args.af_strategy_file)
        return is_rfi

    def protect_mw(self):
        args = self.args
        freq = self.freq
        mw_frange = args.mw_frange

        if mw_frange is None:
            mw_frange = [1419, 1422]
            if (freq[0] > mw_frange[1]) or (freq[-1] < mw_frange[0]):
                # Do not need to protect mw.
                pass
            else:
                print(f"mw_frange redicting to {mw_frange}, you should check it again.")

        protect_use = (freq>mw_frange[0])&(freq<mw_frange[1])
        return protect_use

    def check_rms_range(self):
        args = self.args
        rms_frange = args.rms_frange
        if (rms_frange is None) or (rms_frange[1] - rms_frange[0] <= 0):
            from .ripple.markRFI import get_rms_frange
            spec = np.nanmean(np.nanmean(self.s2p, axis = 0),axis = -1)
            self.args.rms_frange = get_rms_frange(spec,self.freq,rms_step=5,)
        else:
            freq = self.freq
            if (rms_frange[0] < freq[0]) or (rms_frange[1] > freq[-1]):
                raise ValueError(f"rms frange {rms_frange} is not in freq range {[freq[0],freq[-1]]}.")

    def _guess_19rfi(self,):
        args = self.args

        from glob import glob
        import re
        subname = args.fpath
        key = re.findall(r'-M[0-1][0-9]', os.path.basename(subname))[0]
        namelist = os.path.basename(subname).split(key)
        name = os.path.join(args.outdir, namelist[0]+'-M*'+namelist[1][:-5]+'-19rfi.hdf5')
#         print(name)
        _19names = glob(name)
        return _19names

    def load_19rfi(self):
        import h5py
        
        f = h5py.File(self._19name,'r')
        fs = f['S']
        freq = fs['freq'][()]
        is_rfi = fs['is_rfi'][()]
        if 'T' in fs.keys():
            infield = 'T'
        elif 'Ta' in fs.keys():
            infield = 'Ta'
        elif 'flux' in fs.keys():
            infield = 'flux'
        else:
            raise(ValueError('can not find spec'))
        s2p = fs[infield][()]
        if s2p.shape[0] == 2 or s2p.shape[0] == 1:
            from .utils.io import PolarMjdChan_to_MjdChanPolar
            s2p = PolarMjdChan_to_MjdChanPolar(s2p)

        f.close()
        
        num_19rfi = s2p[:, :, 0]
        is_19rfi = (num_19rfi >= self.args.b19_min_percent)
        
        print(f"rfi number in 19 beams rfi file: {np.sum(is_19rfi)}")
            
        # protect extended signals in low z
        fuse = (freq > 1410) & (freq < 1421)
        is_19rfi[:, fuse] = False
        self.is_19rfi = is_19rfi

        return is_rfi

    @staticmethod
    def gen_19rfi_file(paras=[]):
        """
        run hifast.rfi_multi
        """
        import subprocess
        command = [sys.executable, '-m', 'hifast.rfi_multi'] + paras
        print('run:')
        print(' '.join(command))
        sys.stdout.flush()
        process = subprocess.Popen(
            command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # process.wait()
        print(*process.communicate())
        if process.returncode != 0:
            raise(ValueError(f'fail to generate the 19 beams rfi file!'))

    def gen_s2p_out(self,):
        args = self.args

        if args.all_beams:
            _19names = self._guess_19rfi()
            if len(_19names) == 0:
                print(f"Try to generate 19 beams rfi file ...")
                argv = ' '.join(sys.argv)
                paras = argv.split()
                paras.pop(0) # delete the first one '... hifast.rfi '
                for s in paras:
                    if s in ['tr', 'pr', 'nr', 'lf', 'pdr']: # except sf
                        s = False

                self.gen_19rfi_file(paras)
                sys.stdout.flush()
                self._19name = os.path.join(args.outdir, os.path.basename(args.fpath).split('.hdf5')[0] + '-19rfi.hdf5')

            elif len(_19names) == 1:
                self._19name = _19names[0]
                print(f"{self._19name} exists. Use it!")
            else:
                raise ValueError("Master Skywalker, there are too many of them. What should we do?")

        self.s2p_out = self.s2p[:]

    def gen_is_rfi(self):
        args = self.args
        self.gen_s2p_out()
        self.s2p = self.s2p[:]
        self.s2p_mean = np.mean(self.s2p,axis = 2)

        # average in 19 beams
        if args.all_beams:
            if os.path.exists(self._19name):
                is_trfi = self.load_19rfi()
            else:
                raise FileNotFoundError(f"Can't find {self._19name}. Did it run hifast.rfi_multi?")

        is_rfi = np.isnan(self.s2p_mean) | np.isinf(self.s2p_mean)
        self.s2p[is_rfi] = np.nan
        self.s2p_mean[is_rfi] = np.nan

        # manual regions
        is_rfi_tmp = self.get_from_regions()
        if is_rfi_tmp is not None:
            is_rfi |= is_rfi_tmp

        self.protect_use = self.protect_mw()
        self.check_rms_range()

        # 19 beams RFI
        if args.all_beams & args.b19:
            is_rfi |= self.is_19rfi

        # long or short freq time RFI
        if args.lf or args.sf:
            if args.all_beams and args.sf_use_time_only:
                    is_rfi |= self.get_time_rfi(is_trfi, is_rfi_tmp)
                    print("Use is_timerfi(sf) in 19 beams rfi :P")
            else:
                if args.all_beams:
                    is_rfi |= is_trfi
                    args.sf = False
                    print("Use is_rfi(sf) in 19 beams rfi :P")
                is_rfi |= self.get_time_rfi(is_rfi = None, is_rfi_tmp = is_rfi_tmp)

        whole_rfi = np.all(is_rfi,axis = 1)
        self.not_rfi_num = np.arange(is_rfi.shape[0])[~whole_rfi]
        self.is_rfi_num = np.arange(is_rfi.shape[0])[whole_rfi]

        # time continuous RFI
        if args.tr:
            print('finding tr')
            is_rfi |= self.get_tr()

        if is_rfi.any():
            self.s2p_mask = deepcopy(self.s2p)
            self.s2p_mask[is_rfi,:] = np.nan
        else:
            self.s2p_mask = self.s2p

        # narrowband RFI
        if args.nr:
            print('finding nr')
            is_rfi |= self.get_nr()

        # 8.1 MHz period RFI
        if args.pdr:
            print('finding period rfi')
            is_rfi |= self.get_pdr()

        # polarized RFI
        if args.pr:
            print('finding pr')
            is_rfi |= self.get_pr()

        if args.aoflagger:
            is_rfi |= self.get_aoflagger()

        # existing RFI
        if 'is_rfi' in self.fs.keys():
            is_rfi |= self.fs['is_rfi'][:]

        return is_rfi


    def __call__(self, save=True):
        args = self.args
        is_rfi = self.gen_is_rfi()
        # change s2p_out if
        if args.nobld:
            self.s2p_out = self._load_s2p_ori()[:]
        if args.replace_rfi:
            self.s2p_out[is_rfi] = np.nan
        self.gen_dict_out(is_rfi = is_rfi)
        # replace outfield as h5py.ExternalLink
        if self.dict_in is None and (not args.replace_rfi) and (not args.nobld):
            self.dict_out[self.outfield] = h5py.ExternalLink(os.path.relpath(
                args.fpath, os.path.dirname(self.fpath_out)), f'/S/{self.infield}')
        if save:
            self.save()


if __name__ == '__main__':
    # for clarity and testing purpose
    dests_hide = [
                 'rms_sigma',
                 's_method_t',
                 's_sigma_t',
                 's_method_freq',
                 's_sigma_freq',
                 'pdr',
                 'rfi_thr',
                 'rfi_width_lim',
                 'ext_sec',
                 'freq_thr',
                 'freq_step',
                 'rfi_groups',
                 'mask_RFI_method',
                 'freq_from_theory',
                 'mask_thr',
                 'ext_edge',
                 'small_rfi_times',
                 'chan_step',
                 'mask_all_theory',
                 'time_coherent_per'
                 ]
    hide_paras(parser, dests_hide)

    args_ = parser.parse_args()
    # print(parser.format_help())
    # print("----------")
    print('#'*35+'Args'+'#'*35)
    args_from = parser.format_values()
    args_from = del_paras_in_string(args_from, dests_hide)
    print(args_from)
    print('#'*35+'####'+'#'*35)

    HistoryAdd = {'args_from': args_from} if args_.my_config is not None else None
    io = IO(args_, HistoryAdd=HistoryAdd)
    io()
