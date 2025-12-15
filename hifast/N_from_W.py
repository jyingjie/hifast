

__all__ = ['IO', 'interp_fun', 'align_calc']


from .utils.io import *


sep_line = '##'+'#'*70+'##'
parser = ArgumentParser(prog=f"python -m hifast.{os.path.basename(sys.argv[0])[:-3]}", formatter_class=formatter_class, allow_abbrev=False,
                        description='', )
add_common_argument(parser)
parser.add_argument('fpath',
                    help='input spectra temperature file path.')
parser.add_argument('--frange', type=float, nargs=2, default=[0, float('inf')],
                    help='Limit frequence range')
parser.add_argument('--no_radec', action='store_true', not_in_write_out_config_file=True,
                    help="if set, don't check or add ra dec")
group = parser.add_argument_group(f'*Flux\n{sep_line}')
# --flux affect the outfield name
group.add_argument('--final_fpath',
                   help='')
group.add_argument('--flux_fpath',
                   help='')
group.add_argument('--flux_fpath_bef',
                   help='')
group.add_argument('--T_fpath',
                   help='')


class IO(BaseIO):
    ver = 'old'

    def _get_fpart(self,):
        """
        need modify this function
        """
        fpart = '-fW'
        return fpart

    def _import_m(self,):
        """
        need modify this function
        """
        super()._import_m()
        global HFSpec, np, interp1d
        from .funcs import HFSpec
        import numpy as np
        from scipy.interpolate import interp1d

    def _process_fpaths(self,):
        args = self.args
        project = get_project(args.fpath)
        date = get_date_from_path(args.fpath)
        nB = get_nB(args.fpath)

        for pname in ['flux_fpath',
                     'flux_fpath_bef',
                     'final_fpath',
                     'T_fpath']:
            path = getattr(args, pname)
            path = sub_patten(path, date=date, nB=f'{nB:02d}', project=project)
            path = os.path.expanduser(path)
            setattr(args, pname, path)

    def gen_s2p_out(self,):
        args = self.args

        self._process_fpaths()
        S_flux = HFSpec(args.flux_fpath)
        S_flux_bef = HFSpec(args.flux_fpath_bef)
        S_final = HFSpec(args.final_fpath)
        S_T = HFSpec(args.T_fpath)

        for S in [S_flux, S_flux_bef, S_final, S_T]:
            assert (S.mjd == self.mjd).all()
        # get gain
        G_freq, Jy_K = align_calc(S_flux.freq, S_flux.flux, S_flux_bef.freq, S_flux_bef.Ta, axis=1, fun=lambda a,b:a/b)
        # get baseline in flux
        freq_tmp, flux_tmp = align_calc(G_freq, Jy_K, S_T.freq, S_T.Ta, axis=1, fun=lambda a,b:a*b)
        bl_freq, bl_flux = align_calc(freq_tmp, flux_tmp, S_final.freq, S_final.flux, axis=1, fun=lambda a,b:a-b)
        del freq_tmp, flux_tmp

        # Temperature cali; only process 'merge' and 'div' for now
        print('Temperature Calibration')
        print(f'Using info from {args.T_fpath}')
        s2p_out = self.s2p / (interp_fun(S_T.freq, S_T.pcals_merged_s, axis=1)(self.freq) * S_T.pcals_amp_diff_interp_values[:][:, None, :])
        # power of W is 16 times of N ??
        s2p_out *= 16
        s2p_out[self.fs['is_on']] -= 1
        s2p_out *= interp_fun(S_T.freq, S_T.Tcal, axis=1)(self.freq)
        # flux cali
        print('Flux Calibration')
        print(f'Using info from {args.flux_fpath} \nand\n{args.flux_fpath_bef} ')
        s2p_out *= interp_fun(G_freq, Jy_K, axis=1)(self.freq)
        self.outfield = 'flux'
        del G_freq, Jy_K
        # Sub baseline
        print('Baseline and SW removal')
        s2p_out -= interp_fun(bl_freq, bl_flux, axis=1)(self.freq)
        self.s2p_out = s2p_out
        print('RFI')
        if 'is_rfi' in S_final.keys():
            is_rfi = interp1d(S_final.freq, S_final.is_rfi, axis=1, kind='nearest', bounds_error=False, fill_value=0)(self.freq)
            is_rfi = is_rfi.astype(bool)
        else:
            is_rfi = None
        self.is_rfi = is_rfi

    def __call__(self, save=True):
        args = self.args
        self.gen_s2p_out()
        self.gen_dict_out()

        if self.is_rfi is not None:
            self.dict_out['is_rfi'] = self.is_rfi
        if save:
            self.save()


def interp_fun(*args, **kwrags):
    from scipy.interpolate import interp1d
    return interp1d(*args, bounds_error=False, fill_value=np.nan, **kwrags)

def align_calc(freq1, arr1, freq2, arr2, axis=-1, fun=None):
    """
    freq1:
    """
    freq_range = [max(freq1.min(), freq2.min()), min(freq1.max(), freq2.max())]
    is_1 = (freq1 > freq_range[0]) & (freq1 < freq_range[1])
    is_2 = (freq2 > freq_range[0]) & (freq2 < freq_range[1])

    assert (freq1[is_1] == freq2[is_2]).all()
    freq_new = freq1[is_1]

    sli1 = [slice(None),] * arr1.ndim
    sli1[axis] = np.where(is_1)[0]
    sli2 = [slice(None),] * arr2.ndim
    sli2[axis] = np.where(is_2)[0]
    sli1 = tuple(sli1)
    sli2 = tuple(sli2)

    return freq_new, fun(arr1[sli1], arr2[sli2])


if __name__ == '__main__':
    args_ = parser.parse_args()
    # print(parser.format_help())
    # print("----------")
    print('#'*35+'Args'+'#'*35)
    args_from = parser.format_values()
    print(args_from)
    print('#'*35+'####'+'#'*35)

    HistoryAdd = {'args_from': args_from} if args_.my_config is not None else None
    io = IO(args_, HistoryAdd=HistoryAdd)
    io()
