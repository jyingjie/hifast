

__all__ = ['IO']


from .utils.io import *


sep_line = '##'+'#'*70+'##'
parser = ArgumentParser(prog=f"python -m hifast.{os.path.basename(sys.argv[0])[:-3]}", formatter_class=formatter_class, allow_abbrev=False,
                        description='Convert the spectra temperature to flux', )
add_common_argument(parser)
parser.add_argument('fpath',
                    help='input spectra temperature file path.')
parser.add_argument('--frange', type=float, nargs=2, default=[0, float('inf')],
                    help='Limit frequence range')
parser.add_argument('--no_radec', action='store_true', not_in_write_out_config_file=True,
                    help="if set, don't check or add ra dec")

group = parser.add_argument_group(f'*Flux\nUse calibrator (`--cbr_store`) or pre-measured flux gain (`--pre_measured`) \n{sep_line}')
# --flux affect the outfield name
group.add_argument('--flux', type=bool_fun, choices=[True], default='True', not_in_write_out_config_file=True,
                   help='It must be True.')
group = parser.add_argument_group(f'*Use Calibrator\n{sep_line}')
group.add_argument('--cbr_store', '--cali_fname', dest='cbr_store', default='none',
                   help='calibrator file name or directory stored multi-files')
group.add_argument('--cbr_name', default='*',
                   help='calibrator name, e.g. 3C48')
group.add_argument('--only_use_19beams', type=bool_fun, choices=[True, False], default='False',
                   help='only use 19beams calibrator')
group.add_argument('--fix_diff_tcal', type=bool_fun, choices=[True, False], default='True',
                   help='If True, fix the tcal differece in spec and calibrator')
group.add_argument('--fix_diff_ZA', type=bool_fun, choices=[True, False], default='False',
                   help='If True, the difference in gain between the spec and the calibrator will be corrected based on the difference in zenith angle.')
group = parser.add_argument_group(f'*Use Pre-measured flux gain if `--cbr_store` is not set\n{sep_line}')
group.add_argument('--pre_measured', default='Jiang2020', choices=['Jiang2020', 'Liu2024'],
                   help='Specify pre-measured flux gain. Support \'Jiang2020\'(arXiv:2002.01786), \'Liu2024\'(arXiv:..)')
## Ambient temperature correction for Liu2024
group.add_argument('--Atemp', '--Tamb', type=float,
                   help='Used for `--pre_measured Liu2024`. Specify the ambient temperature in degree Celsius. If not specified, no ambient temperature correction will be applied.')


class IO(BaseIO):
    ver = 'old'

    def _get_fpart(self,):
        """
        need modify this function
        """
        fpart = '-flux'
        return fpart

    def gen_s2p_out(self,):
        args = self.args
        s2p = self.s2p[:]
        from .core.flux import FluxCali
        print('Flux calibrating ...')
        if args.fix_diff_tcal:
            tcal_spec = self.fs['Tcal'][:]
            if self.is_use_freq is not None:
                tcal_spec = tcal_spec[:, self.is_use_freq]
        else:
            tcal_spec = None
#         if args.cbr_store != 'none' and args.cbr_store is not None:
#             print(f'using {args.cbr_store} ...')

        fcali = FluxCali(self.nB, self.freq, cbr_store=args.cbr_store, cbr_name=args.cbr_name, tcal_spec=tcal_spec,  only_use_19beams=args.only_use_19beams,
                         pre_measured=args.pre_measured, ATemp=args.Atemp,
                         mjd=self.mjd, ra=self.ra, dec=self.dec, fix_diff_ZA=args.fix_diff_ZA)
        self.s2p_out = fcali(s2p)

        if args.cbr_store is not None and args.cbr_store != 'none':
            self.Header['Calibrater_fpath'] = fcali.cbr_fpath
        else:
            self.Header['pre_measured_fluxgain'] = args.pre_measured


if __name__ == '__main__':
    dests_hide = ['flux',]
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
