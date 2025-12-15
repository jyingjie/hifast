
"""
hifast.flux
===========

Flux calibration module.

This module converts the unit of spectra from Antenna Temperature (Ta, Kelvin) to Flux Density (Jy).
It supports two calibration methods:
1.  **Noise Diode Calibration**: Uses a calibrator source (e.g., 3C48) observed with a noise diode to determine the gain.
2.  **Pre-measured Gain**: Uses a standard gain curve (e.g., Jiang et al. 2020) based on zenith angle.
"""

__all__ = ['IO', 'parser']


from .utils.io import *


def create_parser():
    sep_line = '##'+'#'*70+'##'
    parser = ArgumentParser(prog=f"python -m hifast.{os.path.basename(sys.argv[0])[:-3]}", formatter_class=formatter_class, allow_abbrev=False,
                            description='Convert spectra from Antenna Temperature (Ta) to Flux Density (Jy).')
    add_common_argument(parser)

    # --- Input/Output ---
    group = parser.add_argument_group('Input/Output')
    group.add_argument('fpath', metavar='FILE',
                        help='Input spectra file path (HDF5 format, typically with units in Ta).')

    # --- Calibration Method: Noise Diode ---
    group = parser.add_argument_group('Calibration: Noise Diode')
    group.add_argument('--cbr_store', '--cali_fname', dest='cbr_store', default='none', metavar='FILE/DIR',
                       help='Path to the calibrator file or directory containing calibrator files. If specified, noise diode calibration is used.')
    group.add_argument('--cbr_name', default='*', metavar='NAME',
                       help='Name pattern of the calibrator source (e.g., "3C48", "3C286"). Used to filter files if a directory is provided.')
    group.add_argument('--only_use_19beams', type=bool_fun, choices=[True, False], default='False',
                       help='If True, only use calibrator data that has all 19 beams.')
    group.add_argument('--fix_diff_tcal', type=bool_fun, choices=[True, False], default='True',
                       help='Correct for Tcal differences between the target observation and the calibrator observation.')
    group.add_argument('--fix_diff_ZA', type=bool_fun, choices=[True, False], default='False',
                       help='Correct for gain differences due to Zenith Angle (ZA) differences between target and calibrator.')

    # --- Calibration Method: Pre-measured Gain ---
    group = parser.add_argument_group('Calibration: Pre-measured Gain')
    group.add_argument('--pre_measured', default='Jiang2020', choices=['Jiang2020', 'Liu2024'],
                       help='Select the pre-measured gain curve to use if `--cbr_store` is not set. "Jiang2020" (arXiv:2002.01786) is the standard.')
    group.add_argument('--Atemp', '--Tamb', type=float, metavar='DEG_C',
                       help='Ambient temperature [Celsius]. Required for "Liu2024" model correction. Ignored for "Jiang2020".')

    # --- General Settings ---
    group = parser.add_argument_group('General Settings')
    group.add_argument('--frange', type=float, nargs=2, default=[0, float('inf')], metavar=('MIN', 'MAX'),
                        help='Limit the frequency range to process [MHz]. Followed by two values (min, max).')
    group.add_argument('--no_radec', action='store_true', not_in_write_out_config_file=True,
                        help="Skip checking or adding RA/DEC coordinates (use with caution).")
    
    # Internal/Hidden
    parser.add_argument('--flux', type=bool_fun, choices=[True], default='True', not_in_write_out_config_file=True,
                       help=argparse.SUPPRESS)

    return parser

parser = create_parser()


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
