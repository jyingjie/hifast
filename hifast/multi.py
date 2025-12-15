
"""
hifast.multi
============

Multi-step operations for spectral data processing.

This module performs several key processing steps on the input spectra:
1.  **Reference Frame Correction**: Converts velocities from the telescope's Topocentric frame to a standard frame (e.g., LSRK).
2.  **RFI Masking**: Replaces flagged RFI data with NaN.
3.  **Polarization Merging**: Averages XX and YY polarizations.
4.  **Extrapolation Masking**: Optionally masks data points with extrapolated coordinates.
"""

__all__ = ['IO', 'parser']


from .utils.io import *


def create_parser():
    sep_line = '##'+'#'*70+'##'
    parser = ArgumentParser(prog=f"python -m hifast.{os.path.basename(sys.argv[0])[:-3]}", formatter_class=formatter_class, allow_abbrev=False,
                            description='Perform multi-step operations: RFI masking, reference frame correction, and polarization merging.')
    add_common_argument(parser)

    # --- Input/Output ---
    group = parser.add_argument_group('Input/Output')
    group.add_argument('fpath', metavar='FILE',
                        help='Input spectra file path (HDF5 format).')
    group.add_argument('--frange', type=float, nargs=2, default=[0, float('inf')], metavar=('MIN', 'MAX'),
                        help='Limit the frequency range to process [MHz]. Followed by two values (min, max).')

    # --- Frame Correction ---
    group = parser.add_argument_group('Frame Correction')
    group.add_argument('--fc', type=bool_fun, choices=[True, False], default='True', not_in_write_out_config_file=True,
                       help='Enable reference frame correction (Topocentric -> LSRK/Heliocentric).')
    group.add_argument('--frame', choices=['BARYCENT', 'HELIOCEN', 'LSRK', 'LSRD'], default='LSRK',
                       help='Target reference frame. "LSRK" (Local Standard of Rest, Kinematic) is commonly used.')
    group.add_argument('--vtype', choices=['radio', 'optical'], default='optical',
                       help='Velocity definition type. "optical" (cz) or "radio". Note: This is primarily for tracking/single-spectrum analysis. `hifast.cube` calculates velocity from frequency and ignores this column.')

    # --- Data Processing ---
    group = parser.add_argument_group('Data Processing')
    group.add_argument('--replace_rfi', type=bool_fun, choices=[True, False], default='True',
                       help='Replace RFI-contaminated spectra with NaN (requires "is_rfi" dataset in input).')
    group.add_argument('--mask_extrapo', type=bool_fun, choices=[True, False], default='False',
                       help='Mask data points where RA/DEC coordinates were extrapolated (requires "is_extrapo" dataset).')
    group.add_argument('--merge_polar', type=bool_fun, choices=[True, False], default='True',
                       help='Merge (average) the two polarizations (XX and YY) into a single Stokes I intensity.')
    
    return parser

parser = create_parser()


class IO(BaseIO):
    ver = 'old'
    def _get_fpart(self,):
        """
        need modify this function
        """
        fpart = '-fc'
        return fpart

    def __call__(self, save=True):
        import numpy as np
        args = self.args
        s2p = self.s2p[:]
        # try load is_rfi
        if 'is_rfi' in self.fs.keys():
                # if not replace_rfi, frame correct it

                if self.is_use_freq is not None:
                    inds = np.where(self.is_use_freq)[0]
                    is_rfi = self.fs['is_rfi'][:, inds]
                else:
                    is_rfi = self.fs['is_rfi'][:]
                if args.replace_rfi:
                    print('replacing RFI')
                    s2p[is_rfi] = np.nan
        else:
            is_rfi = None

        # mask data points where ra/dec are extrapolated
        if args.mask_extrapo and 'is_extrapo' in self.fs.keys():
            is_extrapo = self.fs['is_extrapo'][:]
            print('masking data points with extrapolated ra/dec')
            # is_extrapo shape: (Mjd,), s2p shape: (Mjd, Chan, Polar)
            s2p[is_extrapo, :, :] = np.nan

        if args.merge_polar and s2p.ndim == 3:
            print('average two polarization...')
            s2p = np.mean(s2p, axis=2, keepdims=True)
        # frame
        if args.fc:
            from .core.corr_vel import freq2vel, correct_spec
            if 'frame' in self.Header.keys():
                raise(ValueError(f'rest frame already corrected... please set --fc as False'))
                sys.exit()
            else:
                self.Header['frame'] = args.frame
                self.Header['vel_type'] = args.vtype
            # correct is_rfi
            if is_rfi is not None and not args.replace_rfi:
                is_rfi, _ = correct_spec(is_rfi, self.freq, self.ra, self.dec, self.mjd,
                                           frame=args.frame, method='interp', interp_kind='nearest')
                is_rfi = np.array(is_rfi, dtype=bool)
            print('frame correcting...')
            s2p, freq = correct_spec(s2p, self.freq, self.ra, self.dec, self.mjd, frame=args.frame, method='interp')
            vel = freq2vel(freq, vtype=args.vtype)
            self.s2p_out = s2p
            self.gen_dict_out(freq=freq, vel=vel)  # replace freq, add vel
            if 'Tcal' in self.dict_out.keys():
                # use ``try`` to be compatible with the file generated with old version having bug
                try:
                    values = self.dict_out['Tcal'] # self.is_use_freq has been applied; still is (Mjd,Chan,Polar)
                    self.dict_out['Tcal'], _ = correct_spec(values, self.freq, self.ra, self.dec, self.mjd, frame=args.frame, method='interp')
                except:
                    pass
        else:
            self.s2p_out = s2p
            self.gen_dict_out()
        # save is_rfi only if replace_rfi is False
        if not args.replace_rfi and is_rfi is not None:
            self.dict_out['is_rfi'] = is_rfi
        else:
            if 'is_rfi' in self.dict_out.keys():
                self.dict_out.pop('is_rfi')
        if save:
            self.save()


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
