

__all__ = ['IO']


from .utils.io import *
import copy


sep_line = '##'+'#'*70+'##'
parser = ArgumentParser(prog=f"python -m hifast.{os.path.basename(sys.argv[0])[:-3]}", formatter_class=formatter_class, allow_abbrev=False,
                        description='downsample.', )
add_common_argument(parser)
parser.add_argument('fpath',
                    help='input spectra file path.')
parser.add_argument('--frange', type=float, nargs=2, default=[0, float('inf')],
                    help='Limit frequence range')
group = parser.add_argument_group(f'*downsample\n{sep_line}')
group.add_argument('--chan_factor', type=int, default = 0,
                    help='down sample on channel')
group.add_argument('--spec_factor', '--t_factor',type=int, default=0,
                    help='down sample on spec')


class IO(BaseIO):
    ver = 'new'

    def __init__(self, args, dict_in=None, inplace_args=False):
        """
        args: class
              including attributes: fpath, outdir, frange
        dict_in: if None, load data from args.fpath, if set, omit data in args.fpath
        """
        args = args if inplace_args else copy.deepcopy(args)
        args.no_radec = True
        super().__init__(args, dict_in=dict_in, inplace_args=inplace_args)

    def _get_fpart(self,):
        """
        need modify this function
        """
        fpart = '-ds'
        return fpart

    def down(self, arr, n, axis=-1, drop=True):
        from .utils.misc import average_every_n
        if n <= 1:
            return arr
        else:
            return average_every_n(arr, n, axis=axis, drop=drop).astype(arr.dtype)

    def __call__(self, save=True):

        from .core.radec import _tight_ra
        import scipy.interpolate as interp

        args = self.args
        type1 = ['is_delay', 'is_extrapo', 'is_on', 'mjd', 'next_to_cal']
        type2 = ['freq', 'vel']
        ## to do: if vel is under optical

        dict_out = {}
        dict_out["Header"] = self.Header
        for key in self.fs.keys():
            if key in type1:
                dict_out[key] = self.down(self.fs[key][()], args.spec_factor, axis=0)
            elif key in type2:
                if self.is_use_freq is not None:
                    dict_out[key] = self.fs[key][()][self.is_use_freq]
                else:
                    dict_out[key] = self.fs[key][()]
                dict_out[key] = self.down(dict_out[key], args.chan_factor, axis=0)
            elif key == 'Tcal':
                if self.is_use_freq is not None:
                    dict_out[key] = self.fs[key][()][:, self.is_use_freq]
                else:
                    dict_out[key] = self.fs[key][()]
                dict_out[key] = self.down(dict_out[key], args.chan_factor, axis=1)
            elif key == self.infield:
                dict_out[key] = self.down(self.down(self.s2p, args.spec_factor, axis=1), args.chan_factor, axis=2)
            elif key == 'is_rfi':
                value = self.fs[key][()]
                if self.is_use_freq is not None:
                    value = value[:, self.is_use_freq]
                dict_out[key] = self.down(self.down(value, args.spec_factor, axis=0), args.chan_factor, axis=1)
        if 'ra' in self.fs.keys():
            if args.spec_factor > 1:
                ra = _tight_ra(self.fs['ra'][()])
                dict_out['ra'] = interp.interp1d(self.fs['mjd'], ra, kind='linear', fill_value ='extrapolate')(dict_out['mjd'])
                dict_out['ra'][dict_out['ra']<0] += 360
            else:
                dict_out['ra'] = self.fs['ra'][()]
        if 'dec' in self.fs.keys():
            if args.spec_factor > 1:
                dict_out['dec'] = interp.interp1d(self.fs['mjd'], self.fs['dec'], kind='linear', fill_value ='extrapolate')(dict_out['mjd'])
            else:
                dict_out['dec'] = self.fs['dec'][()]

        self.dict_out = dict_out
        if save:
            self.save()


if __name__ == '__main__':
    args_ = parser.parse_args()
    print('#'*35+'Args'+'#'*35)
    print(parser.format_values())
    print('#'*35+'####'+'#'*35)
    io = IO(args_)
    io()
