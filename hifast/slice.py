

__all__ = ['IO']


from .utils.io import *
import copy


sep_line = '##'+'#'*70+'##'
parser = ArgumentParser(prog=f"python -m hifast.{os.path.basename(sys.argv[0])[:-3]}", formatter_class=formatter_class, allow_abbrev=False,
                        description='', )
add_common_argument(parser)
parser.add_argument('fpath',
                    help='input spectra temperature file path.')
parser.add_argument('--frange', type=float, nargs=2, default=[0, float('inf')],
                    help='Limit frequence range')
parser.add_argument('--irange', type=int, nargs=2,
                    help="")


class IO(BaseIO):
    ver = 'new'

    def _get_fpart(self,):
        """
        need modify this function
        """
        fpart = '-sli'
        return fpart

    def __call__(self, save=True):
        import numpy as np

        args = self.args
        type1 = ['ra', 'dec',  'is_delay', 'is_extrapo', 'is_on', 'mjd', 'next_to_cal']
        type2 = ['freq', 'vel']

        dict_out = {}
        dict_out["Header"] = self.Header
        for key in self.fs.keys():
            if key in type1:
                dict_out[key] = self.fs[key][args.irange[0]:args.irange[1]]
            elif key in type2:
                if self.is_use_freq is not None:
                    inds = np.where(self.is_use_freq)[0]
                    dict_out[key] = self.fs[key][inds[0]:inds[-1]+1]
                else:
                    dict_out[key] = self.fs[key]
            elif key == 'Tcal':
                if self.is_use_freq is not None:
                    inds = np.where(self.is_use_freq)[0]
                    dict_out[key] = self.fs[key][:, inds[0]:inds[-1]+1]
                else:
                    dict_out[key] = self.fs[key]
            elif key in ['is_excluded', 'is_rfi']:
                if self.is_use_freq is not None:
                    inds = np.where(self.is_use_freq)[0]
                    dict_out[key] = self.fs[key][args.irange[0]:args.irange[1], inds[0]:inds[-1]+1]
                else:
                    dict_out[key] = self.fs[key][args.irange[0]:args.irange[1]]
            elif key == self.infield:
                dict_out[key] = self.s2p[:, args.irange[0]:args.irange[1]]

        self.dict_out = dict_out
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
    args_.no_radec = True
    io = IO(args_, HistoryAdd=HistoryAdd)
    io()
