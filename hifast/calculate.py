__all__ = ['IO']

# Cell
from .utils.io import *
#nbdev_comment _all_ = ['parser']

# Internal Cell
sep_line = '##'+'#'*70+'##'
parser = ArgumentParser(prog=f"python -m hifast.{os.path.basename(sys.argv[0])[:-3]}", formatter_class=formatter_class, allow_abbrev=False,
                        description='one array minus another', )
add_common_argument(parser)
parser.add_argument('fpath',
                    help='input spectra file path (x).')
parser.add_argument('--frange', type=float, nargs=2, default=[0, float('inf')],
                    help='Limit frequence range')
parser.add_argument('--y_fpattern',
                    help='second input path(y)')
parser.add_argument('--fpart',
                    help='append fpart')

parser.add_argument('--calculate', default='+', choices=['+','*','/'],
                    help='calculate type')
parser.add_argument('--A', type=float, default=1,
                    help='(Ax) o (By), where o is +,* or /')
parser.add_argument('--B', type=float, default=-1,
                    help='(Ax) o (By), where o is +,* or /')

# Cell
class IO(BaseIO):
    ver = 'old'
    def _get_fpart(self,):
        """
        need modify this function
        """
        args = self.args
        fpart = args.fpart
        if fpart is None:
            fpart = '-'
        else:
            fpart = '-' + fpart
        return fpart
    
    def _import_m(self,):
        """
        need modify this function
        """
        super()._import_m()
        global h5py, np
        import h5py
        import numpy as np
        
    def _gen_s2p_y_fpath(self,):
        """
        find file path of y
        """
        args = self.args
        if args.y_fpattern is not None:
            print('using the path from --y_fpattern as the spectra file of y')
            y_fpath = self.replace_nB(args.y_fpattern, self.nB)
            if not os.path.exists(y_fpath):
                raise(FileNotFoundError(f'file {y_fpath} not exists'))
            self.s2p_y_fpath = y_fpath
        else:
            raise ValueError("You should input --y_fpattern!")
        
    def _load_s2p_y(self):
        
        self._gen_s2p_y_fpath()
        print(f'load y from: \n{self.s2p_y_fpath}')
        f = h5py.File(self.s2p_y_fpath, 'r')['S']
        freq = f['freq'][:]
        s2p_y = f[self.outfield][:]
        is_use_freq = (freq >= self.freq.min()) & (freq <= self.freq.max())
        if not np.all(is_use_freq):
            inds = np.where(is_use_freq)[0]
            s2p_y = s2p_y[..., inds[0]:inds[-1]+1]  # freq axis at end
        s2p_y = PolarMjdChan_to_MjdChanPolar(s2p_y)
        return s2p_y
    
    def gen_s2p_out(self,):
        args = self.args
        # gen self.s2p_out
        s2p = self.s2p[:]
        s2p_y = self._load_s2p_y()
        
        if s2p.shape != s2p_y.shape:
            raise ValueError(f"x's shape is {s2p.shape}, but y's shape is {s2p_y.shape}!")
        
        calculate = args.calculate
        A = args.A
        B = args.B
        print(f"process:({A} * {s2p.shape}) {calculate} ({B} * {s2p_y.shape})")
        if calculate in ['+','*','/']:
            exec(f"self.s2p_out = (A * s2p) {calculate} (B * s2p_y)")
        else:
            raise ValueError(f"Unsupport calculation type {calculate}.")
        

# Cell
if __name__ == '__main__':
    args_ = parser.parse_args()
    # print(parser.format_help())
    # print("----------")
    print('#'*35+'Args'+'#'*35)
    print(parser.format_values())  # useful for logging where different settings came from
    print('#'*35+'####'+'#'*35)
    io = IO(args_)
    io()
