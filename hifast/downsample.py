__all__ = ['IO']

from .utils.io import *

sep_line = '##'+'#'*70+'##'
parser = ArgumentParser(prog=f"python -m hifast.{os.path.basename(sys.argv[0])[:-3]}", formatter_class=formatter_class, allow_abbrev=False,
                        description='downsample.', )
add_common_argument(parser)
parser.add_argument('fpath',
                    help='input spectra file path.')
parser.add_argument('--frange', type=float, nargs=2, default=[0, float('inf')],
                    help='Limit frequence range')
group = parser.add_argument_group(f'*downsample\n{sep_line}')
group.add_argument('--chan_factor',type=int, default = 0, 
                    help='down sample on channel')
group.add_argument('--spec_factor',type=int, default = 0, 
                    help='down sample on spec')

if __name__ == '__main__':
    args_ = parser.parse_args()
    print('#'*35+'Args'+'#'*35)
    print(parser.format_values())
    print('#'*35+'####'+'#'*35)

import numpy as np
from astropy import log
from .utils.misc import down_sample

def down_sample_3d(data, dfactor, axis = -1):
    if len(data.shape) == 1:
        ret = down_sample(data[None,:,None],dfactor)[0,:,0]
    elif len(data.shape) == 2:
        if axis == 0:
            trans = (1,0)
        elif axis == 1 or axis == -1:
            trans = (0,1)
        ret = down_sample(data.transpose(trans)[:,:,None],dfactor)[:,:,0].transpose(trans)
    elif len(data.shape) == 3:
        if axis == 0:
            trans1 = (2,0,1); trans_ = (1,2,0)
        elif axis == 1 or axis == -2:
            trans = (0,1,2); trans_ = trans
        elif axis == 2 or axis == -1:
            trans = (1,2,0); trans_ = (2,0,1)
        ret = down_sample(data.transpose(trans),dfactor).transpose(trans_)
    return ret

class IO(BaseIO):
    ver = 'old'
    def _get_fpart(self,):
        """
        need modify this function
        """
        fpart = '-ds'
        return fpart
    
    def gen_dict_out(self,chan_factor=0,spec_factor=0):
        fs = self.fs
        dict_out = {}
        if spec_factor > 0:
            length = len(self.mjd)
            for key in fs.keys():
                x = fs[key][:]
                if length in x.shape:
                    n = np.where(np.array(x.shape)==length)[0][0]
                    x1 = down_sample_3d(x, spec_factor, axis = n)
                    dict_out[key] = x1 
                    
            print("time axis: before:", x.shape, "after:", x1.shape)
            
        if chan_factor > 0: 
            length = len(self.freq)
            for key in fs.keys():
                if key in dict_out.keys():
                    x = dict_out[key][:]
                else:
                    x = fs[key][:]
                if length in x.shape:
                    n = np.where(np.array(x.shape)==length)[0][0]
                    x1 = down_sample_3d(x, chan_factor, axis = n)
                    dict_out[key] = x1 
            print("freq axis: before:", x.shape, "after:", x1.shape) 
            
        for key in fs.keys():
            if key not in dict_out.keys():
                dict_out[key] = fs[key][:]

        dict_out['Header'] = self.Header
        self.dict_out = dict_out

    def __call__(self,):
        args = self.args
        
        self.gen_dict_out(args.chan_factor, args.spec_factor)
        self.save()

if __name__ == '__main__':
    io = IO(args_)
    io()
