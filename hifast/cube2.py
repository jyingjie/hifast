

__all__ = []


from .utils.io import *
from .cube import parser
import copy


parser.add_argument('--share_mem_dir', default='/dev/shm',
                   help='Default is `/dev/shm` which is in memory.')
hide_paras(parser, 'share_mem')
for act in parser._actions:
    if act.dest == 'method':
        act.choices += ['mean', 'median', 'wmean', 'wmedian']
        break


if __name__ == '__main__':
    import os
    import tempfile
    from .core.image2 import Imaging2
    args_ = parser.parse_args()
    args_.share_mem = True
    tempfile.tempdir = args_.share_mem_dir

    img2 = Imaging2(args_)
    img2()
