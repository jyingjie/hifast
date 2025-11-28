

__all__ = ['sep_line', 'parser']


from .utils.io import *
from glob import glob
import json


sep_line = '##'+'#'*70+'##'
parser = ArgumentParser(prog=f"python -m hifast.{os.path.basename(sys.argv[0])[:-3]}",
                        formatter_class=formatter_class, allow_abbrev=False,
                        description='cut FAST RAW Data', )

parser.add_argument('fpath',
                    help='only need input the first chunk file path (e.g. XXX_0001.fits) of RAW spectra data')
parser.add_argument('-f', dest='force', action='store_true',
                    help='overwriting file if out file exists')
parser.add_argument('--outdir', required=True,
                   help='the directory to store output files.')
parser.add_argument('-g', is_write_out_config_file_arg=True,
                        help='save config to file path')
parser.add_argument('-c', '--my-config', is_config_file_arg=True,
                        help='config file path')
parser.add_argument('--frange', type=float, nargs=2,
                    help='freq range')
parser.add_argument('--step', type=int, default=1,
                   help='number of files processed every time, default 1')
parser.add_argument('--start', type=int, default=1,
                   help='chunk number of start')
parser.add_argument('--stop', type=int,
                   help='chunk number of stop')
parser.add_argument('--sep_save', type=bool_fun, choices=[True, False], default='False',
                   help='every step save to a file')
parser.add_argument('--h5_compression', default='none',
                   help='`compression` in `h5py create_dataset`; `none`, `lzf`, `gzip`, or a number in range(10). '
                         'Note: compressed file may not be opened by carta.')


if __name__ == '__main__':
    args_ = parser.parse_args()
    print('#'*35+'Args'+'#'*35)
    print(parser.format_values())
    print('#'*35+'####'+'#'*35)

    args = args_

    # check h5_compression
    try:
        args.h5_compression = int(args.h5_compression)
    except ValueError:
        pass

    #record history
    header = rec_his(args=json.dumps(args.__dict__))

    #check input file exits
    if not os.path.exists(args.fpath):
        raise(OSError(f'File {args.fpath} not exists.'))

    # replace patten in outdir
    nB = get_nB(args.fpath)
    project = get_project(args.fpath)
    ## use the dirname of the fits file as "date"
    date = os.path.basename(os.path.dirname(os.path.abspath(args.fpath)))
    args.outdir = sub_patten(args.outdir, date=date, nB=f'{nB:02d}', project=project)
    print(f'outdir: {args.outdir}')
    if args.outdir is not None:
        if not os.path.exists(args.outdir):
            print(f'outdir {args.outdir} not exists. Create it now')
            os.makedirs(args.outdir, exist_ok=True)
    ## check out file
    fname_add = date
    fname_part = re.sub('[0-9]{4}\.fits\Z', '', args.fpath)
    out_name_base = os.path.join(args.outdir, f"{os.path.basename(fname_part)}")

    fileout = out_name_base + rf"0001.hdf5"
    fileout_sep = glob(out_name_base + rf"[0-9][0-9][0-9][0-9].hdf5")


    if os.path.exists(fileout) or len(fileout_sep)>0:
        if args.force:
            print(f"will overwrite the existing out file")
        else:
            print(f"File exists {fileout}")
            print(fileout_sep)
            print('exit... Using -f to overwrite it.')
            sys.exit()
    # run
    from .core.cal2 import FASTRawCut
    Fd = FASTRawCut(args.fpath,
                    start=args.start,
                    stop=args.stop,
                    frange=args.frange,
                    verbose=True)
    Fd(outdir=args.outdir,
       step=args.step,
       header=header,
       sep_save=args.sep_save,
       h5_compression=args.h5_compression,
      )
