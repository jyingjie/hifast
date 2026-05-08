

__all__ = ['sep_line', 'parser']


from .utils.io import *
from .utils.h5compression import H5CompressionConfig
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
                   help=(
                       'HDF5 compression method for /1/DATA. Choices: `none`, `gzip`, `lzf`, '
                       '`bitshuffle_lz4`, `bitshuffle_zstd`, or a number in range(10) as a '
                       'backward-compatible alias for `gzip` level.\n'
                       'Recommended combinations:\n'
                       '  Compatibility/CARTA: `--h5_compression gzip --h5_compression_level 2 '
                       '--h5_chunk_rows 128 --h5_chunk_chans 512`\n'
                       '  Higher compression: `--h5_compression gzip --h5_compression_level 7 '
                       '--h5_chunk_rows 128 --h5_chunk_chans 512`\n'
                       '  Faster plugin mode: `--h5_compression bitshuffle_lz4 '
                       '--h5_chunk_rows 128 --h5_chunk_chans 1024`\n'
                       '  Smaller plugin mode: `--h5_compression bitshuffle_zstd '
                       '--h5_compression_level 5 --h5_chunk_rows 128 --h5_chunk_chans 1024`\n'
                       'Sample benchmark on recent full-frequency 2-pol outputs:\n'
                       '  gzip level 2, 128x512: about 1.50x compression, write ~62-80 s, sampled read ~0.03-0.04 s\n'
                       '  gzip level 7, 128x512: about 1.51x compression, write ~94-104 s, sampled read ~0.09 s\n'
                       '  bitshuffle_lz4, 128x1024: about 1.45x compression, write ~12-18 s, sampled read ~0.02 s\n'
                       '  bitshuffle_zstd level 5, 128x1024: about 1.47x compression, write ~29 s, sampled read ~0.02 s\n'
                       '  Read timings above come from the benchmark script''s sampled read path rather than full-file sequential reads,\n'
                       '  so use them only as a rough comparison between methods.\n'
                       'Note: plugin compression needs `hdf5plugin`; when opening plugin-compressed files '
                       'in other HDF5 programs, you may need to set `HDF5_PLUGIN_PATH` to the plugin directory. '
                       'Prebuilt plugin binaries can be taken from the `hdf5plugin` package: '
                       'https://pypi.org/project/hdf5plugin/'
                   ))
parser.add_argument('--h5_compression_level', type=int,
                   help='compression level for `gzip` or `bitshuffle_zstd`; default `2` for gzip and `5` for bitshuffle_zstd')
parser.add_argument('--h5_chunk_rows', type=int,
                   help='chunk size along the row axis of output DATA; method-dependent default if omitted')
parser.add_argument('--h5_chunk_chans', type=int,
                   help='chunk size along the channel axis of output DATA; method-dependent default if omitted')


if __name__ == '__main__':
    args_ = parser.parse_args()
    print('#'*35+'Args'+'#'*35)
    print(parser.format_values())
    print('#'*35+'####'+'#'*35)

    args = args_

    # normalize h5 compression arguments
    try:
        legacy_level = int(args.h5_compression)
    except ValueError:
        legacy_level = None
    if legacy_level is not None:
        if not 0 <= legacy_level <= 9:
            raise(ValueError('legacy numeric --h5_compression should be in range(10)'))
        if args.h5_compression_level is not None:
            raise(ValueError('do not set both numeric --h5_compression and --h5_compression_level'))
        args.h5_compression = 'gzip'
        args.h5_compression_level = legacy_level
    if args.h5_chunk_rows is not None and args.h5_chunk_rows <= 0:
        raise(ValueError('--h5_chunk_rows should be a positive integer'))
    if args.h5_chunk_chans is not None and args.h5_chunk_chans <= 0:
        raise(ValueError('--h5_chunk_chans should be a positive integer'))

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
    h5_compression_config = H5CompressionConfig(
        compression=args.h5_compression,
        compression_level=args.h5_compression_level,
        chunk_rows=args.h5_chunk_rows,
        chunk_chans=args.h5_chunk_chans,
    )
    Fd(outdir=args.outdir,
       step=args.step,
       header=header,
       sep_save=args.sep_save,
       h5_compression_config=h5_compression_config,
      )
