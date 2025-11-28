

__all__ = []


from .utils.io import *
from glob import glob


sep_line = '##'+'#'*70+'##'
parser = ArgumentParser(prog=f"python -m hifast.{os.path.basename(sys.argv[0])[:-3]}",
                    formatter_class=formatter_class, allow_abbrev=False,
                    description='separate Cal on and Cal off; noise-diode Calibration', )

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

parser.add_argument('-d', '--n_delay', type=int, required=True,
                   help='time of delay divided by sampling time')
parser.add_argument('-m', '--n_on', type=int, required=True,
                   help='time of Tcal_on divided by sampling time')
parser.add_argument('-n', '--n_off', type=int, required=True,
                    help='time of Tcal_off divided by sampling time')
parser.add_argument('--frange', type=float, nargs=2,
                    help='freq range')
parser.add_argument('--ext_frange', type=bool_fun, choices=[True, False], default='True',
                    help="if True and smooth method is gaussian (`--smooth`), then extend 1*s_sigma in the beginning and end of freq")
parser.add_argument('--dfactor',
                    help='Down-sample the spectra freq resolution; str "W" or a int number')
parser.add_argument('--step', type=int, default=1,
                   help='number of files processed every time, default 1')
parser.add_argument('--start', type=int, default=1,
                   help='chunk number of start')
parser.add_argument('--stop', type=int,
                   help='chunk number of stop')
parser.add_argument('--sep_save', type=bool_fun, choices=[True, False], default='False',
                   help='every step save to a file')

group = parser.add_argument_group(f'*Power and Temperature of Cal smoothing\n{sep_line}')
group.add_argument('--smooth', choices=['mean','poly','gaussian'], default='gaussian',
                    help="smooth method: 'mean','poly', gaussian',.. default:gaussian")
group.add_argument('--s_sigma', type=float, default=5,
                    help='sigma for gaussian smooth (unit: MHz), default 5MHz')
group.add_argument('--s_deg', type=int, default=1,
                    help='Degree for `--smooth poly`, 0 equals mean; 1 is linear fit. default 1')
group.add_argument('--med_filter_size', type=int,
                    help='if larger than 0, apply median filter to spectra with kernel size as this nubmer; odd number; default None')
group.add_argument('--noise_mode', default='high', choices=['high','low'],
                    help='noise_mode, high or low')
group.add_argument('--noise_date', default='auto',
                    help='noise obs date, default auto')
group.add_argument('--med_filter_size_cal', type=int, default=5,
                    help='if larger than 0, apply median filter to power of Cal with this kernel size before smoothing; odd number; default 5')

group = parser.add_argument_group(f'*Checking Power of Cal (pcal)\n{sep_line}')
group.add_argument('--check_cal', choices=['none', 'A'], default='A',
                    help='check if Cal is on continuous source or contaminated by RFI, method `A` only support `-n` larger than 6.')
group.add_argument('--freq_step_c', type=float, default=4,
                    help='pcals will be checked on each freq bin with width of `--freq_step_c`. unit: MHz')
group.add_argument('--pcal_vary_lim_bin', type=float, default=0.01,
                    help='If the power change near a Cal is larger than ``pcal_vary_lim_bin``, the Cal in the freq bins will be recorded as bad and masked as nan.')
group.add_argument('--pcal_bad_lim_freq', type=float, default=0.5,
                    help='If bad fraction of a Cal in freq axis is larger than this limit, the entire Cal will be mask as nan.')

group = parser.add_argument_group(f'*Selecting pcal for spectra\n{sep_line}')
group.add_argument('--merge_pcals', type=bool_fun, choices=[True, False], default='False',
                    help='Whether merge pcal to calculate pcal with frequency.')


group = parser.add_argument_group(f'**`--merge_pcals False`: smoothing each pcal separately \n{sep_line}')
group.add_argument('--cal_dis_lim', type=float, default=1.6,
                   help='assign the nearest pcal to each spec, if the nearest pcal has been mask as nan, try to find further pcal with distance (measured with number of spec)'
                        'not exceed \n`--cal_dis_lim`*(`-m`+`-n`).')

group = parser.add_argument_group(f'**`--merge_pcals True`: merging all pcal to get a smoothing shape a pcal with freq \n{sep_line}')
group.add_argument('--method_merge', choices=['median', 'mean'], default='median',
                    help='method used to merge cal along time axis')
group.add_argument('--merge_cal_pre_process', choices=['scale', 'none'], default='scale',
                    help="if `scale`, scale the pcals to similar amplitude before merge")
group.add_argument('--calc_diff_method', choices=['div'], default='div',
                    help="used if '--merge_pcals' is not set as 'not'")
group.add_argument('--squeeze_diff_freq', choices=['median', 'mean',], default='median',
                    help="method applied to the difference of pcals with the 'merged-pcals' along freq to get relative amplitude")
group.add_argument('--squeeze_diff_freq_frange', type=float, nargs=2,
                    help='freq range used in `--squeeze_diff_freq`')
group.add_argument('--squeeze_diff_freq_bad_lim', type=float, default=0.5,
                    help='For an channel, if it fraction of nan in pcals remained after `--check_cal` large than this limit, it will not be used in `--squeeze_diff_freq`.')
group.add_argument('--method_interp', default='quadratic',
                    help="Interpolate the relative amplitude to get pcal for each spectrum."
                         "Choose from 'gaussian', 'slinear', 'quadratic', 'cubic', 'nearest' \nor 'poly1d', poly2d, ..., 'polynd'.\n"
                         "'gaussian' is gaussian smooth; 'polynd' is poly fitting")
group.add_argument('--method_interp_sigma_t', type=float, default=300,
                    help='sigma_t in  ``--method_interp gaussian``')
group.add_argument('--method_interp_edges', choices=['nearest', 'extrapolate'], default='nearest',
                    help="how to deal the edge spetra in interpolating")

# parser.add_argument('--p_cal_fname',
#                     help='input hdf5 power of cal file')
parser.add_argument('--save_pcals', type=bool_fun, choices=[True, False], default='True',
                    help='if True, save power of cal to file')

parser.add_argument('--not_cali', type=bool_fun, choices=[True, False], default='False',
                    help='if True, save power data, no temperature calibration')


if __name__ == '__main__':
    args_ = parser.parse_args()
    print('#'*35+'Args'+'#'*35)
    print(parser.format_values())
    print('#'*35+'####'+'#'*35)

    args = args_

    #check input file exits
    if not os.path.exists(args.fpath):
        raise(OSError(f'File {args.fpath} not exists.'))

    # replace patten in outdir
    nB = get_nB(args.fpath)
    project = get_project(args.fpath)
    ## use the dirname of the fits file as "date"
    date = os.path.basename(os.path.dirname(os.path.abspath(args.fpath)))
    args.outdir = sub_patten(args.outdir, date=date, nB=f'{nB:02d}', project=project)
    # expand '~' as outdir may be a string in bash
    args.outdir = os.path.expanduser(args.outdir)
    print(f'outdir: {args.outdir}')
    if args.outdir is not None:
        if not os.path.exists(args.outdir):
            print(f'outdir {args.outdir} not exists. Create it now')
            os.makedirs(args.outdir, exist_ok=True)

    ## add pre info to the print output
    from .utils.output import set_output
    pre_str = f'[hifast.{os.path.basename(sys.argv[0])[:-3]}]['
    try:
        pre_str += project
    except:
        pass
    try:
        pre_str += f"-M{nB:02d}"
    except:
        pass
    try:
        pre_str += f"-{date}"
    except:
        pass
    pre_str += '] '
    set_output(pre_str)

    ## check out file
    fname_add = date
    fname_part = re.sub('[0-9]{4}\.fits\Z', '', args.fpath)
    fname_part = re.sub('[0-9]{4}\.hdf5\Z', '', fname_part)
    out_name_base = os.path.join(args.outdir, f"{os.path.basename(fname_part)[:-1]}-{fname_add}")
    fileout =  out_name_base + f"-specs_T.hdf5"
    fileout_sep = glob(out_name_base + rf"-specs_T_[0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9].hdf5")

    if os.path.exists(fileout) or len(fileout_sep)>0:
        if args.force:
            print(f"will overwrite the existing out file")
        else:
            print(f"File exists {fileout}")
            print(fileout_sep)
            print('exit... Using -f to overwrite it.')
            sys.exit()

    # check frange
    frange = args.frange
    if frange is not None and args.smooth == "gaussian" and args.ext_frange:
        frange = [frange[0] - args.s_sigma, frange[1] + args.s_sigma]
    args.frange = frange
    # chunk files
    if args.start is None:
        args.start = 1

    # cal smooth
    smooth = args.smooth #"mean", "gaussian","poly"
    s_sigma = args.s_sigma
    s_deg = args.s_deg
    if smooth == "gaussian":
        if args.frange is not None:
            if (args.frange[1] - args.frange[0])/s_sigma <3:
                raise(ValueError('s_sigma is too larger for the input freq range'))
        s_para = {'s_sigma': s_sigma,}
    if smooth == 'mean':
        smooth = 'poly'
        s_deg = 0
    if smooth == 'poly':
        s_para = {'s_deg': s_deg,}

    dfactor = args.dfactor
    if dfactor is not None and dfactor.upper() != "W":
        dfactor = int(dfactor)

    print('processing: ', fname_part)
    print('n_delay, n_on, n_off: ', end='')
    print(args.n_delay, args.n_on, args.n_off, sep=', ')
    print('freq range: ', args.frange)

    # import
    from .core.cal import CalOnOff
    from .core.cal2 import CalOnOffM
    from .core.cal2 import CalOnOff1111
    import json

    #record history
    header = rec_his(args=json.dumps(args.__dict__))

    # cal
    keys = ['n_delay', 'n_on', 'n_off',
            'start', 'stop',
            'frange',
            'dfactor',
            'med_filter_size',
            'noise_mode',
            'noise_date',
            'med_filter_size_cal',
            ]
    paras = {key:getattr(args, key) for key in keys}
    paras['fname_part'] = args.fpath # use full
    paras['verbose'] = True
    paras['smooth'] = smooth
    paras['s_para'] = s_para


    if args.check_cal == 'none' or args.not_cali:
        Cal_cls = CalOnOff
    elif args.check_cal == 'A':
        paras['freq_step_c'] = args.freq_step_c
        paras['pcal_vary_lim_bin'] = args.pcal_vary_lim_bin
        paras['pcal_bad_lim_freq'] = args.pcal_bad_lim_freq
        if args.merge_pcals:
            # check
            if args.method_interp in ['gaussian', 'slinear', 'linear', 'quadratic', 'cubic', 'nearest', 'next', 'previous']:
                pass
            elif args.method_interp.startswith('poly'):
                try:
                    int(args.method_interp[4:-1])
                except ValueError:
                    raise(ValueError(f'`--method_interp {args.method_interp}` not support'))
            else:
                raise(ValueError(f'`--method_interp {args.method_interp}` not support'))
            Cal_cls = CalOnOffM
        else:
            Cal_cls = CalOnOff1111
    if args.merge_pcals and args.check_cal == 'none':
        raise(ValueError('if `--merge_pcals True`, `--check_cal` must not be `none`'))
    # init
    spec = Cal_cls(**paras)
    # set para added
    if args.check_cal == 'A' and (not args.not_cali):
        if args.merge_pcals:
            spec.set_para_pcals(calc_diff_method=args.calc_diff_method,
                                squeeze_diff_freq=args.squeeze_diff_freq,
                                squeeze_diff_freq_bad_lim=args.squeeze_diff_freq_bad_lim,
                                squeeze_diff_freq_frange=args.squeeze_diff_freq_frange,
                                method_merge=args.method_merge,
                                merge_cal_pre_process=args.merge_cal_pre_process,
                                method_interp=args.method_interp,
                                method_interp_edges=args.method_interp_edges,
                                method_interp_sigma_t=args.method_interp_sigma_t,
                               )
        else:
            spec.set_para_pcals(cal_dis_lim=args.cal_dis_lim,
                                )

    spec(outdir=args.outdir, step=args.step, header=header,
         sep_save=args.sep_save, save_pcals=args.save_pcals, cali=(not args.not_cali))
