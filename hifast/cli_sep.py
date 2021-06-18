#!/usr/bin/env python
# coding: utf-8
import os
import sys
import re
from glob import glob

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(allow_abbrev=False)
    parser.add_argument('fname',
                        help='file name; one of the chunk files')
    parser.add_argument('-f', '--force', action='store_true',
                        help='overwriting file if out file exists')
    parser.add_argument('-d', '--n_delay', type=int, required=True,
                       help='time of delay divided by sampling time')
    parser.add_argument('-m', '--n_on', type=int, required=True,
                       help='time of Tcal_on divided by sampling time')
    parser.add_argument('-n', '--n_off', type=int, required=True,
                        help='time of Tcal_off divided by sampling time')
    parser.add_argument('--frange', type=float, nargs=2,
                        help='freq range')
    parser.add_argument('--step', type=int, default=5,
                       help='number of files loaded into memory every time, default 5')
    parser.add_argument('--start', type=int,
                       help='chunk number of start')
    parser.add_argument('--stop', type=int,
                       help='chunk number of stop')
    parser.add_argument('--sep_save', action='store_true',
                       help='every step save to a file')
    parser.add_argument('--outdir', required=True,
                       help='the directory to store output files.')
    
    parser.add_argument('--smooth', default='mean', choices=['mean','poly','gaussian'],
                        help="smooth method: 'mean','poly',gaussian',.. default:mean")
    parser.add_argument('--s_sigma', type=float, default=5,
                        help='sigma for gaussian smooth, default 5MHz')
    parser.add_argument('--s_deg', type=int, default=1,
                        help='Degree of the fitting polynomial, 0 equals mean; 1 is linear fit. default 1')
    parser.add_argument('--dfactor', type=int,
                        help='Down-sample the data before cal Teff')
    parser.add_argument('--med_filter_size', type=int,
                        help='median filter kernel size for power of spec; odd number; default None')
    parser.add_argument('--noise_mode', default='high', choices=['high','low'],
                        help='noise_mode, high or low')
    parser.add_argument('--noise_date', default='20190115',
                        help='noise obs date, default auto')
    parser.add_argument('--med_filter_size_cal', type=int, default=5,
                        help='median filter kernel size for power of cal; odd number; default 5')
    parser.add_argument('--p_cal_fname',
                        help='input hdf5 power of cal file')
    parser.add_argument('--save_p_cal', action='store_true',
                        help='save power of cal to file')
    parser.add_argument('--not_cali', action='store_true',
                        help='save power of cal to file')
    
    args = parser.parse_args()
    
    #check file exits
    if not os.path.exists(args.fname):
        raise(OSError(f'File {args.fname} not exists.'))
    fname_part= re.sub('[0-9]{4}\.fits\Z', '', args.fname)
    #nB= int(re.findall(r'-M[0-1][0-9]',fname_part)[-1][2:])# used to read tcal
    n_delay, n_on, n_off = args.n_delay, args.n_on, args.n_off
    frange = args.frange
    step = args.step
    sep_save =args.sep_save
    outdir = args.outdir
    if outdir is not None:
        os.makedirs(outdir, exist_ok=True)
        
    ## check out file 
    def gen_out_name_base(fname_part, outdir):
        fname_add = os.path.basename(os.path.dirname(os.path.abspath(fname_part)))
        out_name_base = os.path.join(outdir, f"{os.path.basename(fname_part)[:-1]}-{fname_add}")
        return out_name_base
    out_name_base = gen_out_name_base(fname_part, outdir)
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
            
    start_all, stop_all= args.start, args.stop
    if start_all is None:
        start_all = 1
    if stop_all is None:
        stop_all = len(glob(fname_part+'*.fits'))
        
    smooth= args.smooth #"mean", "gaussian","poly"
    s_sigma= args.s_sigma
    s_deg= args.s_deg
    if smooth== "gaussian":
        if frange is not None:
            if (frange[1] - frange[0])/s_sigma <3:
                raise(ValueError('s_sigma is too larger for the input freq range'))
        s_para= {'s_sigma': s_sigma,}
    if smooth=='mean':
        smooth='poly'
        s_deg=0
    if smooth=='poly':
        s_para= {'s_deg': s_deg,}
    dfactor= args.dfactor
    med_filter_size= args.med_filter_size
    #ave_cycle= args.ave_cycle
    noise_mode = args.noise_mode
    noise_date = args.noise_date
    med_filter_size_cal = args.med_filter_size_cal
    p_cal_fname = args.p_cal_fname
    save_p_cal = args.save_p_cal
    
    print('processing: ', fname_part)
    print('n_delay, n_on, n_off: ', n_delay, n_on, n_off)
    print('freq range: ', frange)
    #record history
    from .utils.io import rec_his
    from .core.tcal_onoff import Tcal_onoff
    header = rec_his(args=args)
    spec = Tcal_onoff(fname_part=fname_part, n_delay=n_delay, n_on=n_on, n_off=n_off, 
                      start=start_all, stop=stop_all, 
                      frange=frange, verbose=True, 
                      smooth=smooth, s_para=s_para, dfactor=dfactor,
                      med_filter_size=med_filter_size, noise_mode=noise_mode, noise_date=noise_date,
                      med_filter_size_cal=med_filter_size_cal, p_cal_fname=p_cal_fname)
    spec(outdir=outdir, step=step, header=header, sep_save=sep_save, save_p_cal=save_p_cal, cali=(not args.not_cali))
    