

__all__ = ['parser', 'args', 'fname', 'obsfile', 'frange', 'noise_mode', 'noise_date', 'sigma', 'obsmode', 'modes',
           'CalinBs', 'outdir', 'calname', 'ccrd', 'fpparas', 't_src', 'ncir', 'saveT', 'processing_nB', 'outdata',
           'nBs_have']


from . import FLUXGAIN as cbr
from ..utils.io import save_dict_hdf5, replace_nB
from ..utils.io import bool_fun
from ..utils import obs_log
import numpy as np
import os
import sys
import h5py


parser = cbr.parser
parser.add_argument('--nproc', '--np', type=int, default=3,
                   help='number of process used')
parser.add_argument('--cache', type=bool_fun, choices=[True, False], default='False',
                       help='If set to True, will cache/load the results of all beams except for the initial one.')
parser.add_argument('--plot', action='store_true',
                       help='plot FluxGain with frequency.')
parser.add_argument('--obs_log',
                       help='If set, will use the obs_log to get `-d`, `-m`, `-n`, `--noise_mode`.\
                        It requires a file path ending in `.log` or `.txt`, or a directory path where the log file can be searched for.')


args = parser.parse_args()


if args.obs_log is not None:
    if not os.path.exists(args.obs_log):
        raise ValueError(f'`--obs_log` not exists: `{args.obs_log}`.')
    if os.path.isfile(args.obs_log):
        obs_log_path = args.obs_log
    elif os.path.isdir(args.obs_log):
        print(f'try to search for log file in the directory: `{args.obs_log}`.')
        obs_log_path = obs_log.search_log_fpath(args.fname, args.obs_log)
    else:
        raise ValueError(f'`--obs_log` should be a file path ending in `.log` or `.txt`, or a directory path where the log file can be searched for.\
        \nBut got `{args.obs_log}`.')
    print(f'Use obs_log file: `{obs_log_path}`.')
    paras = obs_log.get_sep_para_from_log(obs_log_path)[0]
    print(f'got paras: {paras}')
    args.d = paras['n_delay']
    args.m = paras['n_on']
    args.n = paras['n_off']
    args.noise_mode = paras['noise_mode']
elif args.d is None or args.m is None or args.n is None or args.noise_mode is None:
    raise ValueError(f'`--obs_log` is not set, but `--d`, `--m`, `--n`, `--noise_mode` are not set either.')


fname= args.fname
obsfile= fname.split('/')[-1]

d, m, n = args.d, args.m, args.n
frange= args.frange
noise_mode = args.noise_mode
noise_date = args.noise_date
sigma= args.smt_sigma
obsmode= args.obsmode

# may remove
modes= ['Drift','MultiBeamCalibration','OnOff','MultiBeamOTF','DriftWithAngle','DecDriftWithAngle']
if obsmode not in modes:
    print('Please check observation mode.')
    os._exit(0)
# end

CalinBs= args.nBs
if CalinBs== None:
    if obsmode in ['MultiBeamCalibration','MultiBeamOTF',]:
        nBs= np.arange(1,20)
    else:
        nBs= [1,]
else:
    nBs= CalinBs
print('Flux Calibration of ' + replace_nB(obsfile, 0), end='')
print(' of Beams ', nBs)


outdir= args.outdir
if outdir==None:
    outdir= '.'
calname= args.calname
ccrd = cbr.load_crd(calname,args.crd)
print('Calibrater position:', ccrd)
fpparas= args.fluxProfilePara

# for tracking
t_src= args.t_src/3600./24#change into days
ncir= args.n_cir
#T_select=args.T_select
#rlim_src = args.rlim_src
saveT= args.saveT



# set cbr.* global variable. Need improve later
cbr.fname = fname
cbr.d = d
cbr.m = m
cbr.n = n
cbr.frange = frange
cbr.noise_mode = noise_mode
cbr.noise_date = noise_date
#cbr.oname = oname
cbr.obsmode = obsmode
cbr.sigma = sigma
cbr.saveT = saveT
cbr.obsmode = obsmode
cbr.nBs = nBs
cbr.ccrd = ccrd
cbr.fpparas = fpparas
cbr.calname = calname
cbr.outdir = outdir
cbr.t_src = t_src
cbr.ncir = ncir
#cbr.T_select = T_select
#cbr.rlim_src = rlim_src
cbr.args = args


def processing_nB(nB, check_cache=False, outname_pattern=None):
    if check_cache:
        cache_name = outname_pattern + f'-M{nB:02d}-cache.h5'
        if os.path.exists(cache_name):
            print(f'Load FluxGain of beam {nB:02d} from cached file: {cache_name}')
            outdata = {k:v[()] for k,v in h5py.File(cache_name,'r').items()}
            for key in outdata.keys():
                if isinstance(outdata[key], bytes):
                    outdata[key] = outdata[key].decode('utf-8')
            return outdata
    outdata = {}
    i = nB
    cbr.nB = nB

    try:
        K_Jy,Tcal_s,tcal_file,mjd,ONp,OFFp,TMAX= cbr.load_data(fname)
    except SystemExit as e:
        print(f'Fail to process Beam {nB:02d}')
        return None

    print(fname,i)
    radec= cbr.radec
    ra,dec,ZD = radec[f'ra{i}'], radec[f'dec{i}'],radec[f'ZD{i}']
    outdata[f'M{nB:02d}']= K_Jy
    outdata[f'Tcal{nB}']= Tcal_s
    outdata[f'ra{i}']= ra
    outdata[f'dec{i}']= dec
    outdata[f'ZD{i}']= ZD
    outdata[f'ONp{i}']= ONp
    outdata[f'OFFp{i}']= OFFp
    outdata[f'Tmax{i}']= TMAX
    print(f'Finished calibration of M{nB:02d}, {cbr.fit_k*100:.3f}% input data has been ignored.')
    print('#'*50)

    # common
    outdata['Tcal_file']= tcal_file
    print('Used Tcal of '+tcal_file)
    outdata['freq']= cbr.freq
    outdata['cal_flux']= cbr.flux
    outdata['mjd']= mjd
    #
    outdata['outname'] = cbr.outname
    if check_cache:
        print(f'Save FluxGain of beam {nB:02d} to cache file {cache_name}')
        save_dict_hdf5(cache_name, outdata)
    return outdata


# process the first one, M01 in general, and the radec are generated
cbr.nproc = args.nproc
print(f'Process Beam {nBs[0]:02d}')
outdata = processing_nB(nBs[0])
if outdata is None:
    print('Fail to process the first beam in ``--nBs``. Abort.')
    sys.exit(1)
nBs_have = 1
if len(nBs) > 1:
    print(f'Process rest Beams with {args.nproc} CPUs')
    if args.cache:
        def processing_nB_with_cache(nB): return processing_nB(nB,True,cbr.outname)
    else:
        processing_nB_with_cache = processing_nB
    from multiprocessing import Pool
    with Pool(args.nproc) as pool:
        result = pool.map(processing_nB_with_cache, np.copy(nBs).tolist()[1:])

    for j, d in enumerate(result):
        if d is not None:
            outdata.update(d)
            nBs_have += 1
        else:
            print(f'Fail to process Beam {nBs[1:][j]:02d}')


if nBs_have == 19:
    print("19 beams have been successfully processed and saving to file with name containing 'All'.")
    outname = outdata['outname']+'-FluxGain-All.hdf5'
else:
    print("Not all 19 beams were processed successfully and saved to a file that doesn't have 'All' in its name.")
    outname= outdata['outname'] + '-FluxGain.hdf5'
save_dict_hdf5(outname,outdata)
print('Saved to '+ outname)
if args.plot:
    from .plot import plot
    fig = plot(outdata, suptitle=outname)
    fig.savefig(outname.split('.', 1)[0] + '.png')
