#!/usr/bin/env python
# coding: utf-8

import sys
import os
import re

if __name__ == '__main__':
    import warnings 
    warnings.filterwarnings("ignore",r'overflow encountered in exp')
    import argparse
    parser = argparse.ArgumentParser(allow_abbrev=False)
    parser.add_argument('fname',
                        help='file name')
    parser.add_argument('--nproc', type=int,
                        help='number process')
    parser.add_argument('-f', '--force', action='store_true',
                        help='overwriting file if out file exists')
    parser.add_argument('--frange', type=float, nargs=2,
                       help='freq range')
    parser.add_argument('--nB_radec', type=int, default=1, 
                       help='Beam number in the radec file name')
    parser.add_argument('--outdir',
                       help='default is same with the input file')
    parser.add_argument('--flux', action='store_true',
                        help='convert temperature to flux before subtract baseline')
    parser.add_argument('-c', '--cali_fname',
                       help='quasar calibration file name')
    
    parser.add_argument('--method', default='arPLS', choices=['arPLS', 'srPLS', 'Chebyshev', 'poly', 'sin_poly', 'sin_poly_2', 'S', 'SP', 'original'],
                       help='method to fit baseline; arPLS, srPLS, 0,')
    
    parser.add_argument('--lam', type=float, default=1.0e8, 
                       help='baseline fit parameters')
    parser.add_argument('--deg', type=int, default=2, 
                       help='baseline fit parameters')
    parser.add_argument('--offset', type=float, default=2, 
                       help='baseline fit parameters')
    parser.add_argument('--ratio', type=float, default=0.01, 
                       help='baseline fit parameters')
    parser.add_argument('--niter', type=int, default=100, 
                       help='baseline fit parameters')
    parser.add_argument('--sin_f', type=float, nargs='+',
                       help='sin freq')
    parser.add_argument('--s_method_freq', choices=['median', 'gaussian', 'boxcar', 'PLS'],
                       help='')
    parser.add_argument('--s_sigma_freq', type=int,
                       help='')
    parser.add_argument('--average_every_freq', type=int,
                       help='')
    
    parser.add_argument('--lam_2', type=float, default=1.0e12, 
                       help='baseline fit parameters')
    parser.add_argument('--deg_2', type=int, default=2, 
                       help='baseline fit parameters')
    parser.add_argument('--offset_2', type=float, default=2, 
                       help='baseline fit parameters')
    parser.add_argument('--ratio_2', type=float, default=0.01, 
                       help='baseline fit parameters')
    parser.add_argument('--niter_2', type=int, default=100, 
                       help='baseline fit parameters')
    parser.add_argument('--sin_f_2', type=float, nargs='+',
                       help='sin freq')
    parser.add_argument('--s_method_freq_2', choices=['median', 'gaussian', 'boxcar', 'PLS'],
                       help='')
    parser.add_argument('--s_sigma_freq_2', type=int,
                       help='')
    parser.add_argument('--average_every_freq_2', type=int,
                       help='')
    
    parser.add_argument('--njoin_t', type=int,
                       help='join nspec along t')
    parser.add_argument('--njoin_t_2', type=int,
                       help='join nspec along t ')
    parser.add_argument('--s_method_t', choices=['gaussian', 'boxcar', 'median_filter'],
                       help='')
    parser.add_argument('--s_sigma_t', type=int, default=5,
                       help='')
    parser.add_argument('--s_method_t_2', choices=['gaussian', 'boxcar', 'median_filter'],
                       help='')
    parser.add_argument('--s_sigma_t_2', type=int, default=5,
                       help='')
    parser.add_argument('-T', '--trans', action='store_true',
                       help='')
    parser.add_argument('--exclude_m', type=int, default=0,
                       help='')
    parser.add_argument('--no_radec', action='store_true',
                       help='')
    
    args = parser.parse_args()
    nproc = args.nproc
    file_spec = args.fname
    outdir = args.outdir
    trans= args.trans
    flux = args.flux
    
    fpart = '-flux_' if args.flux else '-'
    fpart += 'bld'
    if args.trans: fpart += '_T'
    if outdir is None: outdir = os.path.dirname(file_spec)
    fileout = os.path.join(outdir, '.'.join(os.path.basename(file_spec).split('.')[:-1]) + f'{fpart}.hdf5')
    if os.path.exists(fileout):
        if args.force:
            print(f"will overwrite the existing out file {fileout}")
        else:
            print(f"File exists {fileout}")
            print('exit... Using -f to overwrite it.')
            sys.exit()
    exclude_m = args.exclude_m
    no_radec = args.no_radec
    
import numpy as np
import h5py
from scipy import ndimage
from .core.baseline import get_baseline, get_baseline_mp
from .utils.misc import extend_Trues, boxcar_smooth1d, median_filter_1d

def gen_radec_file(file_spec, paras):
    import subprocess
    process = subprocess.Popen([sys.executable, '-m' , 'hifast.cli_radec', file_spec] + paras, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #process.wait()
    print(*process.communicate())
    if process.returncode !=0:
        raise(ValueError(f'fail to generate the radec file of {file_spec}'))
def get_radec(file_spec, nB, nB_radec, mjd):

    Replace_nB= lambda path,nB: os.path.join(os.path.dirname(path), re.sub(r'-M[0-1][0-9]',f"-M{nB:02d}",os.path.basename(path)))
    file_radec= '.'.join(file_spec.split('.')[:-1])+'-radec.hdf5'
    file_radec= Replace_nB(file_radec, nB_radec)

    f = h5py.File(file_radec,'r')
    #check if mjd match
    mjd_match= False if len(mjd)!=len(f['mjd'][:]) else (mjd==f['mjd'][:]).all()
    if not mjd_match:
        print('the mjd of spec file is not same with that in ',file_radec)
        sys.stdout.flush()
        file_radec_2= Replace_nB(file_radec, nB)
        if not os.path.exists(file_radec_2):
            print('try to generate ', file_radec_2)
            sys.stdout.flush()
            header = f['Header'].attrs
            import json
            for key in header.keys():
                argv = json.loads(header[key])['argv']
                if 'cli_radec' not in argv:
                    continue
                paras = []
                for s in argv.split():
                    if s[:1] == '-':
                        paras += [s,]
            gen_radec_file(file_spec, paras)
            sys.stdout.flush()
        f.close()
        f= h5py.File(file_radec_2,'r')
        print('using', file_radec_2)
        sys.stdout.flush()
    #check again
    mjd_match= False if len(mjd)!=len(f['mjd'][:]) else (mjd==f['mjd'][:]).all()
    if not mjd_match:
        raise(ValueError(f'the mjd of spec file is not same with that in {file_radec_2} Abort...'))
    #load     
    ra= f['ra'+'%d'%nB][:]
    dec= f['dec'+'%d'%nB][:]
    if 'is_extrapo' in f.keys():
        is_extrapo = f['is_extrapo'][:]
    else:
        is_extrapo = None
    f.close()
    return ra, dec, is_extrapo

def array_split(arr, size, base=None):
    if base is None:
        return np.array_split(arr, size)
    else:
        len_t = len(arr)
        n_ = int(np.ceil(len_t/base))
        split_ = base*np.cumsum(np.array([n_//size+1,]*(n_%size) + [n_//size,]*(size- n_%size)))[:-1]
        return np.array_split(arr, split_)

class TEMPFile(object):
    def __init__(self, outdir):
        import uuid
        tmpdir = os.path.join(outdir,'tmp')
        os.makedirs(tmpdir, exist_ok=True)
        ftmp_part = os.path.join(tmpdir,f"pid_{os.getpid()}_{uuid.uuid4().hex}")
        ftmp_part += '_rank_'
        self.ftmp_part = ftmp_part
        
    def save(self, rank=None, **kwargs):
        if 'rank' in kwargs.keys():
            raise(ValueError("can't save variable named as 'rank'"))
        for key in kwargs.keys():
            if rank is None:
                [np.save(self.ftmp_part + f'{i}_{key}.npy', data) for i, data in enumerate(kwargs[key])]
            else:
                np.save(self.ftmp_part + f'{rank}_{key}.npy', kwargs[key])
    def load(self, rank, name, rm=True):
        fname = self.ftmp_part + f'{rank}_{name}.npy'
        res = np.load(fname)
        if rm:
            os.remove(fname)
        return res
    
if __name__ == '__main__':
    frange = args.frange
    nB_radec = args.nB_radec
    cali_fname = args.cali_fname
    
    method = args.method
    njoin_t = args.njoin_t
    njoin_t_2 = args.njoin_t_2
    s_method_t = args.s_method_t
    s_sigma_t = args.s_sigma_t
    s_method_t_2 = args.s_method_t_2
    s_sigma_t_2 = args.s_sigma_t_2
    #smooth before fit baseline
    s_method_freq = args.s_method_freq
    s_sigma_freq = args.s_sigma_freq
    average_every_freq = args.average_every_freq
    s_method_freq_2 = args.s_method_freq_2
    s_sigma_freq_2 = args.s_sigma_freq_2
    average_every_freq_2 = args.average_every_freq_2
    
    #fit args
    fit_args={}
    if method in ['arPLS', 'srPLS', 'S', 'SP']:
        fit_args['lam'] = args.lam
    if method == 'sin_poly' or method == '1':
        fit_args['f'] = args.sin_f
        fit_args['rew'] = False  # only fit once
    fit_args['offset'] = args.offset
    fit_args['deg'] = args.deg
    fit_args['ratio'] = args.ratio
    fit_args['niter'] = args.niter
    
    if method in ['S', 'SP']:
        fit_args_2 = {}
        fit_args_2['f'] = args.sin_f_2
        fit_args_2['offset'] = args.offset_2
        fit_args_2['deg'] = args.deg_2
        fit_args_2['ratio'] = args.ratio_2
        fit_args_2['niter'] = args.niter_2
        fit_args_2['rew'] = False  # only fit once
    
     
    nB = int(re.findall(r'-M[0-1][0-9]',file_spec)[-1][2:])

    fs = h5py.File(file_spec,'r')
    mjd = fs['mjd'][()]
    if not no_radec:
        if 'ra' not in fs.keys():
            ra, dec, is_extrapo = get_radec(file_spec, nB, nB_radec, mjd)
        else:
            ra = fs['ra'][()]
            dec = fs['dec'][()]
            is_extrapo = None
    if 'T' in fs.keys():
        T = fs['T'][()]
        outfield = 'Ta'
    elif 'Ta' in fs.keys():
        T = fs['Ta'][()]
        outfield = 'Ta'
    elif 'flux' in fs.keys():
        T = fs['flux'][()]
        outfield = 'flux'
        if flux:
            raise(ValueError(f'flux already exists, please remove \"--flux\"'))
    if flux:
        outfield = 'flux'


    ind_sort = np.argsort(mjd)
    mjd = mjd[ind_sort]
    if not no_radec:
        ra = ra[ind_sort]
        dec = dec[ind_sort]
    T = T[ind_sort]
    freq = fs['freq'][:]

    if frange is not None:
        is_ = (freq >= frange[0]) & (freq <= frange[1])
        freq = freq[is_]
        T = T[:,is_]

    # try to load header
    if 'Header' in fs.keys():
        from collections import OrderedDict
        header_in= OrderedDict(fs['Header'].attrs.items())
        try:
            header_in.update(OrderedDict(fs['Header'].attrs.items()))
        except:
            print('radec file has no header')
    else:
        header_in = None

    #close files
    # fs.close() # don't close, load extra later

    if trans:
        T = T.transpose((1,0,2))
        if flux:
            raise()
#     if smooth_t:
#         T = ndimage.gaussian_filter1d(T, smooth_t_sigma, axis=0)
     
  
    if flux:
        from .flux import cali_src
        print('Flux calibrating ...')
        T = cali_src(T, nB, freq, cali_fname, ra=ra, dec=dec, mjd=mjd)
    
    ## subtract baseline
    verbose = True
    
    def sub(yss, njoin, method, para, exclude_fun=None, s_method_t=None, s_sigma_t=3):
        #global freq, verbose, nproc
        yss_ori = yss
        if s_method_t == 'gaussian':
            yss = ndimage.gaussian_filter1d(yss, s_sigma_t, axis=0)
        elif s_method_t == 'boxcar':
            yss = boxcar_smooth1d(yss, s_sigma_t, axis=0)
        elif s_method_t == 'median_filter':
            yss = median_filter_1d(yss, axis=0, size=s_sigma_t)

        import copy
        para = copy.deepcopy(para)
        def apply(yss):
            if njoin is not None:
                shape_add = yss.shape[1:]
                yss = yss.reshape((-1,njoin) + shape_add)
                yssm = np.mean(yss, axis=1)
            else:
                yssm = yss
            if exclude_fun is not None:
                exclude = exclude_fun(yssm)
            else:
                exclude = None
            if para['s_method'] == 'PLS':
                # use arPLS, lam from s_sigma ( s_sigma_freq)
                yssm = get_baseline_mp(nproc, freq, yssm, axis=1, method='arPLS', bl_para={'lam': para['s_sigma'], "offset":2, 'deg':2}, verbose=verbose)
                para.pop('s_method')
            bls = get_baseline_mp(nproc, freq, yssm, axis=1, method=method, verbose=verbose, exclude=exclude, **para)
            
            if njoin is not None: 
                return (yss - bls[:, None, :, :]).reshape((-1,) + shape_add)
            else:
                return yss_ori - bls
            
        if njoin is not None:
            num = len(yss)//njoin*njoin
            if num != len(yss):
                return np.vstack([apply(yss[:num]), apply(yss[-njoin:])[num-len(yss):]])
            else:
                return apply(yss)
        else:
            return apply(yss)
            
    
    if method == 'S' or method == 'SP':
        method_a = 'arPLS'
        para = {
        's_sigma': s_sigma_freq,
        's_method': s_method_freq,
        'average_every': average_every_freq,
        'bl_para': fit_args,
               }
        
        method_b = 'sin_poly'
        bounds = [(0.,1.), (0.899, 0.961), (0, 2*np.pi), (-1,1)]
        fit_args_2['opt_para'] = {'bounds':bounds,}
        para2 = {
        's_sigma': s_sigma_freq_2,
        's_method': s_method_freq_2,
        'average_every': average_every_freq_2,
        'bl_para': fit_args_2,
               }
        T_bld1 = sub(T, njoin_t, method_a, para=para, s_method_t=s_method_t, s_sigma_t=s_sigma_t)
        #exclude_fun = lambda x:abs(x) > 1.2*np.diff(np.percentile(x, [16,84], axis=1), axis=0)[0][:,None,:]
        if exclude_m == 0:
            exclude_fun = lambda x: extend_Trues(abs(x) > 1.*np.diff(np.percentile(x, [16,84], axis=1), axis=0)[0][:,None,:], axis=1, ext_frac=1/3, ext_add=3)
        elif exclude_m == 1:
            exclude_fun = lambda x: extend_Trues(abs(x) > 2.5*np.min(np.diff(np.percentile(x, [10, 50, 90], axis=1), axis=0), axis=0)[:,None,:], axis=1, ext_frac=1/3, ext_add=3)
        bl_sin = T_bld1 - sub(T_bld1, njoin_t_2, method_b, para=para2, exclude_fun=exclude_fun, s_method_t=s_method_t_2, s_sigma_t=s_sigma_t_2)
        del(T_bld1)
        T = T - bl_sin
        del(bl_sin)
        if method == 'SP':
            T = sub(T, njoin_t, method_a, para=para)
    elif method == '1':
        pass
    else:
        para = {
        's_sigma': s_sigma_freq,
        's_method': s_method_freq,
        'average_every': average_every_freq,
        'bl_para': fit_args,
               }
        T = sub(T, njoin_t, method, para=para, s_method_t=s_method_t, s_sigma_t=s_sigma_t)
        
    if trans:
        T = T.transpose((1,0,2))
    print(f"Saving...")
    dict_out= {}
    dict_out['mjd'] = mjd
    if not no_radec:
        dict_out['ra'] = ra
        dict_out['dec'] = dec
        if is_extrapo is not None:
            dict_out['is_extrapo'] = is_extrapo
    dict_out[outfield] = T.astype('float32')
    dict_out['freq'] = freq
    
    #save file
    from .utils.io import rec_his, save_dict_hdf5, add_extra
    add_extra(fs, dict_out)
    fs.close()
    header=rec_his(args=args)
    if header_in is not None: header.update(header_in)
    save_dict_hdf5(fileout, dict_out, header=header)
    print(f"Saved to {fileout}")
        