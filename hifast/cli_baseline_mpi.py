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
    
    parser.add_argument('--method', default='arPLS', choices=['arPLS', 'srPLS', 'Chebyshev', 'poly', 'sin_poly', 'sin_poly_2', '0', '1'],
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
    parser.add_argument('--s_method_freq', choices=['median', 'gaussian', 'boxcar'],
                       help='')
    parser.add_argument('--s_sigma_freq', type=int, default=3,
                       help='')
    parser.add_argument('--average_every_freq', type=int, default=2,
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
    parser.add_argument('--s_method_freq_2', choices=['median', 'gaussian', 'boxcar'],
                       help='')
    parser.add_argument('--s_sigma_freq_2', type=int, default=3,
                       help='')
    parser.add_argument('--average_every_freq_2', type=int, default=2,
                       help='')
    
    parser.add_argument('--njoin_t', type=int,
                       help='join nspec along t')
    parser.add_argument('--smooth_t',
                       help='')
    parser.add_argument('--smooth_t_sigma', type=int, default=5,
                       help='')
    parser.add_argument('-T', '--trans', action='store_true',
                       help='')
    
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    
    args = parser.parse_args()
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
            if rank == 0:
                print(f"will overwrite the existing out file {fileout}")
            pass
        else:
            if rank == 0:
                print(f"File exists {fileout}")
                print('exit... Using -f to overwrite it.')
            sys.exit()
    
import numpy as np
import h5py
from scipy import ndimage
from .baseline import get_baseline

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
    smooth_t = args.smooth_t
    smooth_t_sigma = args.smooth_t_sigma
    #smooth before fit baseline
    s_method_freq = args.s_method_freq
    s_sigma_freq = args.s_sigma_freq
    average_every_freq = args.average_every_freq
    s_method_freq_2 = args.s_method_freq_2
    s_sigma_freq_2 = args.s_sigma_freq_2
    average_every_freq_2 = args.average_every_freq_2
    
    #fit args
    fit_args={}
    if method in ['arPLS', 'srPLS', '0']:
        fit_args['lam'] = args.lam
    if method == 'sin_poly' or method == '1':
        fit_args['f'] = args.sin_f
        fit_args['rew'] = False  # only fit once
    fit_args['offset'] = args.offset
    fit_args['deg'] = args.deg
    fit_args['ratio'] = args.ratio
    fit_args['niter'] = args.niter
    
    if method == '0':
        fit_args_2 = {}
        fit_args_2['f'] = args.sin_f_2
        fit_args_2['offset'] = args.offset_2
        fit_args_2['deg'] = args.deg_2
        fit_args_2['ratio'] = args.ratio_2
        fit_args_2['niter'] = args.niter_2
        fit_args_2['rew'] = False  # only fit once
    
     
    nB = int(re.findall(r'-M[0-1][0-9]',file_spec)[-1][2:])

    if rank == 0:
        fs = h5py.File(file_spec,'r')
        mjd = fs['mjd'][()]
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
        if smooth_t:
            T = ndimage.gaussian_filter1d(T, smooth_t_sigma, axis=0)
        ## for mpi 
        # store some variable first
        mjd_ori = mjd
        ra_ori = ra
        dec_ori = dec
        freq_ori = freq
        if trans:
            val_range = None
            mjd=np.arange(len(T))
        # split
        T = array_split(T, size, njoin_t)
        mjd = array_split(mjd, size, njoin_t)
        if flux:
            ra= np.array_split(ra, size)
            dec= np.array_split(dec, size)
        print('Scattering T')
        sys.stdout.flush()
        ##save large variable into tmp file
        tempfile = TEMPFile(outdir)
        tempfile.save(T=T)
    else:
        #mpi
        freq = None
        tempfile = None
        mjd = None
        ra = None
        dec = None
    # broadcast val range in case they are changed.
    frange = comm.bcast(frange, root=0)
    flux = comm.bcast(flux, root=0)
    freq = comm.bcast(freq, root=0)
    tempfile = comm.bcast(tempfile, root=0)
    mjd = comm.scatter(mjd, root=0)
    T = tempfile.load(rank, 'T', rm=True)
    
    if flux:
        ra = comm.scatter(ra, root=0)
        dec = comm.scatter(dec, root=0)
        from .flux import cali_src
        if rank == 0: print('Flux calibrating ...')
        T = cali_src(T, nB, freq, cali_fname, ra=ra, dec=dec, mjd=mjd)
        if rank == 0: 
            print('Calibrated')
            sys.stdout.flush()
    ## subtract baseline
    verbose = True if rank==0 else False
    def proc(T):
        if njoin_t is not None:
            shape_add = T.shape[1:]
            T = T.reshape((-1,njoin_t) + shape_add)
            Tm = np.mean(T, axis=1)
        else:
            Tm = T

        if method == '0':
            method_a = 'arPLS'
            method_b = 'sin_poly'
            Tmres = Tm
            bls2 = get_baseline(freq, Tmres, axis=1, s_sigma=s_sigma_freq, s_method=s_method_freq, average_every=average_every_freq, method=method_a, bl_para=fit_args, verbose = verbose)
            Tmres = Tmres-bls2
            exclude = abs(Tmres) > np.diff(np.percentile(Tmres, [16,84], axis=1), axis=0)[0][:,None,:]
            bls2 = get_baseline(freq, Tmres, axis=1, s_sigma=s_sigma_freq_2, s_method=s_method_freq_2, average_every=average_every_freq_2, exclude=exclude, method=method_b, bl_para=fit_args_2, verbose = verbose)
            Tmres = Tmres-bls2
            bls2 = get_baseline(freq, Tmres, axis=1, s_sigma=s_sigma_freq, s_method=s_method_freq, average_every=average_every_freq, method=method_a, bl_para=fit_args, verbose = verbose)
            Tmres = Tmres-bls2
            bls = Tm - Tmres
        elif method == '1':
            exclude = abs(Tm) > np.diff(np.percentile(Tm, [16,84], axis=1), axis=0)[0][:,None,:]
            bls = get_baseline(freq, Tm, axis=1, s_sigma=s_sigma_freq, s_method=s_method_freq, average_every=average_every_freq, exclude=exclude, method='sin_poly', bl_para=fit_args, verbose = verbose)
        else:
            bls = get_baseline(freq, Tm, axis=1, s_sigma=s_sigma_freq, s_method=s_method_freq, average_every=average_every_freq, method=method, bl_para=fit_args, verbose = verbose)
        
        if njoin_t is not None: 
            return (T - bls[:, None, :, :]).reshape((-1,) + shape_add)
        else:
            return T - bls
    if njoin_t is not None:
        num = len(T)//njoin_t*njoin_t
        if num != len(T):
            T = np.vstack([proc(T[:num]), proc(T[-njoin_t:])[num-len(T):]])
        else:
            T = proc(T)
    else:
        T = proc(T)
    res=T
    tempfile.save(rank, res=res)
    del(res, T)
    ##synchronize?
    aft_sub=True
    aft_sub= comm.gather(aft_sub, root=0)
    # save res
    if rank == 0:
        print("\nGathering" )
        sys.stdout.flush()
        res = [tempfile.load(i, 'res', rm=True) for i in range(size)]
        res = np.vstack(res)
        if trans:
            res[1] = res[1].transpose((1,0,2))
        print(f"Saving...")
        dict_out= {}
        dict_out['mjd'] = mjd_ori
        dict_out['ra'] = ra_ori
        dict_out['dec'] = dec_ori
        dict_out[outfield] = res.astype('float32')
        dict_out['freq'] = freq_ori
        if is_extrapo is not None:
            dict_out['is_extrapo'] = is_extrapo
        #save file
        from .util import add_extra
        add_extra(fs, dict_out)
        fs.close()
        from .util import rec_his, save_dict_hdf5
        header=rec_his(args=args)
        if header_in is not None: header.update(header_in)
        save_dict_hdf5(fileout, dict_out, header=header)
        print(f"Saved to {fileout}")
        