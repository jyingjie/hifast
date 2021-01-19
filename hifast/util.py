import os
import warnings

import numpy as np
from scipy import signal
from scipy import ndimage
from scipy.stats import binned_statistic
import h5py


def apply_along_axis(fun):
    #from functools import wraps
    #@wraps(fun)
    def wrapper(arr, axis, *args, **kwargs):
        """
    arr: array like
    axis: which axis to apply
    ----------------
    Others:"""
        return np.apply_along_axis(fun, axis, arr, *args, **kwargs)
    wrapper.__doc__ += fun.__doc__ 
    wrapper.__name__ = fun.__name__
    wrapper.__qualname__= fun.__qualname__
    return wrapper

@apply_along_axis
def mask_Trues(arr, leng_lim=20):
    """
    mask "continuous Trues" (number >= leng_lim) as True
    Parameters
    ----------
    leng_lim: int
    """
    if not np.any(arr):
        return arr
    arr = np.hstack([[False],arr,[False]])
    diff = np.diff(arr.astype('int32'))
    
    ind_neg = np.where(diff==-1)[0]
    ind_posi = np.where(diff==1)[0]

    leng = (ind_neg - ind_posi)
    is_use = leng < leng_lim 
    
    for i,j in zip(ind_posi[is_use],ind_neg[is_use]):
        arr[i+1:j+1]= False
    return arr[1:-1]

@apply_along_axis
def extend_Trues(is_rfi, leng_lim=1, ext_add=0, ext_frac=1./4):
    """
    extend "continuous Trues" to left and right
    
    Parameters
    ----------
    leng_lim: int; only extend "continuous Trues" which length >= leng_lim
    ext_add: int; fix number to extend
    frac: float; entend length*frac
    """
    if not np.any(is_rfi):
        return is_rfi
    is_rfi = np.hstack([[False],is_rfi,[False]])
    diff= np.diff(is_rfi.astype('int32'))
    
    ind_neg = np.where(diff==-1)[0]
    ind_posi= np.where(diff==1)[0]
#     if ind_neg[0]==0:
#         ind_neg= ind_neg[1:]
#     if ind_posi[-1] == len(diff)-1:
#         ind_posi = ind_posi[:-1]

    leng = (ind_neg- ind_posi)
    is_use= leng>= leng_lim # not needed
    l_ext= np.ceil(leng[is_use]*ext_frac).astype('int') + ext_add
    #l_ext= np.full(len(leng),5)
    for i,n in zip(ind_neg[is_use],l_ext):
        is_rfi[i+1:i+n+1]= True
    for i,n in zip(ind_posi[is_use],l_ext):
        is_rfi[max(0,i-n+1):i+1]= True
    return is_rfi[1:-1]

@apply_along_axis
def median_filter_1d(arr, *args, **kwargs):
    """
    see scipy.ndimage.median_filter
    """
    return ndimage.median_filter(arr, *args, **kwargs)

def median_filter_axis1_d3(arr, kernel_size=5):
    return median_filter_1d(arr, axis=1, size=kernel_size)


def average_every_n(arr, n, axis=-1, drop=True):
    """
    drop: if drop arr.shape[axis]%n in end
    """
    len_ = arr.shape[axis]
    num = len_//n*n
    drop = num==len_ or drop
        
    def _fun(arr_):
        res= np.mean(arr_[:num].reshape((-1,n)), axis=1, dtype='float64')
        if not drop:
            res2 = np.mean(arr_[-n:].reshape((-1,n)), axis=1, dtype='float64')
            res=np.hstack([res,res2])
        return res
    return np.apply_along_axis(_fun, axis=axis, arr=arr)

def down_sample(data,dfactor):
    """
    average every n values along axis 1, n equal dfactor
    Parameters
    ----------
    data : array_like
        shape is (x,N,y)

    axis : int, optional
        The axis in the result array along which the input arrays are stacked.

    dfactor : int
        number of values that take the average

    Returns
    -------
    data_d : ndarray
        shape is (x,N//dfactor,y)

    """
    length=data.shape[1]
    bins=np.arange(0,length+1,dfactor)
    data_d= [binned_statistic(np.arange(length)+0.5,data[:,:,i],'mean',bins=bins)[0] for i in range(data.shape[2])]
    data_d= np.stack(data_d,axis=2)
    return data_d

def save_dict_hdf5(fname, dict_in, header=None, mode='w'):
    """
    fname: str
    dict_in: dict; keys of the dict_in are str, value are numpy array like.
    header: dict
    """
    try:
        f= h5py.File(fname,mode)
    except OSError:
        if mode=='w':
            from datetime import datetime
            os.rename(fname,fname+'.bak.empty.'+datetime.now().strftime("%Y%m%d-%H%M%S"))
            f= h5py.File(fname,mode)
    if header is not None:
        # sometimes hdf5 raises error if track_order = True
        f.create_group('Header',track_order=False)
        for key in header.keys():
            f['Header'].attrs[key]= header[key]
    for key in dict_in.keys():
        f[key]= dict_in[key]
    f.close()
    
def freq2vlsr(freq, ra, dec, mjd):
    '''
    HI freq to velocity under LSR
    ------------
    ra, dec: deg; array or scalar
       len(ra)==len(mjd)
    '''
    from .ugdopplerfast import ugdopplerfast
    mjd= np.asarray([mjd]) if np.isscalar(mjd) else np.asarray(mjd)
    ra= np.asarray([ra]) if np.isscalar(ra) else np.asarray(ra)
    dec= np.asarray([dec]) if np.isscalar(dec) else np.asarray(dec)
    
    restfreq = np.array([1420.405751])
    c        =  2.99792458e5  # km/s
    velo     =  c*(restfreq[0]- freq)/ restfreq[0] 
    
    jd = mjd+2400000.5

    vlsrcor  = ugdopplerfast(ra,dec,jd)
    vlsr  = np.vstack([velo-ii for ii in vlsrcor])
    if len(vlsrcor)==1:
        vlsr= vlsr[0]
    return vlsr
    

def gaussian_smooth1d(vals, sigma, axis=0):
    '''
    One-dimensional Gaussian smooth.

    Parameters
    ----------
    vals : array_like
        The input array.
    sigma : scalar
        standard deviation for Gaussian kernel
    axis : int, optional
        The axis of `input` along which to calculate. Default is 0.
    '''
    N= min(vals.shape[axis], int(6*sigma))
    #N= vals.shape[axis]
    win_shape= [1 if i!=axis else N  for i in range(len(vals.shape))]
    win_g= signal.windows.gaussian(win_shape[axis],sigma).reshape(win_shape)
    # here,  np.sum(win_g) is same effect with np.sum(win_g,axis=axis) coz other axis shape is 1. 
    res= signal.convolve(vals, win_g, method='fft',mode='same') / np.sum(win_g)      
    return res

def boxcar_smooth1d(vals, sigma, axis=0):
    '''
    One-dimensional moving mean (boxcar) smooth.

    Parameters
    ----------
    vals : array_like
        The input array.
    sigma : scalar
        Number of points = 2*sigma+1
    axis : int, optional
        The axis of `input` along which to calculate. Default is 0.
    '''
    
    win_shape= [1 if i!=axis else 2*int(sigma)+1  for i in range(len(vals.shape))]
    win_g= np.ones(win_shape[axis]).reshape(win_shape)
    #print(win_g)
    res= signal.convolve(vals, win_g, method='fft',mode='same') / np.sum(win_g)
    return res

def smooth_axis1_d3(vals, x=None, method=None, sigma=None, deg=None):
    """
    s_para: dict
    """
    if len(vals) == 0:
        warnings.warn("input array is empty in smooth_axis1_d3")
        return vals
    if vals.ndim !=3:
        raise(ValueError('input array should 3 dim'))
    if method=='poly':
        res= np.zeros_like(vals)
        for i in range(vals.shape[0]):
            for j in range(vals.shape[2]):
                res[i,:,j]= np.poly1d(np.polyfit(x, vals[i,:,j], deg=deg))(x)
        return res
    
    if method=='gaussian':
        extend= min((int(sigma*4), vals.shape[1]))

    #     if vals.shape[1] < extend:
    #         num= int(np.ceil(extend/vals.shape[1]))
    #         tail= np.concatenate([vals,]*1,axis=1)[:,-extend:,:][:,::-1,:]
    #         head= np.concatenate([vals,]*1,axis=1)[:,:extend,:][:,::-1,:]
    #     else:
        tail= vals[:,-extend:,:][:,::-1,:]
        head= vals[:,:extend,:][:,::-1,:]

        vals_new= np.concatenate((head, vals, tail),axis=1)
        res= gaussian_smooth1d(vals_new,sigma,axis=1)
        res= res[:,extend:-extend,:]

        return res

    
def read_tcal_sav(nB, s_type='w', tcal_dir=None, mode='high', date='20190115'):
    """
    nB: int
      beam number
    """
    from scipy.io.idl import readsav
    if tcal_dir is None:
        tcal_dir= os.path.expanduser("~")+'/Tcal/'
    fname = tcal_dir+f'{date}/median_{date}.Tcal-results.HI_{s_type}.{mode}.sav'
    tc_info= readsav(fname)[f'{mode}_{s_type}'][0]
    tc_freq= tc_info['freq']
    return tc_freq, tc_info['M%02d_TC'%nB], fname

def read_tcal_fits(nB, s_type='w', tcal_dir=None, mode='high', date=''):
    """
    nB: int
      beam number
    """
    from astropy.io import fits
    import os
    if tcal_dir is None:
        tcal_dir= os.path.expanduser("~")+'/Tcal/'
    fname = tcal_dir + f'{date}/CAL.{date}.{mode}.{s_type.upper()}.fits'
    f = fits.open(fname)
    tc_freq = f[1].data['FREQ'][0]
    tc_T = f[1].data['TCAL'][0,nB-1].T
    return tc_freq, tc_T, fname

def read_tcal(nB, s_type='w', tcal_dir=None, mode='high', date='auto', mjd=None):
    """
    nB: int
       beam number
    s_type: str
       type, w or n
    tcal_dir: str
       if not set, use ~/Tcal/
    mode: str
       high or low
    date: str
       example: 20190115 or 20200531, default is None
    mjd:
       if date is 'auto', using the nearest date of tcal
    """
    from glob import glob
    import os
    if tcal_dir is None:
        tcal_dir= os.path.expanduser("~")+'/Tcal/'
    dates_have = [os.path.basename(i) for i in glob(tcal_dir + '/20[0-9][0-9][0-9][0-9][0-9][0-9]')]
    if len(dates_have) == 0:
        raise(ValueError('can not find tcal file'))
    if date == 'auto':
        from astropy.time import Time
        if mjd is None:
            raise(ValueError('need input mjd if date is auto'))
        mjds_h = Time([f'{s[:4]}-{s[4:6]}-{s[6:8]} 00:00:00.000' for s in dates_have], format='iso').mjd
        #print(mjds_h)
        date = dates_have[np.argmin(abs(mjds_h - mjd))]
    else:
        if date not in dates_have:
            raise(ValueError(f'can not find tcal file in {date}'))
    ## read tcal
    if date == '20190115':
        return read_tcal_sav(nB, s_type, tcal_dir, mode, date='20190115')
    else:
        return read_tcal_fits(nB, s_type, tcal_dir, mode, date=date)
    
def rec_his(**kwargs):
    import sys
    import json
    from datetime import datetime
    from collections import OrderedDict
    from ._version import get_versions
    
    history= OrderedDict()
    history['version']= get_versions()['version']
    history['cwd']= os.getcwd()
    history['argv']= ' '.join(sys.argv)
    for key in kwargs.keys():
        history[key]= repr(kwargs[key])
    current_time= datetime.now().strftime("%Y%m%d-%H:%M:%S")
    return OrderedDict({"HISTORY-"+current_time : json.dumps(history,indent=2)})

def add_extra(fin, _dict=None, fields_add=[]):
    """
    add some fields of fin to _dict
    """
    fields = ['is_on', 'next_to_cal', 'is_delay', 'Tcal', 'is_extrapo', 'vel']
    fields += fields_add
    out_add = {}
    for field in fields:
        if field in fin.keys():
            out_add[field] = fin[field][:]
    if _dict is None:
        return out_add
    else:
        _dict.update(out_add)
