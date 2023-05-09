import os
import warnings

import numpy as np
from scipy import signal
from scipy import ndimage
from scipy.stats import binned_statistic

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

    
def convolve1d(arr, kernel1d, axis=-1, mode='reflect', method='auto'):
    '''
    One-dimensional convolving kernel1d along axis of arr

    Parameters
    ----------
    arr : array_like
        The input array.
    kernel1d : array_like
        dimension = 1
    axis : int, optional
        The axis of `input` along which to calculate. Default is -1.
    mode : str{'reflect', 'none'}, optional
        The mode parameter determines how "arr" is extended beyond its boundaries.
    method : str{'auto', 'direct', 'fft'} 
        see scipy.signal.convolve
    '''
    # boundary
    if mode == 'reflect':
        extend = int(len(kernel1d)/2)
        
        sli_h = [slice(None)]*arr.ndim
        sli_h[axis] = slice(None, extend)
        sli_h = tuple(sli_h)

        sli_t = [slice(None)]*arr.ndim
        sli_t[axis] = slice(-extend,None)
        sli_t = tuple(sli_t)

        sli_flip = [slice(None)]*arr.ndim
        sli_flip[axis] = slice(None, None, -1)
        sli_flip = tuple(sli_flip)

        if extend > arr.shape[axis]:
            n_add_ker = int(np.ceil(extend/arr.shape[axis]))//2*2+1  # odd number
            if n_add_ker < 3:
                n_add_ker = 3
            if n_add_ker > 9:
                raise(ValueError('kernel length is too larger that arr'))
            arr_ = np.concatenate([arr,] + [arr[sli_flip], arr,]*((n_add_ker-1)//2), axis=axis)
        else:
            arr_ = arr    
        tail = arr_[sli_t][sli_flip]
        head = arr_[sli_h][sli_flip]
        arr = np.concatenate((head, arr, tail),axis=axis)
        del(arr_, tail, head)
    elif mode == 'none':
        extend = None
    else:
        raise(ValueError("mode not support"))
    # change kernel shape    
    ker_shape = [1,]*arr.ndim
    ker_shape[axis] = len(kernel1d)
    kernel = kernel1d.reshape(ker_shape)
    # here,  np.sum(kernel) is same effect with np.sum(win_g,axis=axis) coz other axis shape is 1. 
    res = signal.convolve(arr, kernel, method=method, mode='same') / np.sum(kernel)
    if extend is not None:
        sli_res = [slice(None)]*arr.ndim
        sli_res[axis] = slice(extend,-extend)
        sli_res = tuple(sli_res)
        res = res[sli_res]
    return res

def smooth1d(arr, method, sigma, axis=-1):
    '''
    One-dimensional smooth arr along axis using different method

    Parameters
    ----------
    arr : array_like
        The input array.
    method : str{{'gaussian_fft', 'boxcar', 'gaussian', 'median'}}
        'gaussian_fft': use scipy.signal.convolve; sigma is std; extend boundaries using 'reflect'.
        'boxcar': use scipy.signal.convolve; average every int(2*sigma+1); extend boundaries using 'reflect'.
        'gaussian': use scipy.ndimage.gaussian_filter1d; sigma is std
        'median': use scipy.ndimage.median_filter; calcute the median value every int(2*sigma+1)
    axis : int, optional
        The axis of `input` along which to smooth. Default is -1.
    '''
    if 0 in arr.shape:
        warnings.warn("input array is empty")
        return arr
    if method == 'gaussian_fft':
        kernel1d = signal.windows.gaussian(2*(4*sigma)+1, sigma)
        return convolve1d(arr, kernel1d, axis=axis)
    elif method == 'boxcar':
        kernel1d = np.ones(int(2*sigma+1))
        return convolve1d(arr, kernel1d, axis=axis)
    elif method == 'gaussian':
        return ndimage.gaussian_filter1d(arr, sigma, axis=axis)
    elif method == 'median':
        return median_filter_1d(arr, size=int(2*sigma+1), axis=axis)
    else:
        raise(ValueError("method not supports"))

def smooth1d_fft(x, ys, Tcut, axis=-1):
    """
    x : array_like shape:(m,)
    ys : array_like
        The input array.
    Tcut : components with period smaller than Tcat will be removed.
    """
    if len(x) != ys.shape[axis]:
        raise ValueError(f"len(x) needs equal to and ys.shape[axis]")
    n = len(x)
    Fy = np.fft.rfft(ys, n=n, norm='ortho', axis=axis)
    f = np.fft.rfftfreq(len(x), x[1]-x[0])
    Fy = np.moveaxis(Fy, axis, -1)
    Fy[...,f > 1/Tcut] *= 0
    Fy = np.moveaxis(Fy, -1, axis)
    return np.fft.irfft(Fy, n=n, norm='ortho', axis=axis) # need specify n when ifft or irfft
