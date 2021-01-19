#!/usr/bin/env python
# coding: utf-8
import scipy.sparse as sparse
from scipy.sparse import linalg
from numpy.linalg import norm
import numpy as np
import scipy.interpolate as interp
from scipy import ndimage
from scipy import optimize
from .util import average_every_n, boxcar_smooth1d

import os
import warnings 
warnings.filterwarnings("ignore",r'overflow encountered in exp')

def get_baseline(x, ys, axis=None, *, s_method=None, s_sigma=None, average_every=None, exclude=None, verbose=False, method='arPLS', bl_para=None, return_f=False, check=True):
    """ 
    subtract the baseline of spectra with XX and YY polar.
    Parameters：
    ---------------
    x: array, shape (n,); freq or vel. 
    ys: array, T or Jy
    axis: axis along fitting baseline
    s_sigma, s_method: if s_sigma is not None, smoothing ys along axis using s_method ('gaussian' or 'boxcar')
    """
    if exclude is not None:
        if exclude.shape != ys.shape:
            raise(ValueError('exclude should has same shape with ys'))
    # check if ys is one dim      
    if ys.ndim == 1:
        ys = ys[:, None]
        if exclude is not None:
            exclude = exclude[:, None]
        ch = True
        axis = 0
    else:
        ch = False
    
    # smooth ys along axis
    if s_sigma is not None and s_method is not None:
        if s_method == 'gaussian':
            ys = ndimage.gaussian_filter1d(ys, s_sigma, axis=axis)
        elif s_method == 'boxcar':
            ys = boxcar_smooth1d(ys, s_sigma, axis=axis)
        elif s_method == 'median':
            ys = np.apply_along_axis(ndimage.median_filter, axis, ys, s_sigma)
        else:
            raise(ValueError('not support the method'))
    else:
        ys = ys
    
    if average_every is not None and average_every !=0:
        drop = False
        x_ori = x
        x = average_every_n(x, average_every, drop=drop)
        ys = average_every_n(ys, average_every, axis=axis, drop=drop)
        if exclude is not None:
            exclude = average_every_n(exclude, average_every, axis=axis, drop=drop).astype('bool')
    # after smooth(average), check nan, posinf, neginf
    if check:
        is_finite = np.isfinite(ys)
        if not is_finite.all():
            # replace value as 0 and add to exclude
            ys = np.copy(ys)
            ys[~is_finite] = 0.
            if exclude is not None:
                exclude = exclude | (~is_finite)
            else:
                exclude = ~is_finite
        del is_finite
    # select baseline method
    if method == 'arPLS':
        if 'sym' in bl_para.keys():
            bl_para = bl_para.copy() # shallow copy
            del bl_para['sym']
        BL = BL_arPLS(**bl_para)
        use_x = False
    elif method == 'srPLS':
        if 'sym' in bl_para.keys():
            bl_para = bl_para.copy() # shallow copy
            del bl_para['sym']
        BL = BL_srPLS(**bl_para)
        use_x = False
    elif method == 'Chebyshev':
        BL = BL_Chebyshev(**bl_para)
        use_x = True
    elif method == 'poly':
        BL = BL_poly(**bl_para)
        use_x = True
    elif method == 'sin_poly':
        BL = BL_sin_poly(**bl_para)
        use_x = True
    elif method == 'sin_poly_2':
        BL = BL_sin_poly_2(**bl_para)
        use_x = True
    else:
        raise(ValueError('method does not support'))
    
    # move axis to last
    ys = np.moveaxis(ys, axis, -1)
    if exclude is not None:
        exclude = np.moveaxis(exclude, axis, -1)
    shape_bak = ys.shape
    #ys = ys.reshape((-1,shape_bak[-1])) #will copy the data
    #or np.apply_along_axis
    bls = np.zeros_like(ys)
    if verbose:
        from tqdm import tqdm
        iter_ = tqdm(np.ndindex(*(shape_bak[:-1])), total=ys.size/shape_bak[-1], desc='CPU 0: ', mininterval=2)
    else:
        iter_= np.ndindex(*(shape_bak[:-1]))
 
    for ii in iter_:
        y = ys[ii+ np.s_[:,]]
        if exclude is not None:
            _exclude = exclude[ii+ np.s_[:,]]
        else:
            _exclude = None
        if use_x:
            bl = BL.fit(x, y, _exclude)
        else:
            bl = BL.fit(y, _exclude)
        bls[ii+ np.s_[:,]] = bl
    # move back
    bls = np.moveaxis(bls,-1,axis)
    
    if average_every is not None:
        bls = interp.interp1d(x, bls, axis=axis, kind='linear', fill_value ='extrapolate')(x_ori)
    if return_f:
        # move back
        ys = np.moveaxis(ys, -1, axis)
        if ch:
            return bls[:,0], ys[:,0], x
        else:
            return bls, ys, x
    else:
        if ch:
            return bls[:,0]
        else:
            return bls

def get_baseline_mp(n, x, ys, *, exclude=None,  **kwargs):
    """
    testing...
    only support axis=1
    """
    from functools import partial
    from multiprocessing import Process, Queue
    def get_baseline_q(q, *args, **kwargs):
        res = get_baseline(*args, **kwargs)
        q.put(res)
    ps = []
    qs = []
    if exclude is not None:
        exclude_list = np.array_split(exclude, n)
    else:
        exclude_list = [None,] * n
    ys_list = np.array_split(ys, n)
    for ys, exclude in zip(ys_list, exclude_list):
        q = Queue()
        p = Process(target=get_baseline_q, args=(q, x, ys), kwargs={'exclude':exclude, **kwargs})
        if 'verbose' in kwargs.keys():
            kwargs['verbose']=False

        ps+= [p] 
        qs+= [q] 
    for p in ps:
        p.start()
    res = np.vstack([q.get()for q in qs])
    for p in ps:
        p.join() # need after q.get()
    return np.stack(res)

    
class BL_arPLS(object):
    """
    (automatic) Baseline correction using asymmetrically reweighted penalized least squares smoothing. 
                     Baek et al. 2015, Analyst 140: 250-257
    Parameters
    ----------
    lam:
    offset:
    deg:
    ratio:
    niter:
    """
    def __init__(self, lam=10**12, offset=2, deg=3, ratio=0.01, niter=100, rew=True):
        self.lam = lam
        self.offset = offset
        self.deg = deg
        self.ratio = ratio
        self.niter = niter
        if not rew : self.niter = 1
    def _reweight(self, d):
        # make d- and get w^t with m and s
        dn = d[d<0]
        m = np.mean(dn)
        s = np.std(dn)
        wt = 1.0/(1 + np.exp( 2* (d-(self.offset*s-m))/s ))
        return wt
    
    def fit(self, y, exclude=None):
        # Adaptation of the Matlab code in Baek et al 2015 and python code in https://github.com/charlesll/rampy
        N = len(y)
        D= sparse.csc_matrix(sparse.eye(N))
        #D= self._diff(D)
        diff_fun= lambda x: x[:,1:]- x[:,:-1] # x is a csc sparse matrix
        for i in range(self.deg):
            D= diff_fun(D)
        w = np.ones(N)
        if exclude is not None:
            w[exclude] = 0
        for i in range(self.niter):
            W = sparse.spdiags(w, 0, N, N)
            Z = W + self.lam * D.dot(D.transpose())
            z = linalg.spsolve(Z, w*y)
            d = y - z
            wt = self._reweight(d)
            if exclude is not None:
                wt[exclude] = 0
            # check exit condition and backup
            ratio_fit= norm(w-wt)/norm(w)
            if ratio_fit < self.ratio:
                break
            w = wt

        bl = z
        self.success = True if ratio_fit <= self.ratio else False
        
        return bl

class BL_srPLS(BL_arPLS):
    def _reweight(self, d):
        dn = abs(d)
        m = np.mean(dn)/3
        s = np.std(dn)/3
        wt = 1.0/(np.exp( 2* (abs(d)-(self.offset*s-m))/s ))
        return wt
    

class BL_base(object):
    """
    (automatic) Baseline fit
    Parameters
    ----------
    offset:
    deg:
    ratio:
    niter:
    """
    def __init__(self, offset=2, deg=3, ratio=0.01, niter=100, rew=True, sym=False):
        self.offset = offset
        self.deg = deg
        self.ratio = ratio
        self.niter = niter
        self.rew = rew
        self._reweight = self._reweight_s if sym else self._reweight_a
           
    def _reweight_a(self, d):
        # make d- and get w^t with m and s
        dn = d[d<0]
        m = np.mean(dn)
        s = np.std(dn)
        wt = 1.0/(1 + np.exp( 2* (d-(self.offset*s-m))/s ))
        return wt
    
    def _reweight_s(self, d):
        dn = abs(d)
        m = np.mean(dn)/3
        s = np.std(dn)/3
        wt = 1.0/(np.exp( 2* (abs(d)-(self.offset*s-m))/s ))
        return wt
    
    def _fit(self, x, y, w):
        return y
    
    def fit(self, x, y, exclude=None):
        N = len(y)
        w = np.ones(N)
        if exclude is not None:
            w[exclude] = 0
        for i in range(self.niter):
            z = self._fit(x, y, w)
            if not self.rew:
                break
            d = y - z
            wt = self._reweight(d)
            if exclude is not None:
                wt[exclude] = 0
            # check exit condition and backup
            ratio_fit= norm(w-wt)/norm(w)
            if ratio_fit < self.ratio:
                break
            w = wt
        bl = z
        if self.rew:
            self.success = True if ratio_fit <= self.ratio else False
        return bl
    
class BL_Chebyshev(BL_base):
    """
    (automatic) Baseline correction using Chebyshev fit
    Parameters
    ----------
    offset:
    deg:
    ratio:
    niter:
    """
    def pred(self, x):
        return np.polynomial.chebyshev.chebval(x, self.para)
    
    def _fit(self, x, y, w):
        para = np.polynomial.chebyshev.chebfit(x, y, self.deg, w=w)
        self.para = para
        return self.pred(x)

class BL_poly(BL_base):
    """
    (automatic) Baseline correction using polynomial fit
    Parameters
    ----------
    offset:
    deg:
    ratio:
    niter:
    """
    def pred(self, x):
        return np.polynomial.polynomial.polyval(x, self.para)
    
    def _fit(self, x, y, w):
        para = np.polynomial.polynomial.polyfit(x, y, self.deg, w=w)
        self.para = para
        return self.pred(x)

class BL_sin_poly(BL_base):
    """
    (automatic) Baseline correction using sin plus poly
    Parameters
    ----------
    offset:
    deg:
    ratio:
    niter:
    """
    def __init__(self, f, ptype='poly', opt_para={}, **kwarg):
        super().__init__(**kwarg)
        if np.isscalar(f):
            f = [f]
        self.f = f
        self.ptype = ptype
        self.set_opt_para(**opt_para)
    def _fun(self, x, *arg):
        arg = np.asarray(arg)
        coef = arg[-(self.deg + 1):]
        y = np.zeros_like(x)
        for _A, _f, _p in zip(*np.split(arg[:-(self.deg + 1)],3)):
            y += _A*np.sin(2*np.pi*_f*x + _p)
        if self.ptype == 'poly':
            y += np.polynomial.polynomial.polyval(x, coef)
        elif 'cheb' in self.ptype.lower():
            y += np.polynomial.chebyshev.chebval(x, coef)
        return y
    
    def _err_func(self, para, x, y, w):
        residuals = y-self._fun(x, *para)
        return np.sum((residuals*w)**2, dtype='float64')
    
    def set_opt_para(self, method='Powell', **kwarg):
        """
        scipy.optimize.minimize parameter
        """
        self.opt_method = method
        self.opt_kwarg = kwarg
    
    def _fit(self, x, y, w):
        nsin = len(self.f)
        self.w=w
        std = np.std(y[w>=0.99])
        p0 = [std] + [std/2]*(nsin-1)
        p0 += self.f + [0,]*nsin
        p0 += [0]*(self.deg+1)
        
        #err_func = self._err_func(x, y, w, *arg)
        res = optimize.minimize(self._err_func, p0, args=(x,y,w), method=self.opt_method, **self.opt_kwarg)
        self.para = res.x
        self.fit_res = res
        return self._fun(x, *res.x)
    
    def pred(self, x):
        return self._fun(x, *self.para)

    
class BL_sin_poly_2(BL_sin_poly):
    """
    (automatic) Baseline correction using sin plus poly
    Parameters
    ----------
    offset:
    deg:
    ratio:
    niter:
    """
    
    def set_opt_para(self, **kwarg):
        """
        scipy.optimize.curve_fit parameter
        """
        self.opt_kwarg = kwarg
    
    def _fit(self, x, y, w):
        nsin = len(self.f)
        self.w = w
        std = np.std(y[w>=1])
        p0 = [std] + [std/2]*(nsin-1)
        p0 += self.f + [0,]*nsin
        p0 += [0]*(self.deg+1)
        
        res = optimize.curve_fit(self._fun, x, y, p0, sigma=1./w, **self.opt_kwarg)
        
        self.para = res[0]
        self.fit_res = res
        return self._fun(x, *self.para)
    