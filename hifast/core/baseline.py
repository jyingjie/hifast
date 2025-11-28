

__all__ = ['get_baseline', 'get_baseline_mp', 'BL_base', 'BL_PLS', 'BL_arPLS', 'BL_Gauss', 'BL_Chebyshev', 'BL_poly',
           'BL_sin_poly', 'BL_sin_poly_2', 'piecewise_polyfit', 'BL_knpoly', 'BL_spline', 'BL_knspline', 'BL_asPLS',
           'BL_masPLS', 'get_exclude_fun', 'sub_baseline']


from functools import lru_cache
import sys

import scipy.sparse as sparse
from scipy.sparse import linalg
from numpy.linalg import norm
import numpy as np
import scipy.interpolate as interp
from scipy import ndimage
from scipy import optimize
from threadpoolctl import ThreadpoolController

from astropy.stats import sigma_clipped_stats

from ..utils.misc import average_every_n, smooth1d, smooth1d_fft, extend_Trues

import os
import warnings
warnings.filterwarnings("ignore", r'overflow encountered in exp')
warnings.filterwarnings("ignore", r'The fit may be poorly conditioned')

try:
    import bottleneck as bn
    MEDIAN = bn.median
    MEAN = bn.nanmean
    STD = bn.nanstd
except ImportError:
    MEDIAN = np.median
    MEAN = np.mean
    STD = np.std


def get_baseline(x, ys, axis=None, *,
                 s_method='none', s_sigma=3, average_every=None, exclude=None,
                 exclude_add='none',
                 method='arPLS', bl_para=None,
                 verbose=False, return_f=False, check=True,
                 interp_nan=False,
                 ):
    """
    subtract the baseline of spectra with XX and YY polar.
    Parameters：
    ---------------
    x: array, shape (n,); freq or vel.
    ys: array, T or Jy
    axis: axis along fitting baseline
    s_sigma, s_method: if s_sigma is not 'none', smoothing ys along axis using s_method ('gaussian' or 'boxcar')
    average_every: int; if set as n and n>1, then average ys along axis every n points
    return_f: If Fasle, return baseline, else return (basline, ys_processed, x_processed, weights)
    check: if True, check if there are nan or inf in ys, if there are, add them to exclude if intep_nan is False.
    intep_nan: if True and check is True, replace nan values with interp.interp1d but not add to exclude.
    """
    if x is not None:
        x = np.array(x, dtype='float')
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
    if s_method is not None and s_method != 'none':
        if s_method in ['boxcar', 'gaussian', 'median']:
            ys = smooth1d(ys, s_method, s_sigma, axis=axis)
        elif s_method == 'fft':
            ys = smooth1d_fft(x, ys, s_sigma, axis=axis)
        else:
            raise(ValueError('not support the method'))
    else:
        ys = ys

    if average_every is not None and average_every > 1:
        drop = False
        x_ori = x
        x = average_every_n(x, average_every, drop=drop)
        ys = average_every_n(ys, average_every, axis=axis, drop=drop)
        if exclude is not None:
            exclude = average_every_n(
                exclude, average_every, axis=axis, drop=drop).astype('bool')

    if 'PLS' in method:
        bl_para = bl_para.copy()  # shallow copy
        if 'sym' in bl_para.keys():
            del bl_para['sym']
        if bl_para['deg'] > 3:
            warnings.warn('deg is too large in PLS methods and has been set to 3.')
            bl_para['deg'] = 3  # bl_para has been shallow copyed
    if 'spline' in method:
        #  also for 'knspline'
        bl_para = bl_para.copy() # shallow copy
        if 'sym' in bl_para.keys():
            del bl_para['sym']
        if bl_para['deg'] > 5:
            warnings.warn('deg is too large in spline methods and has been set to 5.')
            bl_para['deg'] = 5  # bl_para has been shallow copyed
        if method.startswith('spline'):
            if 'knots' in bl_para.keys():
                del bl_para['knots']
        if method.startswith('knspline'):
            if 'knots' not in bl_para.keys() or bl_para['knots'] is None:
                raise(ValueError('need knots for knspline- method'))
    if method.startswith('knpoly'):
        if 'knots' not in bl_para.keys() or bl_para['knots'] is None:
            raise(ValueError('need knots for knpoly- method'))
    # select baseline method

    if '-' in method:
        p1, p2 = method.split('-')
        if p1 == 'PLS':
            BL = BL_PLS(**bl_para, rew_type=p2)
        elif p1 == 'spline':
            BL = BL_spline(**bl_para, rew_type=p2)
        elif p1 == 'knspline':
            BL = BL_knspline(**bl_para, rew_type=p2)
        elif p1 == 'Chebyshev':
            BL = BL_Chebyshev(**bl_para, rew_type=p2)
        elif p1 == 'poly':
            BL = BL_poly(**bl_para, rew_type=p2)
        elif p1 == 'knpoly':
            BL = BL_knpoly(**bl_para, rew_type=p2)
        elif p1 == 'Gauss':
            BL = BL_Gauss(**bl_para, rew_type=p2)
        elif p1 == 'sin_poly':
            BL = BL_sin_poly(**bl_para, rew_type=p2)
        elif p1 == 'masPLS':
            BL = BL_masPLS(**bl_para, rew_type=p2)
        else:
            raise(ValueError('not supported method: {method}'))
    elif method == 'arPLS':
        BL = BL_PLS(**bl_para, rew_type='asym')
    elif method == 'srPLS':
        BL = BL_PLS(**bl_para, rew_type='sym')
    elif method == 'masPLS':
        BL = BL_masPLS(**bl_para, rew_type='asym')
    elif method == 'asPLS':
        BL = BL_asPLS(**bl_para)
    elif method == 'sin_poly':
        BL = BL_sin_poly(**bl_para, rew_type='asym')
    elif method == 'original':
        return ys
    else:
        raise(ValueError('not supported method: {method}'))


    # move axis to last
    ys = np.moveaxis(ys, axis, -1)
    if exclude is not None:
        exclude = np.moveaxis(exclude, axis, -1)
    shape_bak = ys.shape
    # ys = ys.reshape((-1,shape_bak[-1])) #will copy the data
    # or np.apply_along_axis
    bls = np.zeros_like(ys)
    if return_f:
        weis = np.zeros_like(ys, dtype=np.float64)
    if verbose:
        from tqdm import tqdm
        iter_ = tqdm(np.ndindex(
            *(shape_bak[:-1])), total=ys.size/shape_bak[-1], desc='CPU 0: ', mininterval=2)
    else:
        iter_ = np.ndindex(*(shape_bak[:-1]))

    for ii in iter_:
        y = ys[ii + np.s_[:, ]]
        if exclude is not None:
            _exclude = exclude[ii + np.s_[:, ]]
        else:
            _exclude = None
        if check:
            is_finite = np.isfinite(y)
            if not is_finite.all():
                y_finite = y[is_finite]
                # if all nan
                if 0 in y_finite.shape:
                    bls[ii + np.s_[:, ]] = np.nan
                    continue
                if not interp_nan:
                    # add nan to exclude and replace nan with max or any value
                    y = np.copy(y)
                    y[~is_finite] = np.max(y_finite)
                    if exclude is not None:
                        _exclude = _exclude | (~is_finite)
                    else:
                        _exclude = ~is_finite
                else:
                    # not add to exclude but repalce nan values with interp.interp1d
                    y = np.copy(y)
                    y[~is_finite] = interp.interp1d(
                        x[is_finite], y_finite, kind='linear', fill_value='extrapolate')(x[~is_finite])

            if _exclude is not None and _exclude.all():
                warnings.warn('All points in this spectrum are excluded, it\'s baseline will set as all zeros.')
                bls[ii + np.s_[:, ]] = 0
                continue
        try:
            bl = BL.fit(x=x, y=y, exclude=_exclude, exclude_add=exclude_add)
        except Exception as err:
            print('fit baseline fail, return zero array. Try other method, or there are bugs and ...')
            print('Errors output:\n')
            print(err)
            print('')
            sys.stdout.flush()
            continue # bls initialized as zeros
        bls[ii + np.s_[:, ]] = bl
        if return_f:
            weis[ii + np.s_[:, ]] = BL.wei
    # move back
    bls = np.moveaxis(bls, -1, axis)
    if return_f:
        weis = np.moveaxis(weis, -1, axis)
    if average_every is not None and average_every > 1:
        bls = interp.interp1d(x, bls, axis=axis, kind='linear',
                              fill_value='extrapolate')(x_ori)
    if return_f:
        # move back
        ys = np.moveaxis(ys, -1, axis)
        if ch:
            return bls[:, 0], ys[:, 0], x, weis[:, 0]
        else:
            return bls, ys, x, weis
    else:
        if ch:
            return bls[:, 0]
        else:
            return bls


controller = ThreadpoolController()



class get_baseline_mp(object):
    def __init__(self, n):
        # n process
        self.n = n

    @staticmethod
    def get_baseline_q(q, *args, **kwargs):
        res = get_baseline(*args, **kwargs)
        q.put(res)

    @controller.wrap(limits=1, user_api='blas')
    def __call__(self, x, ys, *, exclude=None,  **kwargs):
        """
        testing...
        only support axis=1,
        not support return_f=True
        """
        from functools import partial
        from multiprocessing import Process, Queue

        if 'axis' in kwargs.keys():
            if kwargs['axis'] != 1:
                raise(ValueError('get_baseline_mp only support axis=1 for now'))
        else:
            kwargs['axis'] = 1
        n = min(self.n, len(ys))

        ps = []
        qs = []
        if exclude is not None:
            exclude_list = np.array_split(exclude, n)
        else:
            exclude_list = [None, ] * n
        ys_list = np.array_split(ys, n)
        for ys, exclude in zip(ys_list, exclude_list):
            q = Queue()
            p = Process(target=self.get_baseline_q, args=(q, x, ys),
                        kwargs={'exclude': exclude, **kwargs})
            if 'verbose' in kwargs.keys():
                kwargs['verbose'] = False

            ps += [p]
            qs += [q]
        for p in ps:
            p.start()
        res = np.vstack([q.get()for q in qs])
        for p in ps:
            p.join()  # need after q.get()
        return np.stack(res)


_rew_fun_dict = {}

#@classmethod
def _reweight_asym1(self, d):
    """
    same to ``Baek et al. 2015, Analyst 140: 250-257``
    """
    if self.exclude is None:
        dn = d[d < 0]
    else:
        dn = d[(d < 0) & (~self.exclude)]
    m = MEAN(dn)
    s = STD(dn)
    wt = 1./(1 + np.exp(2 * (d-(self.offset*s-m))/s))
    return wt

_rew_fun_dict['asym1'] = _reweight_asym1

def _reweight_asym2(self, d):
    """
    similar to ``Baek et al. 2015, Analyst 140: 250-257``
    """
    if self.exclude is None:
        dn = d[d < 0]
    else:
        dn = d[(d < 0) & (~self.exclude)]
    m = MEAN(dn)
    s = STD(dn)
    wt = 1./(1 + np.exp(2 * (abs(d)-(self.offset*s-m))/s))
    return wt

_rew_fun_dict['asym2'] = _reweight_asym2

def _reweight_asym3(self, d):
    """
    similar to ``Baek et al. 2015, Analyst 140: 250-257``
    """
    if self.exclude is None:
        dn = d[d < 0]
    else:
        dn = d[(d < 0) & (~self.exclude)]
    m = MEAN(dn)
    s = STD(dn)
    wt = np.zeros_like(d)
    wt[abs(d) < self.offset*s-m] = 1
    return wt

_rew_fun_dict['asym3'] = _reweight_asym3

def _reweight_sym1(self, d):
    """
    """
    dn = d if self.exclude is None else d[~self.exclude]
    mean, median, stddev = sigma_clipped_stats(dn)
    m = median
    s = stddev/2.
    wt = 1./(1 + np.exp(2 * (abs(d)-(self.offset*s-m))/s))
    return wt

_rew_fun_dict['sym1'] = _reweight_sym1

## set default
_rew_fun_dict['asym'] = _reweight_asym1
_rew_fun_dict['sym'] = _reweight_sym1


class BL_base(object):
    """
    (automatic) Baseline fit
    Parameters
    ----------
    offset:
    deg:
    rew_type:
    ratio:
    niter:
    rew: deprecated, use niter=1 instead
    """

    def __init__(self, *, offset=2, deg=3, ratio=0.01, niter=100, rew_type='asym', rew=True, sym=False):

        self.offset = offset
        self.deg = deg
        self.ratio = ratio
        self.niter = niter
        BL_base._reweight = _rew_fun_dict[rew_type]
        if not rew:
            self.niter = 1

    def _fit(self, x, y, w):
        ...
        return y

    def fit(self, *, x=None, y=None, exclude=None, exclude_add='none', wei=None):
        """
        iteration

        Parameters
        ----------
        x: None or array (N,)
        y: array (N,); nan value need be replaced as normal nubmer in advance, and put them in ``exclude``
        exclude: bool array (N,); masked points that will never be used in fitting baseline
        exclude_add: str; additional exclude at last step, 'none', 'auto', 'auto1', 'auto2'
        wei: array; initial weigth for each point

        """
        # Adaptation of the code in Baek et al 2015 and https://github.com/charlesll/rampy
        self.exclude = exclude
        N = len(y)
        wt = np.ones(N) if wei is None else np.copy(wei)
        if self.exclude is not None:
            wt[self.exclude] = 0
        for i in range(self.niter):
            w = wt
            z = self._fit(x, y, w)
            d = y - z
            wt = self._reweight(d)
            # also used to deal replaced nan value, need applied always
            if self.exclude is not None:
                wt[self.exclude] = 0
            # check exit condition and backup
            ratio_fit = norm(w-wt)/norm(w)
            if ratio_fit < self.ratio:
                break
        exclude2 = None
        if exclude_add == 'auto' or exclude_add == 'auto1':
            is_ = w<0.01
            if self.exclude is not None:
                is_[self.exclude] = False # not extend the input exclude range
            exclude2 = extend_Trues(is_, axis=-1, ext_frac=1/3)
        elif exclude_add == 'auto2':
            ds = ndimage.gaussian_filter1d(d, 3)
            is_ = ds > 3*STD(d[d<0])
            if self.exclude is not None:
                is_[self.exclude] = False # not extend the input exclude range
            exclude2 = extend_Trues(is_, axis=-1, ext_frac=1/3)
        if exclude2 is not None:
            if self.exclude is not None:
                exclude2 |= self.exclude
            w = wt # use wt
            w[exclude2] = 0
            z = self._fit(x, y, w)

        bl = z
        self.success = True if ratio_fit <= self.ratio else False
        self.wei = w
        self.i = i
        return bl


class BL_PLS(BL_base):
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

    def __init__(self, lam=10**8, **kwargs):
        """
        lam:
        offset:
        deg:
        ratio:
        niter:
        """
        self.lam = lam
        super().__init__(**kwargs)

    @lru_cache(maxsize=20)
    def get_H(self, N):
        """
        """
        D = sparse.csc_matrix(sparse.eye(N))
        for i in range(self.deg):
            D = D[:, 1:] - D[:, :-1] # D is a csc sparse matrix, don't use np.diff
        H = D.dot(D.transpose()).multiply(self.lam)
        return H

    def _fit(self, x, y, w):
        """
        """
        N = len(y)
#         H = getattr(self, 'H', self.get_H(N))
#         self.H = H
        H = self.get_H(N)
        W = sparse.spdiags(w, 0, N, N)
        Z_m = W + H
        z = linalg.spsolve(Z_m, w*y)
        return z

BL_arPLS = BL_PLS


class BL_Gauss(BL_base):
    """
    """
    def __init__(self, lam=200, **kwargs):
        """
        lam: here as gaussian sigma
        """
        self.lam = int(lam)
        super().__init__(**kwargs)

    def _fit(self, x, y, w):
        z = smooth1d(y*w, 'gaussian_fft', self.lam, axis=-1)
        z /= smooth1d(w, 'gaussian_fft', self.lam, axis=-1)
        return z



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
        for _A, _f, _p in zip(*np.split(arg[:-(self.deg + 1)], 3)):
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
        self.w = w
        std = np.std(y[w >= 0.99])
        p0 = [std] + [std/2]*(nsin-1)
        p0 += self.f + [0, ]*nsin
        p0 += [0]*(self.deg+1)

        #err_func = self._err_func(x, y, w, *arg)
        res = optimize.minimize(self._err_func, p0, args=(
            x, y, w), method=self.opt_method, **self.opt_kwarg)
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
        std = np.std(y[w >= 1])
        p0 = [std] + [std/2]*(nsin-1)
        p0 += self.f + [0, ]*nsin
        p0 += [0]*(self.deg+1)

        res = optimize.curve_fit(
            self._fun, x, y, p0, sigma=1./w, **self.opt_kwarg)

        self.para = res[0]
        self.fit_res = res
        return self._fun(x, *self.para)



def piecewise_polyfit(x, y, weights, knots, deg, xnew):
    # Find indices at which to split data based on knots
    split_indices = np.append(np.searchsorted(x, knots), len(x))
    split_indices_new = np.append(np.searchsorted(xnew, knots), len(xnew))
    # Fit a polynomial to each section
    ynew = np.zeros(len(xnew), dtype=y.dtype)
    start = 0
    start_new = 0
    for split_index, split_index_new in zip(split_indices, split_indices_new):
        # Slice data for this section
        xs = x[start:split_index]
        ys = y[start:split_index]
        ws = weights[start:split_index]

        # Fit polynomial to this section
        p = np.polynomial.polynomial.polyfit(xs, ys, deg=deg, w=ws)

        # Evaluate polynomial for this section
        ynew[start_new:split_index_new] = np.polynomial.polynomial.polyval(xnew[start_new:split_index_new], p)

        # Update start indices for next section
        start = split_index
        start_new = split_index_new

    return ynew


class BL_knpoly(BL_base):
    """
    (automatic) Baseline correction using
    Parameters
    ----------
    offset:
    deg:
    ratio:
    niter:
    """
    def __init__(self, knots=None, **kwargs):
        """
        """
        self.knots = knots
        super().__init__(**kwargs)

    def _fit(self, x, y, w):
        return piecewise_polyfit(x, y, weights=w, knots=self.knots, deg=self.deg, xnew=x)


class BL_spline(BL_base):
    """
    (automatic) Baseline correction using
    Parameters
    ----------
    offset:
    deg:
    ratio:
    niter:
    """
    def __init__(self, lam=10**8, **kwargs):
        """
        """
        self.lam = lam
        super().__init__(**kwargs)

    def _fit(self, x, y, w):
        spl = interp.UnivariateSpline(x, y, w, k=self.deg, s=self.lam)
        return spl(x)


class BL_knspline(BL_base):
    """
    (automatic) Baseline correction using
    Parameters
    ----------
    offset:
    deg:
    ratio:
    niter:
    """
    def __init__(self, knots=None, **kwargs):
        """
        """
        self.knots = knots
        super().__init__(**kwargs)

    def _fit(self, x, y, w):
        spl = interp.LSQUnivariateSpline(x, y, t=self.knots, w=w, k=self.deg)
        return spl(x)



def _reweight_asPLS(self, d):
        # make d- and get w^t with m and s
        dn = d[d < 0]
        s = STD(dn, ddof=1)
        wt = 1.0/(1 + np.exp(self.offset * (d-s)/s))
        return wt

class BL_asPLS(BL_arPLS):
    """
    (automatic) Baseline correction using a modified version of adaptive smoothness parameter penalized least squares method.
        Feng Zhang et al. 2020, DOI: 10.1080/00387010.2020.1730908

    Parameters
    ----------
    lam:
    offset:
    deg:
    ratio:
    niter:
    """
    def __init__(self, *args, **kwrags):
        super().__init__(*args, **kwrags)
        BL_asPLS._reweight = _reweight_asPLS

    def _update_alpha(self, d):
        d_abs = np.abs(d)
        return d_abs/np.max(d_abs)

    def fit(self, *, x=None, y=None, exclude=None, wei=None):
        self.exclude = exclude
        N = len(y)
        D = sparse.csc_matrix(sparse.eye(N))
        #D= self._diff(D)

        def diff_fun(x): return x[:, 1:] - \
            x[:, :-1]  # x is a csc sparse matrix
        for i in range(self.deg):
            D = diff_fun(D)
        w = np.ones(N) if wei is None else np.copy(wei)
        alpha = np.ones(N)
        if self.exclude is not None:
            w[self.exclude] = 0
        H = D.dot(D.transpose()).multiply(self.lam)
        for i in range(self.niter):
            W = sparse.spdiags(w, 0, N, N)
            Z = W + H.multiply(alpha[:, None])
            z = linalg.spsolve(Z, w*y)
            d = y - z
            wt = self._reweight(d)
            if self.exclude is not None:
                wt[self.exclude] = 0
            # check exit condition and backup
            ratio_fit = norm(w-wt)/norm(w)
            if ratio_fit < self.ratio:
                break
            w = wt
            alpha = self._update_alpha(d)

        bl = z
        self.success = True if ratio_fit <= self.ratio else False
        self.wei = w
        return bl


class BL_masPLS(BL_asPLS, BL_arPLS):
    def __init__(self, *args, **kwrags):
        BL_arPLS.__init__(self, *args, **kwrags)
    def _update_alpha(self, d):
        d_abs = np.abs(d)
        return 1+d_abs/np.max(d_abs)


def get_exclude_fun(exclude_m):
    if exclude_m == 0:
        def exclude_fun(x):
            return extend_Trues(abs(x) > 1.*np.diff(np.percentile(x, [16, 84], axis=1), axis=0)[0][:, None, :], axis=1, ext_frac=1/3, ext_add=3)
    elif exclude_m == 1:
        def exclude_fun(x):
            return extend_Trues(abs(x) > 2.5*np.min(np.diff(np.percentile(x, [10, 50, 90], axis=1), axis=0), axis=0)[:, None, :], axis=1, ext_frac=1/3, ext_add=3)
    return exclude_fun


def sub_baseline(freq, yss, *, subtract=True, nproc=1, exclude_fun=None, is_excluded=None, exclude_add='none', inplace=False,
                 njoin=1, s_method_t='none', s_sigma_t=None,
                 method='arPLS', s_method_freq='none', s_sigma_freq=None, average_every_freq=None,
                 lam=1e8, deg=2, offset=2, ratio=0.01, niter=100, sin_f=[0.925, ], rew=True, opt_para=None,
                 verbose=True,
                 knots=None,
                 ):
    """
    freq: array, shape:(n,) ; Frequency.
    yss: array, shape:(m,n,2)

    njoin: int
    nproc:
    exclude_fun:
    s_method_t: smooth method along axis 0
    s_sigma_t:
    s_method_freq:
    s_sigma_freq:
    average_every_freq:
    method: str; baseline fitting method
    bl_para: baseline fitting parameter
    verbose:
    """
    # check
    if yss.ndim != 3:
        raise(ValueError('yss should has 3-dim'))
    if yss.shape[1] != len(freq):
        raise(ValueError('the 2nd dimension of yss should equal to the length of freq'))

    # baseline fit parameter
    bl_para = {
        'deg': deg,
        'offset': offset,
        'ratio': ratio,
        'niter': niter,
        'rew': rew,
    }
    if 'PLS' in method or 'Gauss' in method:
        bl_para['lam'] = lam
    if method.startswith('spline'):
        bl_para['lam'] = lam
    if method.startswith('knspline') or method.startswith('knpoly'):
        if knots is None:
            raise(ValueError('*-knspline methods need knots input'))
        bl_para['knots'] = knots
    if 'sin' in method:
        bl_para['f'] = sin_f
        if opt_para is not None:
            bl_para['opt_para'] = opt_para

    yss_ori = yss
    # smooth spectra along t
    if s_method_t is not None and s_method_t != 'none':
        if s_method_t in ['boxcar', 'gaussian', 'median']:
            yss = smooth1d(yss, s_method_t, s_sigma_t, axis=0)
        elif s_method_t == 'fft':
            yss = smooth1d_fft(np.arange(yss.shape[0]), yss, s_sigma_t, axis=0)
        else:
            raise(ValueError('not support the method'))
    # join
    if njoin is not None and njoin > 1:
        yss = average_every_n(yss, njoin, axis=0, drop=False)
    # add exclude
    if is_excluded is not None:
        exclude = is_excluded
        if njoin is not None and njoin > 1:
            # todo: using np.logical_and.reduce
            exclude = average_every_n(exclude, njoin, axis=0, drop=False).astype('bool')
    elif exclude_fun is not None:
        exclude = exclude_fun(yss)
    else:
        exclude = None

    if s_method_freq == 'PLS':
        # use arPLS, lam from s_sigma ( s_sigma_freq)
        yss = get_baseline_mp(nproc)(freq, yss, axis=1, method='arPLS', bl_para={
            'lam': para['s_sigma'], "offset": 2, 'deg': 2}, verbose=verbose)
        s_method_freq = None
    # get baseline
    bls = get_baseline_mp(nproc)(freq, yss, axis=1, exclude=exclude,
                                 exclude_add=exclude_add,
                                 s_method=s_method_freq, s_sigma=s_sigma_freq, average_every=average_every_freq,
                          method=method, bl_para=bl_para, verbose=verbose)

    if njoin is not None and njoin > 1:
        bls = bls[np.searchsorted(np.arange(0, yss_ori.shape[0], njoin), 1+np.arange(yss_ori.shape[0]))-1]
    if subtract:
        return yss_ori - bls
    else:
        return bls
