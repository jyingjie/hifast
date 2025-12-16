

__all__ = ['IS_Interactive', 'h5py_write', 'write_header', 'binned_statistic_e', 'split_continues', 'CheckPCal',
           'CalOnOffSav', 'CalOnOffM', 'CalOnOff1111', 'FASTRawCut']


import numpy as np
from scipy import interpolate as interp
from scipy import ndimage
import h5py

import os
import re
from glob import glob

from .cal import *
from ..utils.io import MjdChanPolar_to_PolarMjdChan, gen_carta_group, save_specs_hdf5

from astropy.stats import sigma_clip

import warnings

import __main__
IS_Interactive = not hasattr(__main__, '__file__')
if not IS_Interactive:
    import matplotlib
    matplotlib.use('Agg')
from matplotlib import pyplot as plt


import h5py
import os
def h5py_write(fname, mode='w'):
    try:
        f = h5py.File(fname, mode)
    except OSError:
        if mode == 'w':
            from datetime import datetime
            os.rename(fname, fname+'.del.empty.'+datetime.now().strftime("%Y%m%d-%H%M%S"))
            f = h5py.File(fname, mode)
    return f


def write_header(f, header):
    f.create_group('Header', track_order=False)
    for key2 in header.keys():
        f['Header'].attrs[key2] = header[key2]


def mean_a(arr):
    """
    calculate the mean value of 5th to 95th peercent of the input arr
    """
    b, e = np.percentile(arr, [5,95],)
    return np.mean(arr[(arr <=e ) & (arr >= b)], dtype='float64')

# def mean_a(arr):
#     """
#     calculate the mean value of 5th to 95th peercent of the input arr
#     """
#     return np.median(arr)


from scipy.stats import binned_statistic

def binned_statistic_e(x, values, axis=-1, **kwargs):
    """
    values: 3-Dim, stat at axis=1, stored in C-order
    """
    from collections import namedtuple
    res_tuple = namedtuple('BinnedStatisticResult',
                                   ('statistic', 'bin_edges', 'binnumber'))
    # move axis to last
    values = np.moveaxis(values, axis, -1)
    res = binned_statistic(x, values.reshape((-1, values.shape[-1])), **kwargs)
    res_stat = res[0].reshape(values.shape[:-1] + (-1,))
    # move back
    res_stat = np.moveaxis(res_stat, -1, axis)
    res = res_tuple(res_stat, *res[1:])
    return res


def test_binned_statistic_e():
    x = np.random.rand(1000)
    values = np.random.random((6, 1000, 2))
    assert (binned_statistic_e(x, values, axis=1, statistic='mean')[0][...,0] ==
            binned_statistic(x, values[...,0], statistic='mean')[0]).all()
    assert (binned_statistic_e(x, values, axis=1, statistic='mean')[0][...,1] ==
            binned_statistic(x, values[...,1], statistic='mean')[0]).all()

def test_binned_statistic_e_2():
    x = np.random.rand(1000)
    values = np.random.random((2, 6, 1000, 2))
    assert (binned_statistic_e(x, values, axis=2, statistic='mean')[0][0, ...,0] ==
            binned_statistic(x, values[0, ...,0], statistic='mean')[0]).all()
    assert (binned_statistic_e(x, values, axis=2, statistic='mean')[0][0, ...,1] ==
            binned_statistic(x, values[0, ...,1], statistic='mean')[0]).all()
    assert (binned_statistic_e(x, values, axis=2, statistic='mean')[0][1, ...,0] ==
            binned_statistic(x, values[1, ...,0], statistic='mean')[0]).all()
    assert (binned_statistic_e(x, values, axis=2, statistic='mean')[0][1, ...,1] ==
            binned_statistic(x, values[1, ...,1], statistic='mean')[0]).all()


def split_continues(bools):
    """
    bools: shape (m,); dtype bool
    """
    if bools.size == 0:
        raise('empty arr')
    diff = np.diff(bools.astype('int16'))
    ind_inv = np.where(diff != 0)[0] + 1

    ind_b = np.insert(ind_inv, 0, 0)
    ind_e = np.append(ind_inv, len(bools))-1

    return zip(bools[ind_b], ind_b, ind_e)


def test_split_continues():
    for arr_bools in [np.array([False] + [True]*3 + [False]*5 + [True]),
              np.array([True] + [False]*5 + [True] + [False]*3),
              np.array([False]*10),
             ]:
        segs = []
        for val, b, e in split_continues(arr_bools):
            seg = arr_bools[b:e+1]
            assert (seg == val).all()
            segs += [seg,]
            print(val, seg)
        assert (np.hstack(segs) == arr_bools).all()


class CheckPCal(CalOnOff):
    """
    """
    def __init__(self, *args, freq_step_c=2,
                              pcal_vary_lim_bin=0.01,
                              pcal_bad_lim_t=0.5,
                              pcal_bad_lim_freq=0.7,
                              **kwargs):
        # self.pcal_vary_frac = pcal_vary_frac
        self.pcal_vary_lim_bin = pcal_vary_lim_bin
        self.pcal_bad_lim_t = pcal_bad_lim_t # not used now
        self.pcal_bad_lim_freq = pcal_bad_lim_freq
        self.freq_step_c = freq_step_c
        super().__init__(*args, **kwargs)

    def _get_vary_frac(self, inds, freq_bins_c, bin_stat='mean'):
        """
        gen change fraction
        """
        if inds.ndim != 2:
            raise(ValueError('inds should be 2D.'))
        if bin_stat == 'mean_a':
            bin_stat = mean_a
        p = self.get_field(inds)
        ps, _, _ = binned_statistic_e(self.freq_use, p, axis=2, bins=freq_bins_c, statistic=bin_stat)
        ps_max, ps_min = np.max(ps, axis=1), np.min(ps, axis=1)
        df = (ps_max - ps_min)/ps_min
        return df


    def check_pcal(self, inds_ton, freq_step_c=2):
        """
        check pcal in freq bins
        """
        n_c = 6
        assert n_c - n_c//2 <= self.n_off

        n_cbef = n_c//2
        n_caft = n_c - n_cbef
        inds_cbef = inds_ton[:,:1] - np.arange(1,n_cbef+1)[::-1]
        inds_caft = inds_ton[:,-1:] + np.arange(1,n_caft+1)

        # fix first
        is_exceed = inds_cbef[0] < 0
        n_exceed = len(np.where(is_exceed)[0])
        if n_exceed >0:
            inds_cbef[0][is_exceed] = inds_caft[0][-1] + np.arange(1,n_exceed+1)
            inds_cbef[0][np.isin(inds_cbef[0], self.inds_on)] = inds_caft[0][-1]
            inds_cbef[0].sort()
        # fix last
        is_exceed = inds_caft[-1] > self.inds[-1]
        n_exceed = len(np.where(is_exceed)[0])
        if n_exceed >0:
            inds_caft[-1][is_exceed] = inds_cbef[-1][0] - np.arange(1,n_exceed+1)[::-1]
            inds_caft[-1][np.isin(inds_caft[-1], self.inds_on)] = inds_cbef[-1][0]
            inds_caft[-1].sort()
        inds_coff = np.hstack([inds_cbef, inds_caft])
        # been sorted !
        inds_coff[0].sort()
        inds_coff[-1].sort()

        # check freq bins
        freq_range = self.freq_use[-1] - self.freq_use[0]
        freq_delta = self.freq_use[1] - self.freq_use[0]
        bins_num = 1 + int(np.round(freq_range/freq_step_c))
        if bins_num < 2:
            bins_num = 2
        freq_bins_c = np.linspace(self.freq_use[0]-freq_delta/10, self.freq_use[-1]+freq_delta/10, bins_num)
        # check off
        is_bad_fbins = self._get_vary_frac(inds_coff, freq_bins_c) > self.pcal_vary_lim_bin
        # check on
        if inds_ton.shape[1] >= 2:
            is_bad_fbins = is_bad_fbins | (self._get_vary_frac(inds_ton, freq_bins_c) > self.pcal_vary_lim_bin)
        return freq_bins_c, is_bad_fbins

    def gen_pcals(self,):
        freq_bins_c, is_bad_fbins = self.check_pcal(self.inds_ton, self.freq_step_c)
        is_aband_whole = np.sum(is_bad_fbins, axis=1) / is_bad_fbins.shape[1] > self.pcal_bad_lim_freq
        ind_in_bins = np.digitize(self.freq_use, freq_bins_c) - 1

        pcals = self._get_cal_power(self.inds_ton, self.inds_toff_bef, self.inds_toff_aft)
        if getattr(self, 'plot_pcals', False):
            from ..waterfall import plot_im
            fig, axs = plt.subplots(2, 2, figsize=(16,10))
            plot_im(pcals[..., 0], x=self.freq_use, ax=axs[0][0], interpolation='none')
            plot_im(pcals[..., 1], x=self.freq_use, ax=axs[0][1], interpolation='none')
            axs[0][0].set_title('pcals ori, polar 0')
            axs[0][1].set_title('pcals ori, polar 1')
        fw = getattr(self, 'fw_pcals', None)
        if fw is not None:
            fw['pcals_ori'] = pcals.astype('float32')
        # replace as nan
        pcals[is_bad_fbins[:, ind_in_bins]] = np.nan
        for i in range(is_aband_whole.shape[1]):
            pcals[is_aband_whole[:,i], :, i] = np.nan
        if getattr(self, 'plot_pcals', False):
            plot_im(pcals[..., 0], x=self.freq_use, ax=axs[1][0], interpolation='none')
            plot_im(pcals[..., 1], x=self.freq_use, ax=axs[1][1], interpolation='none')
            axs[1][0].set_title('pcals masked, polar 0')
            axs[1][1].set_title('pcals masked, polar 1')
            fig.suptitle(os.path.basename(self.out_name_base + '-pcals.png'))
            fig.tight_layout()
            fig.savefig(self.out_name_base + '-pcals.png', dpi=111)
        if fw is not None:
            fw['pcals'] = pcals.astype('float32')
        self.pcals = pcals
        self.freq_bins_c = freq_bins_c
        self.is_bad_fbins = is_bad_fbins
        self.ind_in_bins = ind_in_bins
        self.is_aband_whole = is_aband_whole


class CalOnOffSav():
    def __call__(self, outdir='./', step=None, header=None, sep_save=False, save_pcals=False, cali=True):
        """
        get T, mjd etc, and save in hdf5 file

        Parameters
        ----------
        outdir : str
            output directory
        step : int
            number of chunk files to process each time
        header : dict
            add in the output file
        sep_save : bool
            if True, save file every step.
        """
        #
        self.gen_out_name_base(outdir)
        n_step = self.lens[0]*step if step is not None else len(self.inds)
        inds_range = np.append(np.arange(0, self.inds[-1], n_step), self.inds[-1]+1)
        # inds_splited= np.array_split(self.inds, len(self.inds)//n_step)
        mjds = []
        extra = self.get_extra()
        # tcal
        tc_inter = self.get_Tcal_s() if cali else 'none'
        # update and write header
        header = {} if header is None else header.copy()
        if hasattr(self, 'tcal_file'):
            header.update({'tcal_file': self.tcal_file})
        # prepare for self.get_count_tcal

        self.plot = True
        self.plot_pcals = True

        if cali:
            if save_pcals:
                outname = self.out_name_base + '-pcals.hdf5'
                fwp = h5py_write(outname)
                write_header(fwp, header)
                fwp_S = fwp.create_group('S')
                self.fw_pcals = fwp_S
            self.prepare_pcals()
#             if len(inds_ton_not_use)/(len(inds_ton_use) + len(inds_ton_not_use)) > 1/3.:
#                 warnings.warn('More than one-third of the Cals are abandoned. Please check your data and ``--pcal_vary_frac``.')
            if save_pcals:
                print(f"saving pcals in {outname}")
                fwp_S['inds_ton'] = self.inds_ton
                fwp_S['is_aband_whole'] = self.is_aband_whole
                for attr in ['pcals_merged',
                             'pcals_merged_s',
                             'pcals_amp_diff_interp_values',
                            ]:
                    if hasattr(self, attr):
                        fwp_S[attr] = getattr(self, attr)
                fwp.close()
                try:
                    del self.pcals
                except AttributeError:
                    pass

        for i in range(len(inds_range)-1):
            print('part', i)
            b, e = inds_range[i:i+2]
            inds_on = self.inds_on[(self.inds_on >= b) & (self.inds_on < e)]
            inds_off = self.inds_off[((self.inds_off >= b) & (self.inds_off < e))]
            inds = np.hstack([inds_on, inds_off])
            sort = np.argsort(inds)
            mjd = self.get_field(inds[sort], field='UTOBS')  # have sorted

            if cali:
                count_tcal_res = self.get_count_tcal(inds_on, inds_off)
            else:
                count_tcal_res = [self.get_field(inds_on, 'DATA',), self.get_field(inds_off, 'DATA', close_file=True)]

            T = np.vstack(count_tcal_res[:2])
            del count_tcal_res
            if cali:
                T = T*tc_inter
            T = T[sort]  # sort T
            T = T.astype('float32')  # finally convert to float32

            if sep_save:
                res = {}
                for key in extra.keys():
                    res[key] = extra[key][inds[sort]]
                res['mjd'] = mjd
                res['freq'] = self.freq_use
                res['Ta'] = MjdChanPolar_to_PolarMjdChan(T)
                res['Tcal'] = tc_inter
                res['inds_ton'] = self.inds_ton
                res['is_aband_whole'] = self.is_aband_whole
                for attr in ['pcals_merged',
                             'pcals_merged_s',
                            ]:
                    if hasattr(self, attr):
                        res[attr] = getattr(self, attr)
                if hasattr(self, 'pcals_amp_diff_interp_values'):
                    res['pcals_amp_diff_interp_values'] = self.pcals_amp_diff_interp_values[inds[sort]]

                # print(res)
                outname = self.out_name_base + f"-specs_T_{i:04d}_{i+1:04d}.hdf5"
                res['Header'] = header
                save_specs_hdf5(outname, res, wcs_data_name='Ta')
                print(f"Saved to {outname}")
                del res
            else:
                T = MjdChanPolar_to_PolarMjdChan(T)
                if i == 0:
                    outname = self.out_name_base + f"-specs_T.hdf5"
                    fout = h5py_write(outname)
                    #
                    fout.create_group('S')
                    g = fout['S']
                    # prepare writing spec
                    d_shape = list(T.shape)
                    d_shape[1] = len(self.inds)
                    g.create_dataset('Ta', shape=d_shape, dtype=T.dtype, chunks=True)
                g['Ta'][:, b:e, :] = T
                fout.flush()
                mjds += [mjd, ]
            self.plot = False  # only plot for the first loop
        # save other small data
        if not sep_save:
            write_header(fout, header)
            for key in extra.keys():
                g[key] = extra[key]
            g['mjd'] = np.hstack(mjds)
            g['freq'] = self.freq_use
            g['Tcal'] = tc_inter
            g['inds_ton'] = self.inds_ton
            g['is_aband_whole'] = self.is_aband_whole
            for attr in ['pcals_merged',
                         'pcals_merged_s',
                         'pcals_amp_diff_interp_values',
                        ]:
                if hasattr(self, attr):
                    g[attr] = getattr(self, attr)

            print(f"Saved to {outname}")
            gen_carta_group(fout, g['Ta'].shape, wcs_data_name='Ta', axis1=g['freq'], axis2=g['mjd'])
            fout.close()


class CalOnOffM(CheckPCal, CalOnOffSav, CalOnOff):

    def gen_squeeze_freq_is_use(self, ):
        """selecting channel used to calculate relative amp"""
        frange = self.squeeze_diff_freq_frange
        self.squeeze_freq_is_use = []
        for i in range(self.is_bad_fbins.shape[-1]):
            # first exclude the whole abandoned
            is_bad_fbins = self.is_bad_fbins[..., i][~self.is_aband_whole[:,i]]
            # select fbins with bad frac less than ...
            frac = np.sum(is_bad_fbins, axis=0)/is_bad_fbins.shape[0]
            bin_is = frac <= self.squeeze_diff_freq_bad_lim
            #bin_is = ~np.logical_or.reduce(self.is_bad_fbins[..., i][~self.is_aband_whole[:,i]], axis=0)
            is_ = bin_is[self.ind_in_bins]
            if frange is not None:
                is_ &= (self.freq_use > frange[0]) & (self.freq_use < frange[1])
            self.squeeze_freq_is_use += [is_,]

    def squeeze_freq(self, arr):
        """
        squeeze freq axis
        """
        method = self.squeeze_diff_freq

        if arr.ndim !=3 :
            raise(ValueError('need 3d'))
        if method == 'median':
            squ_fun = np.nanmedian
        elif method == 'mean':
            squ_fun = np.nanmean
        else:
            raise(ValueError(f'method {method} not supported'))
        arr_s = np.zeros((arr.shape[0], arr.shape[2]), dtype='float64')
        for i in range(arr.shape[2]):
            is_ = self.squeeze_freq_is_use[i]
            arr_s[:,i] = squ_fun(arr[:, is_, i], axis=1)
        return arr_s

    def merge_pcals(self,):
        """
        merging pcals at all time and smoothing it
        """
        method_merge = self.method_merge
        # 'scale' or 'none'
        pre_process = self.merge_cal_pre_process

        pcals = self.pcals
        if pre_process == 'scale':
            # use the power not is_bad_fbins to sacle
            amp = self.squeeze_freq(pcals)
            # not change self.pcals
            if self.calc_diff_method == 'div':
                pcals = pcals * (np.nanmedian(amp, axis=0) / amp)[:,None,:]
            elif self.calc_diff_method == 'sub':
                pcals = pcals + (np.nanmedian(amp, axis=0) - amp)[:,None,:]
            # keep scaled
            self.pcals_scaled = pcals
            self.pcals_scaled_amp = amp

        if method_merge == 'median':
            pcals_merged = np.nanmedian(pcals, axis=0, keepdims=True)
        elif method_merge == 'mean':
            pcals_merged = np.nanmean(pcals, axis=0, dtype='float64', keepdims=True)

        self.pcals_merged = pcals_merged
        self.pcals_merged_s = self._get_smoothed(pcals_merged, check_nan=True)

    def gen_pcals_amp_diff(self,):
        """
        run self.merge_pcals first to get self.pcals_merged and self.pcals_merged_s
        gen self.pcals_amp_diff

        Parameters
        ----------
        squeeze_diff_freq: str; 'median', 'mean', 'sigclip_median', 'sigclip_mean'
        exclued_franges: str; freq ranges which are not used in calculating amplitude offset; e.g. [[1375, 1385],]
        """
        pcals = self.pcals
        calc_diff_method = self.calc_diff_method

        if calc_diff_method == 'div':
            pcals_diff = pcals / self.pcals_merged_s
        elif calc_diff_method == 'sub':
            pcals_diff = pcals - self.pcals_merged_s
        else:
            raise(ValueError(f"calc_diff_method:{calc_diff_method}"))
        self.pcals_amp_diff = self.squeeze_freq(pcals_diff)
        self.pcals_amp_diff_inds = np.mean(self.inds_ton, axis=1)
        self.calc_diff_method = calc_diff_method

    def gen_pcals_amp_diff_interp_values(self,):
        """
        interplate pcals_amp_diff
        gen self.pcals_amp_diff_interp_values: the amplitude offset relate to mergred pcals, shape is (inds, polar)

        Parameters
        ----------
        method_interp: str; 'gaussian' 'CubicSpline', 'CubicHermiteSpline', kind in scipy.interpolate.interp1d
                            if 'polyXd'('X'is a integer), use poly fitting with deg=X
        edges: 'nearest', 'extrapolate'
        """
        pcals_amp_diff = self.pcals_amp_diff
        method_interp = self.method_interp
        sigma_t = self.method_interp_sigma_t
        edges = self.method_interp_edges

        self.pcals_amp_diff_interp_values = np.zeros((len(self.inds), pcals_amp_diff.shape[1]), dtype='float64')
        for i in range(pcals_amp_diff.shape[1]):
            x = self.pcals_amp_diff_inds
            y = pcals_amp_diff[:, i]
            is_ = np.isfinite(y)
            x = x[is_]
            y = y[is_]
            if method_interp == 'gaussian':
                t_sample = self.hduls[0][1].data['EXPOSURE'][0]
                sigma = sigma_t/((self.n_on+self.n_off)*t_sample)/(t_sample*self.n_on)
                y = ndimage.gaussian_filter1d(y, sigma, axis=0)
                fun = interp.interp1d(x, y, fill_value="extrapolate", kind='quadratic', axis=0)
            elif method_interp == 'CubicSpline':
                fun = interp.CubicSpline(x, y, extrapolate=None)
            elif method_interp == 'CubicHermiteSpline':
                fun = interp.CubicHermiteSpline(x, y, extrapolate=None)
            elif method_interp.startswith('poly'):
                fun = np.poly1d(np.polyfit(x, y, deg=int(method_interp[4:-1])))
            else:
                # method_interp is one of "kind"
                fun = interp.interp1d(x, y, fill_value="extrapolate", kind=method_interp)
            self.pcals_amp_diff_interp_values[:, i] = fun(self.inds)
            if edges == "nearest":
                if (self.inds < x [0]).any():
                    self.pcals_amp_diff_interp_values[self.inds < x [0], i] = fun(x[0])
                if (self.inds > x [-1]).any():
                    self.pcals_amp_diff_interp_values[self.inds > x [-1], i] = fun(x[-1])



    def set_para_pcals(self, calc_diff_method='div',
                       squeeze_diff_freq='median',
                       squeeze_diff_freq_frange=None,
                       squeeze_diff_freq_bad_lim=1,
                       method_interp='quadratic',
                       method_interp_edges='nearest',
                       method_interp_sigma_t=600,
                       method_merge='median',
                       merge_cal_pre_process='scale',
                       ):
        """
        para used in self.prepare_pcals

        calc_diff_method: str; 'div', 'sub'
        """
        self.calc_diff_method = calc_diff_method # in self.pcals_amp_interp & self.get_count_tcal
        self.squeeze_diff_freq = squeeze_diff_freq
        self.squeeze_diff_freq_frange = squeeze_diff_freq_frange
        self.squeeze_diff_freq_bad_lim = squeeze_diff_freq_bad_lim
        self.method_interp = method_interp
        self.method_interp_edges = method_interp_edges
        self.method_interp_sigma_t = method_interp_sigma_t
        self.method_merge = method_merge
        self.merge_cal_pre_process = merge_cal_pre_process

    def prepare_pcals(self,):
        self.gen_pcals()
        self.gen_squeeze_freq_is_use()
        self.merge_pcals()
        self.gen_pcals_amp_diff()
        self.gen_pcals_amp_diff_interp_values()
        self.plot_pcal_merged()
        if not IS_Interactive:
            if hasattr(self, 'pcals_scaled'):
                del self.pcals_scaled

    def plot_pcal_merged(self,):
        if not getattr(self, 'plot_pcals', False):
            return
        fig, axs = plt.subplots(3,2,figsize=(15,15))

        if hasattr(self, 'pcals_scaled'):
            pcals = self.pcals_scaled
        else:
            pcals = self.pcals
        try:
            bins = self.freq_bins_c
        except:
            try:
                bwidth = 1.2 * self.s_para['s_sigma']
            except:
                bwidth = 3.5
            bins = np.arange(self.freq_use[0], self.freq_use[-1] + bwidth, bwidth)
        pcals_d, _, _ = binned_statistic_e(self.freq_use, pcals, bins=bins, axis=1, statistic='mean')
        freq_d, _, _ = binned_statistic_e(self.freq_use, self.freq_use, bins=bins, axis=0, statistic='mean')
        ax = axs[0][0]
        ax.plot(freq_d, pcals_d[:,:,0].T, alpha=0.5)
        ax.set_title(f'pcals[_scaled] downsample, polar 0')
        ax = axs[0][1]
        ax.plot(freq_d, pcals_d[:,:,1].T, alpha=0.5)
        ax.set_title(f'pcals[_scaled] downsample, polar 1')

        ax = axs[1][0]
        ax.plot(self.freq_use, self.pcals_merged[0,:,0], 'r', alpha=0.5)
        ax.plot(self.freq_use, self.pcals_merged_s[0,:,0], 'k', lw=2)
        ax.set_title(f'pcals merged, polar 0')
        ax.set_xlabel(f'Freq')
        if hasattr(self, 'squeeze_freq_is_use'):
            ax = ax.twinx()
            ax.plot(self.freq_use, self.squeeze_freq_is_use[0])
            ax.set_ylim(-0.3, 14)
        ax = axs[1][1]
        ax.plot(self.freq_use, self.pcals_merged[0,:,1], 'r', alpha=0.5)
        ax.plot(self.freq_use, self.pcals_merged_s[0,:,1], 'k', lw=2)
        ax.set_title(f'pcals merged, polar 1')
        ax.set_xlabel(f'Freq')
        if hasattr(self, 'squeeze_freq_is_use'):
            ax = ax.twinx()
            ax.plot(self.freq_use, self.squeeze_freq_is_use[1])
            ax.set_ylim(-0.3, 14)

        ax = axs[2][0]
        ax.plot(self.inds, self.pcals_amp_diff_interp_values[:,0])
        ax.scatter(self.pcals_amp_diff_inds, self.pcals_amp_diff[:,0], color='r')
        ax.set_title('pcal amp diff, polar 0')
        ax.set_xlabel(f'Index')
        ax = axs[2][1]
        ax.plot(self.inds, self.pcals_amp_diff_interp_values[:,1])
        ax.scatter(self.pcals_amp_diff_inds, self.pcals_amp_diff[:,1], color='r')
        ax.set_title('pcal amp diff, polar 1')
        ax.set_xlabel(f'Index')
        fig.suptitle(os.path.basename(self.out_name_base))
        fig.tight_layout()
        fig.savefig(self.out_name_base + '-pcals-merged.png')


    def get_count_tcal(self, inds_on, inds_off):
        """
        run self.prepare_pcals first
        get count of tcal
        """
        #
        # interp cal power amp
        amp_interp_on = self.pcals_amp_diff_interp_values[inds_on]
        amp_interp_off = self.pcals_amp_diff_interp_values[inds_off]
        # load on and off power
        c_on = self.get_field(inds_on, 'DATA',).astype('float64')
        c_off = self.get_field(inds_off, 'DATA', close_file=True).astype('float64')
        #print(p_on.dtype, p_off.dtype, p_cal_s.dtype)
        if getattr(self, 'plot', False):
            figname = self.out_name_base + "-sep.pdf"
            plot_sep(inds_on, inds_off, c_on, c_off, figname=figname)
        # count as pcal
        if self.calc_diff_method == 'div':
            c_on /= (self.pcals_merged_s * amp_interp_on[:, None, :])
            c_on -= 1  # have subtracted cal
            c_off /= (self.pcals_merged_s * amp_interp_off[:, None, :])
        elif self.calc_diff_method == 'sub':
            c_on /= (self.pcals_merged_s + amp_interp_on[:, None, :])
            c_on -= 1  # have subtracted cal
            c_off /= (self.pcals_merged_s + amp_interp_off[:, None, :])

        return c_on, c_off


class CalOnOff1111(CheckPCal, CalOnOffSav, CalOnOff):
    @staticmethod
    def split_along_a1d2(bools_2d):
        """
        """
        if bools_2d.ndim !=2 :
            raise('should input 2d array')

        ind_0, ind_1 = np.where(np.diff(bools_2d.astype('int16'), axis=1) !=0)
        ind_inv = np.unique(ind_1) + 1
        ind_1_b = np.insert(ind_inv, 0, 0)
        ind_1_e = np.append(ind_inv, bools_2d.shape[1])-1

        return ind_1_b, ind_1_e


    def plot_pcal_s(self,):
        if not getattr(self, 'plot_pcals', False):
            return
        from ..waterfall import plot_im
        fig, axs = plt.subplots(1, 2, figsize=(16,6))
        plot_im(self.pcals_s[..., 0], x=self.freq_use, ax=axs[0], interpolation='none')
        plot_im(self.pcals_s[..., 1], x=self.freq_use, ax=axs[1], interpolation='none')
        axs[0].set_title('pcals smoothed, polar 0')
        axs[1].set_title('pcals smoothed, polar 1')
        fig.suptitle(os.path.basename(self.out_name_base))
        fig.tight_layout()
        fig.savefig(self.out_name_base + '-pcals-smoothed.png')

    def set_para_pcals(self, cal_dis_lim=2):
        self.delat_t_lim = cal_dis_lim

    def prepare_pcals(self, ):
        self.gen_pcals()

        pcals_s = self._get_smoothed(self.pcals, check_nan=True)
        del self.pcals

        segs_be_list = []
        is_use_in_segs_list = []
        for polar in range(pcals_s.shape[-1]):
            bools_2d = np.isnan(pcals_s[:, :, polar])
            segs_be = self.split_along_a1d2(bools_2d)
            if len(segs_be[0]) > 2*len(self.freq_bins_c):
                raise(ValueError('something is wrong'))
            segs_be_list.append(segs_be)
            is_use_in_segs = []
            for b, e in zip(*segs_be):
                seg = bools_2d[:, b:e+1]
                nan_num = np.sum(seg, axis=1)
                if len(np.unique(nan_num)) > 2:
                    raise(ValueError('nan values splitting failed'))
                is_use_in_segs.append(nan_num == 0)
            is_use_in_segs_list.append(is_use_in_segs)
        self.pcals_s = pcals_s
        self.segs_be_list = segs_be_list
        self.is_use_in_segs_list = is_use_in_segs_list
        self.inds_ton_m = np.mean(self.inds_ton, axis=1)
        self.plot_pcal_s()
        fw = getattr(self, 'fw_pcals', None)
        if fw is not None:
            fw['pcals_s'] = self.pcals_s.astype('float32')


    def get_count_tcal(self, inds_on, inds_off):
        delat_t_lim = self.delat_t_lim
        pcals_s = self.pcals_s

        c_on = self.get_field(inds_on, 'DATA',).astype('float64')
        c_off = self.get_field(inds_off, 'DATA', close_file=True).astype('float64')

        if getattr(self, 'plot', False):
            figname = self.out_name_base + "-sep.pdf"
            plot_sep(inds_on, inds_off, c_on, c_off, figname=figname)

        inds_ton_m = self.inds_ton_m

        for polar in range(pcals_s.shape[-1]):
            is_use_in_segs = self.is_use_in_segs_list[polar]
            segs_be = self.segs_be_list[polar]
            for i, (b, e) in enumerate(zip(*segs_be)):
                inds_ton_m_use = inds_ton_m[is_use_in_segs[i]]
                # on
                inds_close = np.argmin(abs(inds_ton_m_use - inds_on[:, None]), axis=1)
                c_on[:, b:e+1, polar] /= pcals_s[:, b:e+1, polar][is_use_in_segs[i]][inds_close]
                ## mask as nan if cal and spec not close enough
                is_too_far = abs(inds_on - inds_ton_m_use[inds_close]) > (self.n_on + self.n_off)*delat_t_lim
                c_on[is_too_far, b:e+1, polar] = np.nan

                # on
                inds_close = np.argmin(abs(inds_ton_m_use - inds_off[:, None]), axis=1)
                c_off[:, b:e+1, polar] /= pcals_s[:, b:e+1, polar][is_use_in_segs[i]][inds_close]
                ## mask as nan if cal and spec not close enough
                is_too_far = abs(inds_off - inds_ton_m_use[inds_close]) > (self.n_on + self.n_off)*delat_t_lim
                c_off[is_too_far, b:e+1, polar] = np.nan
        c_on -= 1  # have subtracted cal
        return c_on, c_off


class FASTRawCut(FastRawData):
    def __init__(self,
                 fname_part,
                 start=1,
                 stop=None,
                 frange=None,
                 verbose=False,):
        """
        fname_part: str; one of chunk file name
        start: int; chunk file start
        stop: int; chunk file stop
        frange: frequency range
        """
        # generate fits filename list
        fname_part = re.sub('[0-9]{4}\.fits\Z', '', fname_part)
        if stop is None:
            stop = len(glob(fname_part+'*.fits'))
        self.fname_part = fname_part
        filenames = [fname_part+'%04d.fits'%i for i in range(start,stop+1)]
        if len(filenames) == 0:
            raise(OSError(f"can not find file, please check fname_part:{fname_part}"))
        super().__init__(filenames, frange=frange, verbose=verbose)

        self.nB = int(re.findall(r'-M[0-1][0-9]', fname_part)[-1][2:])

        self.field_names = self.hduls[0][1].data.dtype.names

        self.fits_header = dict(self.hd1s[0].items())

    def _get_freq(self,):
        # make sure not correct to process 'FREQ_new'
        super()._get_freq(center_corr=None)

    def gen_out_name_base(self, outdir):
        self.out_name_base = os.path.join(outdir, f"{os.path.basename(self.fname_part)}")

    def gen_new_freq_axis(self, ):
        is_ = (self.freq >= self.frange[0]) & (self.freq <= self.frange[1])
        self.NCHAN_new = len(np.where(is_)[0])
        self.FREQ_new = self.freq[is_][0]
        self.is_use_chan = is_

    def write_fits_header(self, g, **kwargs):
        """
        g: h5py group
        kwargs: header field pair: name:value
        """
        for key in kwargs.keys():
            g.attrs[key] = kwargs[key]


    def load_all_fields(self, inds, chan_fix=True, polar_trans=True):
        dict1 = {} # small
        dict2 = {} # large
        for key in self.field_names:
            if key != 'DATA':
                dict1[key] = self.get_field(inds, key)
                if dict1[key].dtype.char == 'U':
                    dict1[key] = dict1[key].astype('S')
            else:
                dict2[key] = self.get_field(inds, key, close_file=True)

        if chan_fix:
            dict1['NCHAN'][:] = self.NCHAN_new
            dict1['FREQ'][:] = self.FREQ_new
            assert dict2['DATA'].shape[1] == self.NCHAN_new
        if polar_trans:
            dict2['DATA'] = MjdChanPolar_to_PolarMjdChan(dict2['DATA'])

        return dict1, dict2

    def dict_stack(self, *args):
        """
        stack item in dicts, keys must be same
        """
        res = {}
        if len(args) == 0:
            return ()
        for key in args[0].keys():
            res[key] = np.hstack([d[key] for d in args])
        return res

    def add_S_group(self, f):
        """
        add a hdf5 group for spcetra type. add 'freq', hdf5 softlink to 'DATA', 'UTOBS' in '/1'
        f: file handle
        """
        g = f.create_group('S')
        # for W,N,F. also see FastRawSpec._get_freq
        center_corr = 0.000476837158203125/2.
        g['freq'] = self.freq_use + center_corr
        g['DATA'] = h5py.SoftLink('/1/DATA')
        g['Ta'] = h5py.SoftLink('/1/DATA')
        g['mjd'] = h5py.SoftLink('/1/UTOBS')


    def __call__(self, outdir='./', step=1, header=None, sep_save=False, h5_compression='none'):
        """

        Parameters
        ----------
        outdir : str
            output directory
        step : int
            number of chunk files to process each time
        header : dict
            add in the output file
        sep_save : bool
            if True, save file every step.
        """
        if h5_compression == 'none':
            h5_compression = None
        self.gen_out_name_base(outdir)
        self.gen_new_freq_axis()

        n_step = self.lens[0]*step if step is not None else len(self.inds)
        inds_range = np.append(np.arange(0, self.inds[-1], n_step), self.inds[-1]+1)
        # inds_splited= np.array_split(self.inds, len(self.inds)//n_step)
        mjds = []
        # update and write header
        header = {} if header is None else header.copy()

        dict1_list = []
        stype = 'DATA'
        for i in range(len(inds_range)-1):
            print('part', i)
            b, e = inds_range[i:i+2]
            inds = self.inds[((self.inds >= b) & (self.inds < e))]
            dict1, dict2 = self.load_all_fields(inds, chan_fix=True, polar_trans=True)
            # large data
            T = dict2[stype]
            if sep_save:
                outname = self.out_name_base + f"{i+1:04d}.hdf5"
                fout_sep = h5py_write(outname)
                write_header(fout_sep, header)
                g_sep_0 = fout_sep.create_group('0') # table 0
                self.write_fits_header(g_sep_0, **dict(self.hd0s[i*step].items()))
                g_sep = fout_sep.create_group('1') # table 1
                for key in dict1.keys():
                    g_sep[key] = dict1[key]
                g_sep.create_dataset(stype, data=T, chunks=True, compression=h5_compression)
                # NAXIS1=self.NCHAN_new*self.2 ?
                self.write_fits_header(g_sep, **self.fits_header)
                # update new
                self.write_fits_header(g_sep,
                                       NAXIS2=len(inds),
                                       TDIM21=f"(2, {self.NCHAN_new})",
                                      )
                # waterfall
                mjd = dict1['UTOBS']
                gen_carta_group(fout_sep, g_sep[stype].shape, wcs_data_name=stype, axis1=self.freq_use, axis2=mjd, wcs_data_group='1')
                # add S group
                self.add_S_group(fout_sep)
                fout_sep.close()
                print(f"Saved to {outname}")
            else:
                # open file
                if i == 0:
                    outname = self.out_name_base + f"0001.hdf5"
                    fout = h5py_write(outname)
                    g_0 = fout.create_group('0') # table 0
                    self.write_fits_header(g_0, **dict(self.hd0s[0].items()))
                    fout.create_group('1')
                    g = fout['1']
                    # prepare writing spec
                    d_shape = list(T.shape)
                    d_shape[1] = len(self.inds)
                    g.create_dataset(stype, shape=d_shape, dtype=T.dtype, chunks=True, compression=h5_compression)
                g[stype][:, b:e, :] = T
                fout.flush()
                dict1_list.append(dict1)
        # save other small data
        if not sep_save:
            write_header(fout, header)
            dict1 = self.dict_stack(*dict1_list)
            for key in dict1.keys():
                g[key] = dict1[key]
            self.write_fits_header(g, **self.fits_header)
                # update new
            self.write_fits_header(g,
                                   NAXIS2=len(self.inds),
                                   TDIM21=f"(2, {self.NCHAN_new})",
                                   )
            # waterfall
            mjd = dict1['UTOBS']
            gen_carta_group(fout, g[stype].shape, wcs_data_name=stype, axis1=self.freq_use, axis2=mjd, wcs_data_group='1')
            # add S group
            self.add_S_group(fout)
            fout.close()
            print(f"Saved to {outname}")
