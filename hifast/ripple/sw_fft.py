__all__ = ['replace_spec', 'save_rep_fixed', 'repalce_near', 'replace_margin_side', 'get_margin_num', 'replace_rfi_lower', 'replace_rfi',
           'find_loc', 'SW_FFT', 'get_sw_conf', 'fit_sw_fft',]

# author: Xu Chen, Li Fujia, 2021.06
# code: Xu Chen

import numpy as np
#import numpy.fft as fft
from scipy import fft
from astropy import log
from copy import deepcopy
from tqdm import tqdm
from matplotlib import pyplot as plt

# Cell
def replace_spec(tn, spec, freq, exceed_use, restrict_use, rfi_width_lim, ext_sec, ext, STD, MAX=None,
                 test=False, sm = None, fill = 'original', shift = 0,):
    """
    replace mw, rfi or others by near ripple section in one spectrum
    exceed_use: spec > thr | rfi
    restrict_use: restrict replace area in exceed_use
    ext: extend freq range to replace (channel number)
    rfi_width_lim: rfi should contain more channels than limit
    ext_sec: extend channel number of start and end of each section
    """
    if MAX is None:
        MAX = np.nanmax(spec)
    from .markRFI import rms, get_startend
    from ..utils.misc import smooth1d
    if np.sum(exceed_use) == 0:
        return spec

    start, end = get_startend(exceed_use, rfi_width_lim, ext_sec)
    if len(start) == 0:
#         print(f"tn={tn} raised a warning")
        return spec
    N = len(freq)
    newspec = deepcopy(spec)
    newspec_ = deepcopy(spec)
    if test:
        rep_from = np.zeros_like(spec)
        rep_to = np.zeros_like(spec, dtype='bool')

    for s, e in zip(start, end):
#         print(freq[s],freq[e])
        if fill == 'zero':
            newspec[s:e] = 0
            if test: rep_to[s:e] = 1
        else:
            # used to find two min values as sin valleys
            spec1 = np.zeros_like(freq) + MAX * 20
            s0 = s - ext
            s0 = max(s0, 0)
            e0 = e + ext
            e0 = min(e0, N)
    #         print(freq[s0],freq[e0])
            spec1[s0:e0] = sm[s0:e0]

            spec1l = deepcopy(spec1)
            spec1r = deepcopy(spec1)

            spec1l[s:] = MAX * 20
            l1 = np.argmin(spec1l) if s > 0 else 0
            spec1r[:e] = MAX * 20
            r1 = np.argmin(spec1r) if e < N else N

            if test:
                rep_to[l1:r1] = True
            # [l1:r1] will be replaced, length = L
            # print(l1,r1)
            L = r1 - l1
            if l1 - L - int(L * shift) < 0:
                direc = 'right'
            elif r1 + L + int(L * shift)> N - 1:
                direc = 'left'
            else:
                rmsl = rms(sm[l1 - L:l1], freq[l1 - L:l1])
                rmsr = rms(sm[r1:r1 + L], freq[r1:r1 + L])
                if rmsl < rmsr:
                    direc = 'left'
                else:
                    direc = 'right'

            if direc == 'left':
                s1 = l1 - L     - int(L * shift)
                e1 = l1         - int(L * shift)
            elif direc == 'right':
                s1 = r1         + int(L * shift)
                e1 = r1 + L     + int(L * shift)

            # print(direc,s1,e1)
            # [s1:e1] will be copied into [l1:r1]
            # fill what kinds of noise is not important, because they are high frequency compoents.
            if fill == 'original':
                ripple = spec[s1:e1]
            elif fill == 'smooth':
                ripple = sm[s1:e1]
            elif fill == 'artificial':
                ripple = sm[s1:e1] + np.random.normal(scale=STD, size=int(e1-s1))
            if test: rep_from[s1:e1] = 1
            try:
                newspec_[l1:r1] = ripple
                newspec[s:e] = newspec_[s:e] # only s to e will be changed
#                 print(s, e, l1, r1)
            except ValueError:
                log.info(f"tn = {tn} has a ValueError, will replace with noise")
                # import traceback
                # traceback.print_exc()
                newspec[s:e] = np.random.normal(scale=STD, size=int(e-s)) + np.nanmedian(spec)
    # restrict bound
    if not restrict_use.all(): newspec[~restrict_use] = spec[~restrict_use]
    if test:
        return newspec, rep_from, rep_to
    else:
        return newspec
    
def save_rep_fixed(sm, freq, exceed_use, rfi_width_lim, ext_sec, step, freq_d, freq_u):
    """
    save the extended replaced area, fixed freq step
    """
    from .markRFI import get_startend
    
    is_excluded = np.zeros_like(sm, dtype='bool') 
    if np.sum(exceed_use) == 0:
        return is_excluded

    start, end = get_startend(exceed_use, rfi_width_lim, ext_sec,)
    if len(start) == 0:
        return is_excluded
    
    N = len(freq)
    for s, e in zip(start, end):
        if e > freq_d and s < freq_u:
            if step > 0:
                s0 = s - step
                s0 = max(s0, 0)
                e0 = e + step
                e0 = min(e0, N)
            else:
                s0 = s 
                e0 = e
            is_excluded[s0:e0] = True
        
    return is_excluded
    
def save_rep_2sides(sm, freq, exceed_use, rfi_width_lim, ext_sec, ext, thr, freq_d, freq_u):
    """
    save the extended replaced area, from peak to 2 sides
    """
    from .markRFI import get_startend, find_edge_2sides

    is_excluded = np.zeros_like(sm, dtype='bool')
    if np.sum(exceed_use) == 0:
        return is_excluded

    start, end = get_startend(exceed_use, rfi_width_lim, ext_sec,)
    if len(start) == 0:
        return is_excluded

    N = len(freq)
    for s, e in zip(start, end):
        if e > freq_d and s < freq_u:
            s0 = s - ext * 5
            s0 = max(s0, 0)
            e0 = e + ext * 5
            e0 = min(e0, N)
            peak_position = freq[s:e][np.argmax(sm[s:e])]
#             if e0 < N: print(s,e,s0,e0,freq[s],freq[e],freq[s0],freq[e0], peak_position)
            mask_frange = find_edge_2sides(sm[s0:e0],freq[s0:e0],peak_position = peak_position,
                    small_rfi_times = 1, step=1*ext*fdelta, rms_thresh=thr,Print=False)
            if not np.isnan(mask_frange[0]):
                dist = np.max(np.abs(mask_frange - peak_position))
                mask_frange = peak_position + np.array([-1,1]) * dist
#                 print(mask_frange)
                is_excluded[(freq>=mask_frange[0])&(freq<=mask_frange[1])] = True
        
    return is_excluded

# Cell
def repalce_near(data_in, freq, time_rfi, mw_use=None, times_thr=None, times_s_thr=None, times_s_thr2=None,
                 rms_sigma=None,  ext_freq=None, rms_frange=None, rfi_width_lim=None,
                 ext_sec=None, data_find = None, data_trough = None, data_restrict = None,
                 save_is_excluded = False,ex_step = 0, verbose = True, **kwargs):
    """
    replace mw, rfi or others by near ripple section

    times_thr: above ~ times of rms will be set noise (unsmooth)
    times_s_thr: above ~ times of rms will be replaced (smoothed once)
    times_s_thr2: above ~ times of rms will be replaced (smooth twice)
    ext_freq: extend freq range to replace (mhz)
    rfi_width_lim: rfi should contain more channels than limit
    ext_sec: extend channel number of start and end of each section
    """
    global fdelta
    fdelta = np.abs(freq[1]-freq[0])
    ext = int(np.around(ext_freq / fdelta))
    
    data = deepcopy(data_in)
    
    from .markRFI import real_std, real_rms, std, rms
    if time_rfi is not None:
        whole_rfi = np.all(time_rfi, axis=1)
        is_rfi_num = np.arange(data.shape[0])[whole_rfi]
        not_rfi_num = np.arange(data.shape[0])[~whole_rfi]
    else:
        not_rfi_num = np.arange(data.shape[0])
        time_rfi = np.full(data.shape, False)

    if mw_use is None:
        mw_use = np.full(freq.shape, False)

    save_is_excluded = True if save_is_excluded and (times_s_thr2 is not None) else False

    data_rep = deepcopy(data)
    if data_find is None: data_find = data
    if times_s_thr2 is None: times_s_thr2 = times_s_thr
    if data_trough is None: data_trough = data_find
    if data_restrict is None:
        RESTRICT = False
    else:
        RESTRICT = False if data_restrict.shape != data.shape else True
    
    if save_is_excluded:
        is_excluded = np.zeros_like(data,dtype=bool)        
        Lfreq = len(freq)
        freq_d = int(Lfreq * 0.1)
        freq_u = int(Lfreq * 0.9)
        save_ex_step = int(np.around(ex_step / fdelta))
    else:
        is_excluded = np.array([None])

    from ..utils.misc import smooth1d
    MAX = np.nanmax(data)*20
    
    if np.sum(np.isnan(data)) > 0:
        data[np.isnan(data)] = MAX
        
    if verbose:
        iter_ = tqdm(range(data.shape[0]), desc='CPU 0: ', mininterval=2)
    else:
        iter_ = range(data.shape[0])
    
    for tn in iter_:
        if tn in not_rfi_num:
            spec = deepcopy(data[tn, :])
            include = mw_use #| time_rfi[tn]

            find = ex_sm = data_find[tn]
            # find used to define replace area
            RMS = rms(find, freq, rms_vrange=rms_frange)
            thr_s = RMS*times_s_thr2
            exceed_use = (np.abs(find) > thr_s) | time_rfi[tn]
            if np.sum(include) > 0:
                spec[include] = MAX * 20
                exceed_use |= include
            
            if RESTRICT:
                restrict = ex_sm = data_restrict[tn]
                # restrict used to restrict replace area in one spectra
                restrict_use = (restrict > thr_s) | time_rfi[tn]
            else:
                restrict_use = np.full(freq.shape, True)
            
            if save_is_excluded:
                # save sources or RFI positions
                is_excluded[tn,:] = save_rep_fixed(ex_sm, freq, exceed_use, rfi_width_lim, 
                                             ext_sec, save_ex_step, freq_d,freq_u)
#                 is_excluded[tn,:] = save_rep_2sides(ex_sm, freq, exceed_use, rfi_width_lim, 
#                                              ext_sec, ext, thr_s, freq_d,freq_u)

            # trough used to find ripples valleys
            trough = data_trough[tn]
            STD = real_std(spec, freq, sigma = rms_sigma, vrange = rms_frange)
            newspec = replace_spec(tn, spec, freq, exceed_use, restrict_use, rfi_width_lim, ext_sec,
                                   ext, STD, MAX, sm = trough, **kwargs)

            RMS = real_rms(spec, freq, sigma=rms_sigma, rms_vrange=rms_frange)
            strange = np.where(np.abs(newspec) > times_thr * RMS)[0]
            newspec[strange] = np.random.normal(scale=RMS, size=len(strange))

            data_rep[tn, :] = newspec

    if save_is_excluded:
        # two edges do not mask 
        is_excluded[:,:freq_d] = False
        is_excluded[:,freq_u:] = False
        
    return data_rep, is_excluded

# Cell

def replace_margin_side(data, freq, is_rfi, mw_use, margin_width=20, ext_freq=1.3, side='both'):
    margin = int(len(freq)/(freq[-1]-freq[0])*margin_width)

    N = data.shape[0]
    MAX = np.nanmax(data)*20
    pads = np.full((N, margin), MAX)

    fdelta = freq[1]-freq[0]
    fext1 = np.arange(freq[0]-margin*fdelta, freq[0], fdelta)
    fext2 = np.arange(freq[-1], freq[-1]+margin*fdelta, fdelta)
    print(f"Margin from {side} side(s) ...")

    if side == 'both':
        data_ = np.hstack((pads, data, pads))
        freq_ = np.hstack((fext1, freq, fext2))
    elif side == 'left':
        data_ = np.hstack((pads, data))
        freq_ = np.hstack((fext1, freq))
    elif side == 'right':
        data_ = np.hstack((data, pads))
        freq_ = np.hstack((freq, fext2))

    ext = int(np.around(ext_freq / fdelta))

    from hifas.ripple.markRFI import real_std
    whole_rfi = np.all(is_rfi, axis=1)
    is_rfi_num = np.arange(data.shape[0])[whole_rfi]

    for ti in tqdm(range(N)):
        if ti not in is_rfi_num:
            spec_ = data_[ti]
            exceed_use = (spec_ == np.inf)
            STD = real_std(spec_, freq_, sigma=6, vrange=[1390, 1400])
            restrict_use = np.full(freq.shape, False)
            newspec = replace_spec(ti, spec_, freq_, exceed_use, restrict_use,
                                   rfi_width_lim=0, ext_sec=1, ext=ext, STD = STD, MAX=MAX)
            newspec[-1] = 0
            data_[ti] = newspec

    if mw_use is not None:
        pads = np.full((N, margin), False)
        pads[is_rfi_num, :] = True

        pad = np.full(margin, False)
        if side == 'both':
            is_rfi = np.hstack((pads, is_rfi, pads))
            mw_use = np.hstack((pad, mw_use, pad))
        elif side == 'left':
            is_rfi = np.hstack((pads, is_rfi))
            mw_use = np.hstack((pad, mw_use))
        elif side == 'right':
            is_rfi = np.hstack((is_rfi, pads))
            mw_use = np.hstack((mw_use, pad))

    return data_, freq_, is_rfi, mw_use, margin

# Cell
def get_margin_num(margin, side):
    if side == 'both':
        margin1 = margin
        margin2 = -margin
    elif side == 'left':
        margin1 = margin
        margin2 = 0
    elif side == 'right':
        margin2 = margin
        margin1 = 0
    return margin1, margin2


def replace_rfi_lower(data, freq, method, times_lower_thr=3, times_lower=None,
                      rms_sigma=5, rms_frange=None):

    from .markRFI import real_rms
    RMS = real_rms(data[0, :], freq, sigma=rms_sigma, rms_vrange=rms_frange)

    thr_lower = RMS*times_lower_thr
    exceed_use = (np.abs(data) > thr_lower)
    from ..utils.misc import extend_Trues
    exceed_use = extend_Trues(exceed_use, ext_add=10, leng_lim=20, axis=-1)

    data_low = deepcopy(data)
    if method == 'lower':
        data_low[exceed_use] = data[exceed_use]/times_lower
    elif method == 'set_zeros':
        data_low[exceed_use] = 0
    elif method == 'set_noise':
        data_low[exceed_use] = np.random.normal(scale=RMS, size=data.shape)[exceed_use]
    return data_low


def replace_rfi(data, freq, time_rfi=None, method='near_ripple',
                mw_frange=None, **rep_args):

    log.info(f"Replace RFI with {method} method ...")
    if mw_frange is None:
        if (max(freq) <= 1419) | (min(freq) >= 1422):
            print("don't contain MW")
            mw_use = np.zeros_like(freq, dtype='bool')
        else:
            mw_use = None
    else:
        mw_use = (freq >= mw_frange[0]) & (freq <= mw_frange[1])
    is_excluded = None
    if method == 'near_ripple':
        data_rep, is_excluded = repalce_near(data, freq, time_rfi=time_rfi, mw_use=mw_use, **rep_args)
    elif (method == 'lower') or (method == 'set_zeros') or (method == 'set_noise'):
        data_rep = replace_rfi_lower(data, freq, method=method, **rep_args)
    else:
        raise ValueError("Unsupport replace RFI method!")

    return data_rep, is_excluded

# Cell
def find_loc(amp, x, amp_thr, xlim, rip_mhz=True, Print=True):
    """
    amp: fft amptitude
    x: fft xfreq
    amp_thr: amp threshold
    rip_mhz: bool
    xlim: find max between xlim[0] and xlim[1]
    """
    loc = None
    sw_mhz = None
    if rip_mhz:
        u = (amp > amp_thr) & (x > xlim[0]) & (x < xlim[1])

        if np.sum(u) > 0:
            amp_ = np.zeros_like(x)
            amp_[u] = amp[u]

            loc = np.argmax(amp_)
            sw_mhz = 1/x[loc]
            if Print:
                print(f'find ripple {sw_mhz:.5f} MHz,{1/sw_mhz:.5f} mu s,locate in {loc}')
        else:
            rip_mhz = False

    return loc, sw_mhz, rip_mhz

# Cell
class SW_FFT(object):
    def __init__(self, s2p, freq, workers=10):
        if s2p.ndim != 2:
            raise(ValueError('input s2p should be 2d'))
        if len(freq) != s2p.shape[1]:
            raise(ValueError('len(freq) need equals to s2p.shape[1]'))
        self.s2p = s2p
        self.freq = freq
        self.workers = workers

    def do_fft(self,):
        A_data = fft.rfft(self.s2p, n=len(self.freq), axis=1, workers=self.workers)
        #self.complex_num = A_data
        self.x = fft.rfftfreq(self.s2p.shape[1], self.freq[1]-self.freq[0])
        self.period = 1/self.x
        self.amp = np.abs(A_data)
        self.phi = np.angle(A_data,)

    @property
    def normal_amp(self):
        N = len(self.freq)
        data = deepcopy(self.amp) # K
        data[:,0] = data[:,0] / N
        data[:,1:] = data[:,1:] / (N / 2)
        return data * 1e3 # mK

    def gen_amp_thr_s(self, amp_thr_mean_factor=1.05, amp_thr_solo_factor=1.4,
                      is_excluded_mean=None):

        if is_excluded_mean is None:
            is_excluded_mean = np.full(len(self.amp), False, dtype=bool)
        self.amp_mean_t = np.nanmean(self.amp[~is_excluded_mean], axis=0)
        self.amp_thr_mean = np.nanmedian(self.amp_mean_t) * amp_thr_mean_factor

        self.amp_thr_solo = np.nanmedian(self.amp_mean_t) * amp_thr_solo_factor

    def find_sw_loc(self, xlims=[[.90, .95], [1.8, 1.9]]):
        """
        find whether a standing wave exists in the xlim range
        can run many times
        """
        xlims = np.array(xlims)
        if xlims.ndim == 1:
            xlims = xlims[None, :]
        self.bingo_list = []
        self.loc_list = []
        self.sw_mhz_list = []
        for xlim in xlims:
            loc, sw_mhz, bingo = find_loc(self.amp_mean_t, self.x, self.amp_thr_mean, rip_mhz=True, xlim=xlim)
            self.loc_list += [loc, ]
            self.sw_mhz_list += [sw_mhz, ]
            self.bingo_list += [bingo, ]

    def find_sw_chans(self, nchans,):
        """
        find channels between loc-nchan and loc+nchan
        run after self.find_sw_loc, can run many times
        """
        if np.isscalar(nchans):
            nchans = np.full(len(self.bingo_list), nchans)
        elif len(nchans) == 1:
            nchans = np.full(len(self.bingo_list), nchans[0])
        elif len(nchans) != len(self.bingo_list):
            raise(ValueError('the number nchans is not same with xlim'))
        if not hasattr(self, 'is_use_sw_chan'):
            self.is_use_sw_chan = np.full(len(self.x), False)
        for i, bingo in enumerate(self.bingo_list):
            if bingo:
                self.is_use_sw_chan |= (np.abs(np.arange(len(self.x))-self.loc_list[i]) < nchans[i])

    def choose_and_ifft(self, method='interpolate', inplace_amp=True, sw_base=True):
        """
        inplace_amp: make a copy of self.amp
        """
        amp = self.amp if inplace_amp else np.copy(self.amp)
        if not hasattr(self, 'is_use_sw_chan'):
            self.is_use_sw_chan = np.full(len(self.x), False)
        is_use_sw_chan_amp = self.is_use_sw_chan & (amp > self.amp_thr_solo)

        if method == 'interpolate':
            # interpolate
            from scipy.interpolate import interp1d
            for ti in tqdm(range(amp.shape[0])):
                x_mask = self.x[~is_use_sw_chan_amp[ti]]
                amp_solo_mask = amp[ti][~is_use_sw_chan_amp[ti]]
                amp_solo_interp = interp1d(x_mask, amp_solo_mask, kind='linear')  # ,fill_value="extrapolate")
                amp[ti][is_use_sw_chan_amp[ti]] -= amp_solo_interp(self.x[is_use_sw_chan_amp[ti]])
        elif method == 'all':
            # no need to process amp[is_use_sw_chan_amp]
            pass
        if sw_base:
            amp[:, 1:][~is_use_sw_chan_amp[:, 1:]] = 0
        else:
            amp[~is_use_sw_chan_amp] = 0
        self.sw = np.real(fft.irfft(amp*np.exp(1j*self.phi), n=len(self.freq), axis=1, workers=self.workers))
        #self.is_use_sw_chan_amp = is_use_sw_chan_amp

# Cell
def get_sw_conf(chan_wide, chan_narr):
    sw_conf = {
        '1mhz': {'xlims': [[.90, .95], [1.8, 1.9]], 'nchans': [chan_wide, chan_narr]},
        '2mhz': {'xlims': [[.5, .6], ], 'nchans': [chan_wide]},
        '0_04mhz': {'xlims': [[25, 26]], 'nchans': [chan_narr]},
    }
    return sw_conf


def fit_sw_fft(s1p, freq, nproc, amp_thr_mean_factor=1.05, amp_thr_solo_factor=1.4, is_excluded_mean=None,
               chan_wide=5, chan_narr=3, sw_periods=['1mhz', '2mhz', '0_04mhz'],
               sw_base=True, choose_method='all'):
    """
    amp_thr_mean_factor: above mean amptitude threshold will be chosed
    amp_thr_solo_factor: above solo amptitude threshold will be chosed
    is_excluded_mean: bool, exclude large RFI
    chan_wide: channel numbers near 1mhz to be chosed (wide)
    chan_narr: channel numbers near 1mhz to be chosed (narrow)
    sw_periods: remove ripple (1mhz: 1.08mhz, 2mhz:1.92mhz, 0_04mhz: 0.039 mhz)
    choose_method: use 'all' modes or 'interpolate' from nearby modes
    sw_base: if True, remove constant components / the base frequency (0 \mu s)
    """

    sw_conf = get_sw_conf(chan_wide, chan_narr)
    sw = SW_FFT(s1p, freq, nproc)
    sw.do_fft()
    sw.gen_amp_thr_s(amp_thr_mean_factor=amp_thr_mean_factor, amp_thr_solo_factor=amp_thr_solo_factor,
                     is_excluded_mean=is_excluded_mean)
    if 'none' not in sw_periods:
        for key in sw_periods:
            sw.find_sw_loc(xlims=sw_conf[key]['xlims'])
            sw.find_sw_chans(nchans=sw_conf[key]['nchans'])
    sw.choose_and_ifft(method=choose_method, inplace_amp=True, sw_base=sw_base)
    return sw.sw

