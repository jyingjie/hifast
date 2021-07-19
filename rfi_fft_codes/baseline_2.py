#!/usr/bin/env python
# coding: utf-8

# author: Xu Chen, Li Fujia, 2021.06
# code: Xu Chen

import numpy as np
import numpy.fft as fft
from astropy import log
from copy import deepcopy
from tqdm import tqdm

def filter_smooth(spec,fdelta,method = 'savgol',**kwargs):
    if method == 'savgol':
        from util import _round_up_to_odd_integer
        from scipy.signal import savgol_filter
        sg_window = kwargs['sg_window']
        sg_polyorder = kwargs['sg_polyorder']
        window_length = _round_up_to_odd_integer(sg_window/fdelta)
        sm = savgol_filter(spec,window_length = window_length ,polyorder=sg_polyorder,)
    elif method == 'gaussian':
        from hifast.utils.misc import smooth1d
        s_sigma = kwargs['s_sigma']
        sm = smooth1d(spec,axis = 0,sigma = s_sigma, method = 'gaussian')
    return sm

def replace_rfi_substract_old(data,freq,is_rfi,mw_use,sg_window,sg_polyorder, 
                          times_lower,times_lower_thr,rms_sigma,rms_frange):
    """
    replace big RFI to reduce the impact in FFT
    
    Parameters:
    data: np array
    sg_window: scipy.signal.savgol_filter window length (unit MHz)
    sg_polyorder: scipy.signal.savgol_filter polyorder
    times_lower: lower milky way area by divide a number
    """
    fdelta = freq[1]-freq[0] 
    
    whole_rfi = np.all(is_rfi,axis = 1)
    is_rfi_num = np.arange(data.shape[0])[whole_rfi]
    not_rfi_num = np.arange(data.shape[0])[~whole_rfi]

    data_rmrfi = np.full(data.shape,np.nan)

    from util import _round_up_to_odd_integer
    window_length = _round_up_to_odd_integer(sg_window/fdelta)
    
    from scipy.signal import savgol_filter
    for tn in tqdm(range(data.shape[0])):
        if tn in not_rfi_num:
            spec = data[tn,:]
            newspec = deepcopy(spec)
            newspec[mw_use] = 0
            sg = savgol_filter(newspec,window_length = window_length ,polyorder=sg_polyorder,)
            newspec[is_rfi[tn]] = spec[is_rfi[tn]] - sg[is_rfi[tn]]
            data_rmrfi[tn,:] = newspec

    from markRFI import real_rms
    RMS = real_rms(data_rmrfi[0,:],freq,sigma=rms_sigma,rms_vrange=rms_frange)

    thr_lower = RMS*times_lower_thr
    low_use = (np.abs(data_rmrfi) > thr_lower)
    from hifast.utils.misc import extend_Trues
    low_use = extend_Trues(low_use,ext_add =10,leng_lim = 20,axis = -1)

    data_rmrfi_low = deepcopy(data_rmrfi)
    data_rmrfi_low[low_use] = data[low_use]/times_lower
    
    return data_rmrfi_low


def replace_rfi_substract(data,freq,is_rfi,short_rfi,mw_use,sg_window,sg_polyorder, s_sigma,
                          rms_sigma,rms_frange):
    """
    replace big RFI to reduce the impact in FFT
    
    Parameters:
    data: np array
    sg_window: scipy.signal.savgol_filter window length (unit MHz)
    sg_polyorder: scipy.signal.savgol_filter polyorder
    """
    fdelta = freq[1]-freq[0] 
    
    from markRFI import real_rms

    whole_rfi = np.all(is_rfi,axis = 1)
    is_rfi_num = np.arange(data.shape[0])[whole_rfi]
    not_rfi_num = np.arange(data.shape[0])[~whole_rfi]
    
    data_rmrfi = np.full(data.shape,np.nan)
    for tn in tqdm(range(data.shape[0])):
        if tn in not_rfi_num:
            spec = data[tn,:]
            RMS = real_rms(spec,freq,sigma=rms_sigma,rms_vrange=rms_frange)
            
            newspec = deepcopy(spec)
            gauss_sm = filter_smooth(spec,fdelta,method ='gaussian',s_sigma = s_sigma)
            newspec[mw_use] = 0
            sg_sm = filter_smooth(spec,fdelta,method = 'savgol',sg_window=sg_window,
                                  sg_polyorder=sg_polyorder)
            
            # period rfi and source
            newspec[is_rfi[tn]] = spec[is_rfi[tn]] - sg_sm[is_rfi[tn]]
            # mw 
            newspec[mw_use] = spec[mw_use] - gauss_sm[mw_use]
            if np.sum(short_rfi[tn]) > 0:
                # time RFI 
                newspec[short_rfi[tn]] = (spec[short_rfi[tn]] / gauss_sm[short_rfi[tn]] - 1)
                cond = short_rfi[tn] & (np.abs(newspec) > RMS * 3)
                newspec[cond] = spec[cond]
            strange = np.where(np.abs(newspec) > 3 * RMS)[0]
            newspec[strange] = np.random.normal(scale=RMS, size=len(strange))

            data_rmrfi[tn] = newspec

    return data_rmrfi

def replace_rfi_substract2(data,freq,short_rfi,mw_use,sg_window,sg_polyorder, s_sigma,
                          rms_sigma,rms_frange,times_low_thr,):
    """
    replace big RFI to reduce the impact in FFT
    
    Parameters:
    data: np array
    sg_window: scipy.signal.savgol_filter window length (unit MHz)
    sg_polyorder: scipy.signal.savgol_filter polyorder
    
    """
    fdelta = freq[1]-freq[0] 
    
    from markRFI import real_rms

    whole_rfi = np.all(short_rfi,axis = 1)
    is_rfi_num = np.arange(data.shape[0])[whole_rfi]
    not_rfi_num = np.arange(data.shape[0])[~whole_rfi]
    
    data_rmrfi = np.full(data.shape,np.nan)
    from hifast.utils.misc import extend_Trues

    for tn in tqdm(range(data.shape[0])):
        if tn in not_rfi_num:
            spec = data[tn,:]
            RMS = real_rms(spec,freq,sigma=rms_sigma,rms_vrange=rms_frange)
            
            thr = RMS*times_low_thr
            low_use = (spec > thr)    
            low_use = extend_Trues(low_use,ext_add =10,leng_lim = 20,axis = -1)
            
            newspec = deepcopy(spec)
            gauss_sm = filter_smooth(spec,fdelta,method ='gaussian',s_sigma = s_sigma)
            newspec[mw_use] = 0
            sg_sm = filter_smooth(spec,fdelta,method = 'savgol',sg_window=sg_window,
                                  sg_polyorder=sg_polyorder)
            
            # period rfi and source
            newspec[low_use] = spec[low_use] - sg_sm[low_use]
            # mw 
            newspec[mw_use] = spec[mw_use] - gauss_sm[mw_use]
            if np.sum(short_rfi[tn]) > 0:
                # time RFI 
                newspec[short_rfi[tn]] = (spec[short_rfi[tn]] / gauss_sm[short_rfi[tn]] - 1)
                cond = short_rfi[tn] & (np.abs(newspec) > thr * 2)
                newspec[cond] = spec[cond]
            strange = np.where(np.abs(newspec) > 3 * RMS)[0]
            newspec[strange] = np.random.normal(scale=RMS, size=len(strange))

            data_rmrfi[tn] = newspec

    return data_rmrfi

def replace_rfi_lower(data,freq,method,times_lower_thr=3,times_lower=None,
                      rms_sigma = 5,rms_frange=None):

    from markRFI import real_rms
    RMS = real_rms(data[0,:],freq,sigma=rms_sigma,rms_vrange=rms_frange)

    thr_lower = RMS*times_lower_thr
    low_use = (np.abs(data) > thr_lower)
    from hifast.utils.misc import extend_Trues
    low_use = extend_Trues(low_use,ext_add =10,leng_lim = 20,axis = -1)
    
    data_low = deepcopy(data)
    if method == 'lower':
        data_low[low_use] = data[low_use]/times_lower
    elif method == 'set zeros':
        data_low[low_use] = 0
        
    return data_low

def replace_rfi(data,freq,is_rfi = None,short_rfi = None, method='subtract with pd',
                mw_use = None, **rep_args):

    log.info(f"Replace RFI with {method} method ...")
    
    if method == 'subtract trpdr':
        data_rmrfi = replace_rfi_substract(data,freq,is_rfi,short_rfi,mw_use,**rep_args)
    elif method == 'subtract tr':
        data_rmrfi = replace_rfi_substract2(data,freq,short_rfi,mw_use,**rep_args)
    elif (method == 'lower') or (method == 'set zeros'):
        data_rmrfi = replace_rfi_lower(data,freq,method=method,**rep_args)
    elif method == 'subtract rfi':
        data_rmrfi = replace_rfi_substract_old(data,freq,is_rfi,mw_use,**rep_args)
    else:
        raise ValueError("Unsupport replace RFI method!")
    return data_rmrfi  


def fft_fit_ripple(data_rmrfi, freq,is_rfi_num,not_rfi_num,ori_shape,
                   sw_freq= 0.9254,amp_thr= None,sw_n = 5,rfi_8mhz = False,rfi_freq_step = None,
                   plot = False,pdf = None,title = None):
    """
    fit baseline ripple (standing wave) by FFT
    Parameter:
    data_rmrfi: data array (after replace RFI)
    sw_freq: standing wave 'frequence' in Fourier space, unit \mu s
    rfi_freq_step: big RFI residual influence in Fourier space, nearly 1/16.2 \mu s.
                    (period is 16.2 MHz)
    amp_thr: above amptitude threshold will be chosed.
    sw_n: channel number next to sw_freq will be chosed all.
    
    """
    fdelta = freq[1]-freq[0]
    
    if len(is_rfi_num) > 0:
        data_rmrfi = np.delete(data_rmrfi,is_rfi_num,axis = 0)

    fftf = FFT(data_rmrfi, freq)
    x = fftf.x
    amp_data = fftf.amp
    loc1 = np.argmin(np.abs((x - sw_freq)))
    use = (np.abs(np.arange(len(x))-loc1) < sw_n) & (amp_data >= amp_thr) 

    amp_data_inpd = deepcopy(amp_data)

    if rfi_8mhz:
        x_loc16 = np.arange(0,x[-1],rfi_freq_step)

        N = freq.shape[0]
        fs = 1/fdelta
        loc16 = np.around(x_loc16 / (fs/N)).astype('int')
        # two sides
        #loc16_ = np.hstack((loc16 - 1,loc16,loc16 + 1))
        #loc16_.sort()

        use16 = np.full(len(x),False)
        use16[loc16] = True;use16[loc1] = False ; use16[0] = False
        use16 = use16 &  (amp_data >= amp_thr) 

        from scipy.interpolate import interp1d
        x_ = deepcopy(x)
        # interplate
        for ti in tqdm(range(amp_data.shape[0])):
            x_mask = x_[~use16[ti]]
            amp_data_mask = amp_data[ti][~use16[ti]]
            amp_data_interp = interp1d(x_mask,amp_data_mask,kind='linear')#,fill_value="extrapolate")
            amp_data_inpd[ti][use16[ti]] = amp_data_interp(x_[use16[ti]])

    amp_data_inpd = amp_data - amp_data_inpd
    amp_data_inpd[:,0] = amp_data[:,0]

    amp_data_inpd[use] = amp_data[use]
    amp_data_inpd[amp_data_inpd < 0] = 0
    
    A_data_inpd =  amp_data_inpd * np.exp(1j*fftf.phi)
    A_data_ifft = np.real(np.fft.irfft(A_data_inpd,n=data_rmrfi.shape[1]))
    
    if plot:
        tn = not_rfi_num[10]
        if pdf is not None:
            from matplotlib import pyplot as plt
            plt.switch_backend('agg')
        is_ = (x >= 0) & (x <= 3)
        fig = plt.figure(figsize=(22,5))
        ax = fig.add_subplot(111) 
        ax.stem(x[is_],amp_data[tn,is_],linefmt='--',markerfmt ='C0o',label = 'whole fft')
        ax.stem(x[is_],amp_data_inpd[tn,is_],linefmt='--',markerfmt ='C1o',label = 'sw fft')
        ax.text(x[loc1],.1,f'{loc1}')
        ax.plot(x[use[tn]],np.zeros_like(x[use[tn]]),'co',label = f'near {x[loc1]:.4f} $\mu$s ')
        if rfi_8mhz:
            ax.plot(x[loc16],np.zeros_like(loc16),'gs',label = f'step {rfi_freq_step:.4f} $\mu$s ')
        ax.axhline(amp_thr,color = 'k',linestyle='--')
        ax.grid();ax.set_xlim(0,3)
        ax.set_xlabel('k [$\mu$s]')
        plt.legend()
        plt.tight_layout()
        if title is not None:
            ax.set_title(title)
        if pdf is not None:
            pdf.savefig();plt.close()

    if len(is_rfi_num)>0:
        A_data_ifft_ = np.full(ori_shape,np.nan)
        A_data_ifft_[not_rfi_num,:] = A_data_ifft
        A_data_ifft = deepcopy(A_data_ifft_)
    
    return A_data_ifft

def mean_fit_ripple(data, nspec):
    n = nspec // 2
    
    tlen = data.shape[0]
    t1 = np.arange(tlen)-n
    t2 = np.arange(tlen)+n
    t1[t1 < 0] = 0;t2[t2 < 2*n] = 2*n
    t1[t1 > tlen - 2*n -1] = tlen - 2*n -1;t2[t2 > tlen -1] = tlen -1
    
    sw = np.zeros_like(data)
    if len(data.shape) == 2:
        for i in tqdm(range(tlen)):
            sw[i,:] = np.mean(data[t1[i]:t2[i],:],axis = 0)
    elif len(data.shape) == 3:
        for i in tqdm(range(tlen)):
            sw[i,:,0] = np.mean(data[t1[i]:t2[i],:,0],axis = 0)
            sw[i,:,1] = np.mean(data[t1[i]:t2[i],:,1],axis = 0)
            
    return sw
        
def fit_ripple(data,freq=None, method='rfft', is_rfi_num=None,not_rfi_num=None,ori_shape=None,nspec = None,
               plot = False,pdf = None,title = None,**fit_args): 
    
    log.info(f"Fit baseline ripple with {method} method ...")
    if method == 'rfft':
        A_data_ifft = fft_fit_ripple(data, freq,is_rfi_num,not_rfi_num,ori_shape,
                   plot = plot,pdf = pdf,title = title, **fit_args)
        
        return A_data_ifft        
    elif method == 'mean':
        sw = mean_fit_ripple(data,nspec = nspec )
        
        return sw
    

class FFT(object):
    def __init__(self, data,freq):
        # 采样点数
        N = freq.shape[0]
        # 采样周期T
        fdelta = freq[1] - freq[0]
        # 采样频率
        #fs = 1/fdelta
        
        # fft
        if len(data.shape) == 1:
            A_data = fft.rfft(data)
        else:
            A_data = fft.rfft(data,axis = 1)
            
        self.complex_num = A_data
        self.x = fft.rfftfreq(N,fdelta)
        self.amp = np.abs(A_data)
        self.phi = np.angle(A_data)        
        
