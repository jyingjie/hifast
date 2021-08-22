#!/usr/bin/env python
# coding: utf-8

# author: Xu Chen, Li Fujia, 2021.06
# code: Xu Chen

import numpy as np
import numpy.fft as fft
from astropy import log
from copy import deepcopy
from tqdm import tqdm
from matplotlib import pyplot as plt

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
    """
    fdelta = freq[1]-freq[0] 
    
    if time_rfi is not None:
        whole_rfi = np.all(time_rfi,axis = 1)
        is_rfi_num = np.arange(data.shape[0])[whole_rfi]
        not_rfi_num = np.arange(data.shape[0])[~whole_rfi]
    else:
        not_rfi_num = np.arange(data.shape[0])
    data_rep = np.full(data.shape,np.nan)

    from util import _round_up_to_odd_integer
    window_length = _round_up_to_odd_integer(sg_window/fdelta)
    from markRFI import real_rms
    from scipy.signal import savgol_filter
    for tn in tqdm(range(data.shape[0])):
        if tn in not_rfi_num:
            spec = deepcopy(data[tn,:])
            RMS = real_rms(spec,freq,sigma=rms_sigma,rms_vrange=rms_frange)
            newspec = deepcopy(spec)
            newspec[mw_use] = 0
            sg = savgol_filter(newspec,window_length = window_length ,polyorder=sg_polyorder,)
            newspec[is_rfi[tn]] = spec[is_rfi[tn]] - sg[is_rfi[tn]]
            newspec[mw_use] = np.random.normal(scale=RMS,size = np.sum(mw_use)) 
            data_rep[tn,:] = newspec

    
    #RMS = real_rms(data_rep[0,:],freq,sigma=rms_sigma,rms_vrange=rms_frange)
    
    thr_lower = RMS*times_lower_thr
    low_use = (np.abs(data_rep) > thr_lower)
    from hifast.utils.misc import extend_Trues
    low_use = extend_Trues(low_use,ext_add =10,leng_lim = 20,axis = -1)

    data_rep = deepcopy(data_rep)
    data_rep[low_use] = data[low_use]/times_lower
    
    return data_rep


def replace_rfi_substract(data,freq,is_rfi,time_rfi,mw_use,sg_window,sg_polyorder, s_sigma,
                          rms_sigma,rms_frange,times_lower_thr = 5):
    """
    replace big RFI to reduce the impact in FFT
    
    Parameters:
    data: np array
    sg_window: scipy.signal.savgol_filter window length (unit MHz)
    sg_polyorder: scipy.signal.savgol_filter polyorder
    """
    fdelta = freq[1]-freq[0] 
    
    from markRFI import real_rms
    if time_rfi is not None:
        whole_rfi = np.all(time_rfi,axis = 1)
        is_rfi_num = np.arange(data.shape[0])[whole_rfi]
        not_rfi_num = np.arange(data.shape[0])[~whole_rfi]
    else:
        not_rfi_num = np.arange(data.shape[0])
    
    data_rep = np.full(data.shape,np.nan)
    for tn in tqdm(range(data.shape[0])):
        if tn in not_rfi_num:
            spec = deepcopy(data[tn,:])
            RMS = real_rms(spec,freq,sigma=rms_sigma,rms_vrange=rms_frange)
            
            newspec = deepcopy(spec)
            
            newspec[mw_use] = 0
            sg_sm = filter_smooth(spec,fdelta,method = 'savgol',sg_window=sg_window,
                                  sg_polyorder=sg_polyorder)
            # period rfi and source
            newspec[is_rfi[tn]] = spec[is_rfi[tn]] - sg_sm[is_rfi[tn]]
            # mw 
            newspec[mw_use] = np.random.normal(scale=RMS,size = np.sum(mw_use)) 
            if np.sum(time_rfi[tn]) > 0:
                gauss_sm = filter_smooth(spec,fdelta,method ='gaussian',s_sigma = s_sigma)
                # time RFI 
                newspec[time_rfi[tn]] = (spec[time_rfi[tn]] / gauss_sm[time_rfi[tn]] - 1)
                cond = time_rfi[tn] & (np.abs(newspec) > RMS * 2)
                newspec[cond] = spec[cond]
            strange = np.where(np.abs(newspec) > times_lower_thr* RMS)[0]
            newspec[strange] = np.random.normal(scale=RMS, size=len(strange))

            data_rep[tn] = newspec

    return data_rep

def replace_rfi_substract2(data,freq,time_rfi,mw_use,sg_window,sg_polyorder, s_sigma,
                          rms_sigma,rms_frange,times_low_thr,times_lower_thr = 5):
    """
    replace big RFI to reduce the impact in FFT
    
    Parameters:
    data: np array
    sg_window: scipy.signal.savgol_filter window length (unit MHz)
    sg_polyorder: scipy.signal.savgol_filter polyorder
    
    """
    fdelta = freq[1]-freq[0] 
    
    from markRFI import real_rms
    if time_rfi is not None:
        whole_rfi = np.all(time_rfi,axis = 1)
        is_rfi_num = np.arange(data.shape[0])[whole_rfi]
        not_rfi_num = np.arange(data.shape[0])[~whole_rfi]
    else:
        not_rfi_num = np.arange(data.shape[0])

    data_rep = np.full(data.shape,np.nan)
    from hifast.utils.misc import extend_Trues

    for tn in tqdm(range(data.shape[0])):
        if tn in not_rfi_num:
            spec = deepcopy(data[tn,:])
            RMS = real_rms(spec,freq,sigma=rms_sigma,rms_vrange=rms_frange)

            thr = RMS*times_low_thr
            low_use = (np.abs(spec) > thr)    
            low_use = extend_Trues(low_use,ext_add =10,leng_lim = 20,axis = -1)

            newspec = deepcopy(spec)
            newspec[mw_use] = 0
            sg_sm = filter_smooth(spec,fdelta,method = 'savgol',sg_window=sg_window,
                                  sg_polyorder=sg_polyorder)
            # period rfi and source
            newspec[low_use] = spec[low_use] - sg_sm[low_use]
            # mw 
            newspec[mw_use] = np.random.normal(scale=RMS,size = np.sum(mw_use)) 
            if np.sum(time_rfi[tn]) > 0:
                gauss_sm = filter_smooth(spec,fdelta,method ='gaussian',s_sigma = s_sigma)
                # time RFI 
                newspec[time_rfi[tn]] = (spec[time_rfi[tn]] / gauss_sm[time_rfi[tn]] - 1)
                cond = time_rfi[tn] & (np.abs(newspec) > thr * 2)
                newspec[cond] = spec[cond]
            strange = np.where(np.abs(newspec) > times_lower_thr * RMS)[0]
            newspec[strange] = np.random.normal(scale=RMS, size=len(strange))

            data_rep[tn] = newspec

    return data_rep

def repalce_near(data,freq,time_rfi,mw_use=None,times_lower_thr=None,rms_sigma = None,
                 ext_freq = None,rms_frange=None,rfi_width_lim=None, ext_sec=None):
    """
    replace mw, rfi or others by near ripple section
    
    times_lower_thr:above * times of rms will be lowered
    ext_freq: extend freq range to replace (mhz)
    rfi_width_lim:rfi should contain more channels than limit 
    ext_sec: extend channel number of start and end of each section
    """
    fdelta = freq[1]-freq[0] 
    N = len(freq)
    ext = int(np.around(ext_freq / fdelta))
    
    from markRFI import real_rms,rms,get_startend
    if time_rfi is not None:
        whole_rfi = np.all(time_rfi,axis = 1)
        is_rfi_num = np.arange(data.shape[0])[whole_rfi]
        not_rfi_num = np.arange(data.shape[0])[~whole_rfi]
    else:
        not_rfi_num = np.arange(data.shape[0])
        time_rfi = np.full(data.shape,False)

    if mw_use is None:
        mw_use = np.full(freq.shape,False)
        
    data_rep = np.full(data.shape,np.nan)
    from hifast.utils.misc import extend_Trues

    for tn in tqdm(range(data.shape[0])):
        if tn in not_rfi_num:
            spec = deepcopy(data[tn,:])
            RMS = real_rms(spec,freq,sigma=rms_sigma,rms_vrange=rms_frange)

            thr = RMS*times_lower_thr
            low_use = (np.abs(spec) > thr) | time_rfi[tn] | mw_use
            start,end = get_startend(low_use,rfi_width_lim, ext_sec)

            newspec = deepcopy(spec)
            newspec_ = deepcopy(spec)
            #usespec = np.zeros_like(spec)
            for s,e in zip(start,end):
                spec1 = np.zeros_like(freq)
                s0 = s - ext
                if s0 < 0:s0 = 0
                e0 = e + ext
                if e0 > N - 1:e0 = N -1
                spec1[s0:e0] = spec[s0:e0]

                peak = np.max(spec[s:e])
                peak_loc = np.where(spec1 == peak)[0][0]
                spec1l = deepcopy(spec1);spec1l[peak_loc:] = 20
                l1 = np.argmin(spec1l)
                spec1r = deepcopy(spec1);spec1r[:peak_loc] = 20
                r1 = np.argmin(spec1r)
                #print(l1,r1)
                L = r1 - l1
                if l1 - L < 0:
                    direc = 'right'
                elif r1 + L > N - 1:
                    direc = 'left'
                else:
                    rmsl = rms(newspec[l1 - L:l1],freq[l1 - L:l1])
                    rmsr = rms(newspec[r1:r1 + L],freq[r1:r1 + L])
                    if rmsl < rmsr:
                        direc = 'left' 
                    else:
                        direc = 'right'
                if direc == 'left':
                    s1 = l1 - L
                    e1 = l1
                elif direc == 'right':
                    s1 = r1
                    e1 = r1 + L
                #print(direc,s1,e1)

                try:
                    newspec_[l1:r1] = newspec[s1:e1]
                    newspec[s:e] = newspec_[s:e]
                    #usespec[s1:e1] = 1
                except ValueError:
                    log.info(f"tn = {tn} has a ValueError")
                    import traceback
                    traceback.print_exc()  
                    newspec[s:e] = np.random.normal(scale=RMS, size=int(e-s))

            strange = np.where(np.abs(newspec) > 3 * RMS)[0]
            newspec[strange] = np.random.normal(scale=RMS, size=len(strange))

            data_rep[tn,:] = newspec
            
    return data_rep

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
    elif method == 'set noise':
        data_low[low_use] = np.random.normal(scale=RMS, size=data.shape)[low_use]
    return data_low

def replace_rfi(data,freq,is_rfi = None,time_rfi = None, method='subtract with pd',
                mw_use = None, **rep_args):
    
    log.info(f"Replace RFI with {method} method ...")
    if method == 'near ripple':
        data_rep = repalce_near(data,freq,time_rfi=time_rfi,mw_use=mw_use,**rep_args)
    elif method == 'subtract trpdr':
        data_rep = replace_rfi_substract(data,freq,is_rfi,time_rfi,mw_use,**rep_args)
    elif method == 'subtract tr':
        data_rep = replace_rfi_substract2(data,freq,time_rfi,mw_use,**rep_args)
    elif (method == 'lower') or (method == 'set zeros') or (method == 'set noise'):
        data_rep = replace_rfi_lower(data,freq,method=method,**rep_args)
    elif method == 'subtract rfi':
        data_rep = replace_rfi_substract_old(data,freq,is_rfi,mw_use,**rep_args)
    else:
        raise ValueError("Unsupport replace RFI method!")
        
    return data_rep  

############################## fit ripple #####################################

        
def find_loc(amp,x,amp_thr,xlim,rip_mhz=True,Print=True):
    """
    amp: fft amptitude
    x: fft xfreq
    amp_thr: amp threshold
    rip_mhz: bool
    xlim: find max between xlim[0] and xlim[1]
    """
    loc = None; sw_mhz = None
    if rip_mhz:
        u = (amp>amp_thr)&(x>xlim[0])&(x<xlim[1])
        
        if np.sum(u) > 0:
            amp_ = np.zeros_like(x)
            amp_[u] = amp[u]

            loc = np.argmax(amp_)
            sw_mhz = 1/x[loc]
            if Print:
                print(f'find ripple {sw_mhz:.5f} MHz,{1/sw_mhz:.5f} mu s,locate in {loc}')
        else:
            rip_mhz = False

    return loc,sw_mhz,rip_mhz

def fft_fit_ripple(data_rep, freq,is_rfi_num,not_rfi_num,ori_shape,amp_thr_mean_factor = 1.05,
                   amp_thr_factor = 1.5,is_on = None,chan_wide = 5,chan_narr = 3,choose_method = 'all',
                   rip_base = True,rip_1mhz= True,rip_2mhz = False,rip_0_04mhz = False,
                   plot = False,pdf = None,title = None,fft_ylim = None,
                   quick_test = False, fftf = None,amp_data = None, x = None, tn = None):
    """
    fit baseline ripple (standing wave) by FFT
    Parameter:
    data_rep: data array (after replace RFI)
    rip_*mhz: bool, remove standing wave period in Fourier space?
    
    amp_thr_mean_factor: above amptitude threshold will be chosed in mean amp
    
    chan_wide: channel number next to sw_freq will be chosed all.
    chan_narrow:
    
    """
    fdelta = freq[1]-freq[0]
    
    if len(is_rfi_num) > 0:
        data_rep = np.delete(data_rep,is_rfi_num,axis = 0)
        
    is_on_num = np.arange(data_rep.shape[0])[is_on]
    #is_off_num = np.arange(data_rep.shape[0])[~is_on]
    
    # FFT
    if not quick_test:
        fftf = FFT(data_rep, freq)
        x = fftf.x
        amp_data = fftf.amp
    
    amp = np.mean(amp_data,axis=0)
    amp_on = np.mean(amp_data[is_on],axis=0)
    amp_off = np.mean(amp_data[~is_on],axis=0)
    
    amp_thr = np.median(amp) * amp_thr_mean_factor
    amp_thr_on = np.median(amp_on) * amp_thr_factor
    amp_thr_off = np.median(amp_off) * amp_thr_factor
    print(f"mean amp thr:{amp_thr:.3f}, amp thr on:{amp_thr_on:.3f}, amp thr off:{amp_thr_off:.3f}")
    
    # ripples freq in fourier space
    # 1 mhz
    loc1,sw_1mhz,rip_1mhz = find_loc(amp,x,amp_thr,rip_mhz=rip_1mhz,xlim=[.90,.95])
    loc1_,sw_1mhz_,_ = find_loc(amp,x,amp_thr,rip_mhz=rip_1mhz,xlim=[1.8,1.9])
    loc1s = np.array([loc1,loc1_])
    # 2 mhz
    loc2,sw_2mhz,rip_2mhz = find_loc(amp,x,amp_thr,rip_mhz=rip_2mhz,xlim=[.5,.6])
    # 0.04 mhz
    loc0_04,sw_0_04mhz,rip_0_04mhz = find_loc(amp,x,amp_thr,rip_mhz=rip_0_04mhz,xlim=[25,26])
    # 8 mhz rfi
    #if rfi_8mhz:
    #    x_loc16 = np.arange(0,x[-1],rfi_8mhz_step)
    #   loc16 = np.around(x_loc16 / (fs/N)).astype('int')
    
    # select which components?
    use = np.zeros_like(x,dtype = 'bool')
    use_1mhz = np.zeros_like(x,dtype = 'bool')
    use_2mhz = np.zeros_like(x,dtype = 'bool')
    use_0_04mhz = np.zeros_like(x,dtype = 'bool')
    #use_8mhz = np.zeros_like(x,dtype = 'bool')
    
    if rip_1mhz:
        for s1 in range(len(loc1s)):
            use_1mhz |= (np.abs(np.arange(len(x))-loc1s[s1]) < chan_wide)
        use[use_1mhz] = True
    if rip_2mhz:
        use_2mhz = (np.abs(np.arange(len(x))-loc2) < chan_wide)
        use[use_2mhz] = True
    if rip_0_04mhz:
        use_0_04mhz = (np.abs(np.arange(len(x))-loc0_04) < chan_narr)
        use[use_0_04mhz] = True
    #if rfi_8mhz:
    #    use_16mhz = np.full(len(x),False)
    #    use_16mhz[loc16] = True; use16[0] = False
        
    use_amp = use  & (amp_data > amp_thr_off)
    use_amp[is_on,:] = use  & (amp_data[is_on,:] > amp_thr_on)
    
    amp_data_inpd = deepcopy(amp_data)
    from scipy.interpolate import interp1d
    x_ = deepcopy(x)
    if choose_method == 'interpolate':
        # interpolate
        for ti in tqdm(range(amp_data.shape[0])):
            x_mask = x_[~use_amp[ti]]
            amp_data_mask = amp_data[ti][~use_amp[ti]]
            amp_data_interp = interp1d(x_mask,amp_data_mask,kind='linear')#,fill_value="extrapolate")
            amp_data_inpd[ti][use_amp[ti]] = amp_data_interp(x_[use_amp[ti]])  
    elif choose_method == 'median':
        raise ValueError("I don't want to write this method now...")

    amp_data_inpd = amp_data - amp_data_inpd
    if choose_method == 'all':
        # use all 
        amp_data_inpd[use_amp] = amp_data[use_amp]
    # 0 mu s constant
    if rip_base:
        amp_data_inpd[:,0] = amp_data[:,0]
    
    amp_data_inpd[amp_data_inpd < 0] = 0
    
    # IFFT
    A_data_inpd =  amp_data_inpd * np.exp(1j*fftf.phi)
    A_data_ifft = np.real(np.fft.irfft(A_data_inpd,n=data_rep.shape[1]))
    print("Finished ifft.")
    
    if plot:
        if not quick_test:
            tn = not_rfi_num[10]
        if pdf is not None:
            plt.switch_backend('agg')   
            
        if tn in is_on_num:
            amp_thr = amp_thr_on
        else:
            amp_thr = amp_thr_off
       
        fig = plt.figure(figsize=(15,6))
        ax = fig.add_subplot(211) 
        if choose_method != 'all':
            ax.plot(x,amp_data[tn]-amp_data_inpd[tn],'g',label = f'{choose_method}')
        ax.plot(x,amp_data[tn], label = 'whole fft')
        ax.plot(x,amp_data_inpd[tn],label = 'sw fft')
        ax.axhline(amp_thr,color = 'k',linestyle='--',label = f'amp thr = {amp_thr:.2f}',alpha=.5)
        if rip_1mhz:
            ax.text(x[loc1],amp_thr,f'{sw_1mhz:.4f} MHz',color = 'k',size = 16)
            ax.text(x[loc1_],amp_thr,f'{sw_1mhz_:.4f} MHz',color = 'k',size = 16)
        if rip_2mhz:
            ax.text(x[loc2],amp_thr,f'{sw_2mhz:.4f} MHz',color = 'k',size = 16)
            #ax.plot(x[use_i[tn]],np.zeros_like(x[use_i[tn]]),'r^',label = '2 mhz')
        #if rfi_8mhz:
        #    ax.plot(x[loc16],np.zeros_like(loc16),'gs',label = f'step {rfi_8mhz_step:.4f} $\mu$s ')
        if rip_base | rip_1mhz :
            ax.plot(x[use_amp[tn]],np.zeros_like(x[use_amp[tn]]),'r^',label = 'base,1mhz,2mhz')
        ax.set_xlim(-.01,2)
        if fft_ylim is not None:
            ax.set_ylim(fft_ylim[0],fft_ylim[1])
        ax.grid();ax.legend(loc = 2)
        ax.set_xlabel('k [$\mu$s]')
        if title is not None:
            ax.set_title(title)
        
        ax = fig.add_subplot(212)
        if choose_method != 'all':
            ax.plot(x,amp_data[tn]-amp_data_inpd[tn],'g',label = f'{choose_method}') 
        ax.plot(x,amp_data[tn], label = 'whole fft')
        ax.plot(x,amp_data_inpd[tn],label = 'sw fft')
        ax.axhline(amp_thr,color = 'k',linestyle='--',label = f'amp thr = {amp_thr:.2f}',alpha=.5)
        if rip_0_04mhz:
            ax.text(x[loc0_04],amp_thr,f'{sw_0_04mhz:.4f} MHz',color = 'k',size = 16)
            ax.plot(x[use_amp[tn]],np.zeros_like(x[use_amp[tn]]),'ro',label = '0.04mhz')
        #if rfi_8mhz:
        #    ax.plot(x[loc16],np.zeros_like(loc16),'gs',label = f'step {rfi_8mhz_step:.4f} $\mu$s ')
        ax.set_xlim(25,26)
        if fft_ylim is not None:
            ax.set_ylim(fft_ylim[0],fft_ylim[1])
        ax.grid();ax.legend(loc = 2)
        ax.set_xlabel('k [$\mu$s]')
        
        plt.tight_layout()
        if pdf is not None:
            pdf.savefig();plt.close()
        else:
            plt.show()

    if len(is_rfi_num)>0:
        A_data_ifft_ = np.full(ori_shape,np.nan)
        A_data_ifft_[not_rfi_num,:] = A_data_ifft
        A_data_ifft = deepcopy(A_data_ifft_)
    
    return A_data_ifft



def mean_fit_ripple(data, nspec,func = 'iter'):
    
    n = nspec // 2
    if func == 'iter':
        print(f'mean {2*n}')
        import bottleneck as bn
        tlen = data.shape[0]
        t1 = np.arange(tlen)-n
        t2 = np.arange(tlen)+n
        t1[t1 < 0] = 0;t2[t2 < 2*n] = 2*n
        t1[t1 > tlen - 2*n -1] = tlen - 2*n -1;t2[t2 > tlen -1] = tlen -1

        sw = np.zeros_like(data)
        if len(data.shape) == 2:
            for i in tqdm(range(tlen)):
                sw[i,:] = bn.nanmean(data[t1[i]:t2[i],:],axis = 0)
        elif len(data.shape) == 3:
            for i in tqdm(range(tlen)):
                sw[i,:,0] = bn.nanmean(data[t1[i]:t2[i],:,0],axis = 0)
                sw[i,:,1] = bn.nanmean(data[t1[i]:t2[i],:,1],axis = 0)
                
    elif func == 'smooth':
        print(f'mean {2*n+1}')
        from hifast.utils.misc import smooth1d
        s_method_t = 'boxcar'
        print('Smooth ing ...')
        sw = smooth1d(data,axis = 0,sigma = n, method = s_method_t)
    return sw

def med_fit_ripple(data, nspec,func = 'iter'):
    
    n = nspec // 2
    if func == 'iter':
        print(f'median {2*n}')
        import bottleneck as bn
    
        tlen = data.shape[0]
        t1 = np.arange(tlen)-n
        t2 = np.arange(tlen)+n
        t1[t1 < 0] = 0;t2[t2 < 2*n] = 2*n
        t1[t1 > tlen - 2*n -1] = tlen - 2*n -1;t2[t2 > tlen -1] = tlen -1

        sw = np.zeros_like(data)
        if len(data.shape) == 2:
            for i in tqdm(range(tlen)):
                sw[i,:] = bn.nanmedian(data[t1[i]:t2[i],:],axis = 0)
        elif len(data.shape) == 3:
            for i in tqdm(range(tlen)):
                sw[i,:,0] = bn.nanmedian(data[t1[i]:t2[i],:,0],axis = 0)
                sw[i,:,1] = bn.nanmedian(data[t1[i]:t2[i],:,1],axis = 0)
    elif func == 'smooth':
        print(f'median {2*n+1}')
        from hifast.utils.misc import smooth1d
        s_method_t = 'median'
        print('Smooth ing ...')
        sw = smooth1d(data,axis = 0,sigma = n, method = s_method_t)            
    
    return sw

def fit_ripple(data,freq=None, method='rfft', is_rfi_num=None,not_rfi_num=None,ori_shape=None,
               plot = False,pdf = None,title = None,is_on = None,**fit_args): 
    
    log.info(f"Fit baseline ripple with {method} method ...")
    if method == 'rfft':
        A_data_ifft = fft_fit_ripple(data, freq,is_rfi_num,not_rfi_num,ori_shape,
                   plot = plot,pdf = pdf,title = title,is_on = is_on, **fit_args)
        
        return A_data_ifft        
    elif method == 'mean':
        sw = mean_fit_ripple(data,**fit_args)
        return sw
    elif method == 'median':
        sw = med_fit_ripple(data,**fit_args)
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
        