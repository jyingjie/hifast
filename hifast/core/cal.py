

__all__ = ['FastRawData', 'adjust_frange', 'get_str_aft_nB', 'HiFASTData', 'FastRawSpec', 'CalOnOff', 'plot_sep',
           'mean_a', 'CalOnOffA', 'PositionSwitch']


"""
@author: yixiancao
@author: YingjieJing
"""
import os
import re
from glob import glob

import numpy as np
import scipy.interpolate as interp
import h5py
from astropy.io import fits

import __main__
if hasattr(__main__, '__file__'):
    import matplotlib
    matplotlib.use('Agg')
from matplotlib import pyplot as plt


from ..utils.tcal import read_tcal
from ..utils.misc import smooth_axis1_d3, down_sample, median_filter_axis1_d3
from ..utils.misc import smooth1d
from ..utils.io import save_specs_hdf5, MjdChanPolar_to_PolarMjdChan
from ..utils.io import H5FitsRead


class FastRawData(object):
    """ RAW data from FAST observations.

    Attributes
    ----------
    filename: str
        Name of Raw data file (FITS file)
    fileinfo: dict
        Infomation inferred from the filename
    obs: dict
        Basic observation informaiton from observation log (NOW hd0 header).

    Methods
    ----------
    listInfo():
        Display the basic infomation about the data file.

    """


    def __init__(self, filenames, frange=None, dfactor=None, med_filter_size=None, verbose=False):
        """
        Parameters
        ----------
        filenames : list
            FAST fits file path name
        frange : [min max], optional
            Frequency range
        dfactor : str 'W' or int or None, optional
            if 'W', downsampling 'N' or 'F' to 'W' band; if int, downsampling the spectra with a factor of dfactor
        med_filter_size : int, optional
            median filter spectra with size of this value
        verbose : bool, optional
            Default is False
        """
        self.filenames = filenames
        self.verbose = verbose
        if self.verbose:
            print('files:')
            for filename in self.filenames:
                print(filename)
        if self.filenames[0].endswith('.fits'):
            self.hduls = [fits.open(filename) for filename in self.filenames]
            self.ftype = 'fits'
        elif self.filenames[0].endswith('.hdf5'):
            self.hduls = [H5FitsRead(filename) for filename in self.filenames]
            self.ftype = 'hdf5'
        else:
            raise(ValueError('File not support'))

        self.hd0s = [hdul[0].header for hdul in self.hduls]
        self.hd1s = [hdul[1].header for hdul in self.hduls]
        self.lens= [header['NAXIS2'] for header in self.hd1s]
        self.inds = np.arange(np.sum(self.lens))

        if isinstance(dfactor, str):
            if dfactor.upper() == 'W':
                dfactor = 16
            elif dfactor.isdigit():
                dfactor = int(dfactor)
            else:
                raise(ValueError("dfactor should be 'W' or an int"))
        # adjust frange
        if dfactor == 16 and frange is not None:
            frange = adjust_frange(frange)

        self.frange = frange
        self.dfactor = dfactor
        self.med_filter_size = med_filter_size
        self._get_freq()
        if frange is not None:
            is_use = (self.freq<= frange[1] ) & (self.freq>= frange[0])
            self.freq_use = self.freq[is_use]
        else:
            self.freq_use= self.freq
        if len(self.freq_use) == 0:
            raise(ValueError(f"frange {frange} is not in freq range of {self.freq}"))
        if self.dfactor is not None:
            self.freq_use = down_sample(self.freq_use[None,:,None], self.dfactor)[0,:,0]

    def ls_fields(self,):
        """
        list fields in fits hdu1.data
        """
        print(self.hduls[0][1].data.dtype)

    def get_mjds(self):
        mjds = np.hstack([hdul[1].data['UTOBS'] for hdul in self.hduls])
        return mjds

    def _get_freq(self, center_corr=None):
        tdata = [self.hduls[0][1].data]
        nchan = tdata[0]['NCHAN'][0]
        freq0 = tdata[0]['FREQ'][0]
        chanwidth = tdata[0]['CHAN_BW'][0]
        freq=(np.arange(nchan))*chanwidth
        freq=freq+freq0
        if center_corr is not None:
            freq += center_corr # using center frequency
        self.freq= freq

    def get_field(self, _inds, field='DATA', close_file=False):
        """"
        get field with index
        """
        # if input _inds is empty, return empty
        if 0 in _inds.shape:
            if field == 'DATA':
                return np.empty(_inds.shape + (len(self.freq_use),2))
        # flat _inds
        _input_ndim = _inds.ndim
        if _input_ndim > 1:
            _in_shape= _inds.shape
            _inds= _inds.flatten()
        # h5py needs flattened _inds is monotone increasing
        if self.ftype == 'hdf5':
            if (np.diff(_inds) <= 0).any():
                # load all data to the numpy array first
                load_all_first = True
            else:
                load_all_first = False
        else:
            load_all_first = False
        # determine the file and index of the spec in
        lens_cum = np.hstack([0,np.cumsum(self.lens)])
        ifile =  np.searchsorted(lens_cum, _inds + 1, side='left') - 1
        ind_ifile = _inds- lens_cum[ifile]
        ifile_uni, ifile_num = np.unique(ifile, return_counts=True)
        ind_ifile_list= np.split(ind_ifile, np.cumsum(ifile_num)[:-1])
        # load data
        j_list=[]
        for j in ifile_uni:
            try:
                # test if file is opened
                if self.ftype == 'fits':
                    self.hduls[j][1].data.shape
                elif self.ftype == 'hdf5':
                    self.hduls[j][1].data[field].shape
            except:
                if self.ftype == 'fits':
                    self.hduls[j] = fits.open(self.filenames[j], memmap=True, lazy_load_hdus=True)
                elif self.ftype == 'hdf5':
                    self.hduls[j] = H5FitsRead(self.filenames[j])
                j_list += [j,]
        if field == 'DATA':
            if self.frange is not None:
                is_use = (self.freq <= self.frange[1] ) & (self.freq >= self.frange[0]) # can't use self.freq_use
                ind_use= np.where(is_use)[0]
                # freq axis (ind_use) need use "slice" to index, coz ii is already a array; ind_use is continuous
                if not load_all_first:
                    data = np.vstack([self.hduls[i][1].data[field][ii, ind_use[0]:ind_use[-1]+1, :2] for i, ii in zip(ifile_uni,ind_ifile_list)])
                else:
                    data = np.vstack([self.hduls[i][1].data[field][:][ii, ind_use[0]:ind_use[-1]+1, :2] for i, ii in zip(ifile_uni,ind_ifile_list)])
            else:
                if not load_all_first:
                    data = np.vstack([self.hduls[i][1].data[field][ii, :, :2] for i, ii in zip(ifile_uni,ind_ifile_list)])
                else:
                    data = np.vstack([self.hduls[i][1].data[field][:][ii, :, :2] for i, ii in zip(ifile_uni,ind_ifile_list)])
            #
            if self.med_filter_size is not None:
                data = median_filter_axis1_d3(data, self.med_filter_size)
            # Downsampling
            if self.dfactor is not None:
                data = down_sample(data.astype('float64'), self.dfactor)
        else:
            data = np.hstack([self.hduls[i][1].data[field][ii] for i, ii in zip(ifile_uni,ind_ifile_list)])
        #check if the length in the header is correct
        for i_tmp in ifile_uni:
            if self.ftype == 'fits':
                  leng = self.hduls[i_tmp][1].data.shape[0]
            elif self.ftype == 'hdf5':
                  leng = self.hduls[i_tmp][1].data[field].shape[0]
            if not self.lens[i_tmp] == leng:
                raise(ValueError(f"length in the header is not equal to the data shape:", self.filenames[i_tmp]))
        [self.hduls[j].close() for j in j_list]
        if close_file:
            [self.hduls[i].close() for i in ifile_uni]

        if _input_ndim > 1:
            return data.reshape(_in_shape + data.shape[1:])
        else:
            return data

    @property
    def fileinfo(self):
        fitsname = dict(value = self.filenames[0].split('/')[-1], desc = "First FITS file name")
        dummy = fitsname['value']
        tinfo = dummy.split('-')[0]
        tinfo = tinfo.split('_')
        #project = dict(value = 'FAST_M31', desc = "Project name")
        source = dict(value = tinfo[0], desc = "Source name")
        ipos = dict(value  = int(tinfo[1][1:]), desc = "Position number")
        #iobs = dict(value = None, desc = "Observation number")
        track_mode = dict(value = tinfo[2], desc = "Observation mode") # e.g. Tracking, drifting

        finfo = dummy.split('-')[1]
        finfo = finfo.split('_')
        ibeam = dict(value = int(finfo[0][1:]), desc = "Beam number")
        #ifile = dict(value = int(finfo[-1][0:4]), desc = "File number") #ith file for the observation

        filetype = dict(value = 'psr' if len(finfo) < 3 else 'spec' +  finfo[1],
                         desc = "Observation type")

#         fileinfo = dict(fitsname = fitsname, project=project, source = source, ipos = ipos, track_mode = track_mode,
#                         iobs =iobs, ibeam = ibeam, ifile = ifile, filetype = filetype)
        fileinfo = dict(fitsname = fitsname, source = source, ipos = ipos, track_mode = track_mode,
                         ibeam = ibeam, filetype = filetype)

        return fileinfo

    def listInfo(self, list_type = 'all'):
        finfo = self.fileinfo
        print("FAST RAW DATA FILE Summary")
        print("==================")
        for key, value in finfo.items():
            print (value['desc'], ':' , value['value'])

    #@property
    def obs(self,chunk_num=0):
        """ Dictionary of observation parameters.
            From HD0 or from observation log
        """
        hd0 = self.hd0s[chunk_num]
        time = {"t_obs": dict(value  = hd0['DATE'],
                              desc = 'Observation time (starting to record data)') ,
                "t_samp": dict(value = None,
                               desc = 'Sampling time interval (s)')
                }
        track = {"Mode": dict(value = self.fileinfo['track_mode'],
                                    desc =  "Tracking mode (e.g. tracking, drifting, etc)"),
                "RA": dict(value = np.double(0.0), desc = "Tracking RA"),
                "DEC": dict(value = np.double(0.0), desc = "Tracking DEC")}

        noise_diode = {"Mode": dict(value = "",  desc = "Noise diode mode (Modulate, ON, or OFF)"),
                      "Power": dict(value = "",  desc = "Noise diode power (High/Low)"),
                      "Delay": dict(value = 0, desc =  "Noise Delay in unit of tunit"),
                      "On": dict(value = 0,  desc = "Noise diode ON in unit of tunit"),
                      "Off":  dict(value = 0,  desc = "Noise diode OFF in unit of tunit"),
                      "tunit":  dict(value = 4e-9, desc =  "Units of delay/on/off time (s)")
                      } #
        receiver = {"Frontend": dict(value = '', desc =  "Receiver frontend (e.g. 19BEAM)"),
                    "rfgain": dict(value = 0.0, desc =  "Receiver rfgain (dB)"),
                    "dgain": dict(value = 0.0, desc = "Receiver dgain"),
                    "Backend": dict(value = '', desc = "Receiver Backend ID (e.g. 'MB4K')")
                    }

        obs = dict(time = time,
               track = track,
               noise_diode = noise_diode,
               receiver = receiver)

        return obs


def adjust_frange(frange):
    """
    adjust freq range frange to let the F band downsample freq is same with W band
    """

#     # W band
#     nchan = 65536
#     chanwidth = 0.00762939453125
#     freq0 = 1000.0035762786865
#     freq_W = np.arange(nchan)*chanwidth + freq0 + 0.000476837158203125/2.

    # F band
    nchan = np.array(1048576, dtype='int64')
    chanwidth = np.array(0.000476837158203125, dtype='float64')
    freq0 = np.array(1000.0, dtype='float64')

    chanwidth_F = np.array(0.000476837158203125, dtype='float64')
    freq_F = np.arange(nchan)*chanwidth + freq0 + chanwidth_F/2.

    offset = chanwidth_F/4.
    ind_close = np.argmin(abs(freq_F - np.array(frange)[:,None]), axis=1)
    ind_close = ((ind_close+1)//16+1)*16-1
    frange_new = (freq_F[ind_close[0]] + offset, frange[1])
#     width = int(np.round(abs(np.log10(offset)))+1)
#     print(ind_close)
#     print(freq_F)
#     print(f"%.{width}f %.{width}f" % frange_new)
    return frange_new


def get_str_aft_nB(path, nth=1):
    """
    path
    nth: beginning with 1
    """
    basename = os.path.basename(path)
    ind = re.search(r'[0-9][0-1]M-', basename[::-1]).start()
    str_list = basename[-ind:].split('-')[1:]
    if len(str_list) >= nth:
        return str_list[nth-1]
    else:
        return None


class HiFASTData(FastRawData):
    """ data in HiFAST.

    """


    def __init__(self, filenames, frange=None, dfactor=None, med_filter_size=None, verbose=False):
        """
        Parameters
        ----------
        filenames : list
            FAST fits file path name
        frange : [min max], optional
            Frequency range
        dfactor : str 'W' or int or None, optional
            if 'W', downsampling 'N' or 'F' to 'W' band; if int, downsampling the spectra with a factor of dfactor
        med_filter_size : int, optional
            median filter spectra with size of this value
        verbose : bool, optional
            Default is False
        """
        self.filenames = filenames
        self.verbose = verbose
        if self.verbose:
            print('files:')
            for filename in self.filenames:
                print(filename)

        self.ftype = 'hdf5' # must be 'hdf5'
        self.hduls = [H5FitsRead(filename, kT='Power', hdu_names=['S'], DATA_key='Power') for filename in self.filenames]
        for hdul in self.hduls:
            hdul.hdus[1] = hdul.hdus['S']
            hdul.hdus[1].name = '1'

        self.hd0s = [{} for hdul in self.hduls]
        self.hd1s = [hdul[1].header for hdul in self.hduls]
        self.lens= [hdul[1].data['mjd'].shape[0] for hdul in self.hduls]
        self.inds = np.arange(np.sum(self.lens))

        if isinstance(dfactor, str):
            if dfactor.upper() == 'W':
                dfactor = 16
            elif dfactor.isdigit():
                dfactor = int(dfactor)
            else:
                raise(ValueError("dfactor should be 'W' or an int"))
        # adjust frange
        if dfactor == 16 and frange is not None:
            frange = adjust_frange(frange)

        self.frange = frange
        self.dfactor = dfactor
        self.med_filter_size = med_filter_size
        self._get_freq()
        if frange is not None:
            is_use = (self.freq<= frange[1] ) & (self.freq>= frange[0])
            self.freq_use = self.freq[is_use]
        else:
            self.freq_use= self.freq
        if len(self.freq_use) == 0:
            raise(ValueError(f"frange {frange} is not in freq range of {self.freq}"))
        if self.dfactor is not None:
            self.freq_use = down_sample(self.freq_use[None,:,None], self.dfactor)[0,:,0]

    def ls_fields(self,):
        """
        list fields in fits hdu1.data
        """
        print(self.hduls[0][1].data.keys())

    def get_mjds(self):
        mjds = np.hstack([hdul[1].data['mjd'] for hdul in self.hduls])
        return mjds

    def _get_freq(self, center_corr=None):
        # center_corr already in freq
        self.freq = self.get_field(field='freq')

    def get_field(self, _inds=slice(None), field='DATA', close_file=False):
        """
        """
        # if input _inds is empty, return empty
        if (not isinstance(_inds, slice)) and 0 in _inds.shape:
            if field == 'DATA':
                # if input _inds is empty, return empty
                return np.empty(_inds.shape + (len(self.freq_use),2))
        # field used get_field to read
        if field in ['DATA', 'ra', 'dec', 'mjd', 'UTOBS']:
            if field == 'UTOBS':
                field = 'mjd'
            res = super().get_field(_inds, field=field, close_file=close_file)
        elif field in ['freq']:
            # no need _inds
            res = np.hstack([hdul[1].data[field][:] for hdul in self.hduls])
        else:
            raise(NotImplementedError)
        return res

    def fileinfo(self):
        raise(NotImplementedError)

    def listInfo(self, list_type = 'all'):
        raise(NotImplementedError)

    def obs(self,chunk_num=0):
        raise(NotImplementedError)


class FastRawSpec(object):
    def __new__(cls, *args, **kwargs):
        if 'fname_part' in kwargs.keys():
            fname_part = kwargs['fname_part']
        else:
            fname_part = args[0]
        if fname_part.endswith('.hdf5') and str(get_str_aft_nB(fname_part, 2)).startswith('specs'):
            cls_base = HiFASTData
        else:
            cls_base = FastRawData
        class DynamicC(cls, cls_base):
            pass
        return super().__new__(DynamicC)

    def __init__(self, fname_part, start=1, stop=None,
                 frange=None,
                 dfactor=None, med_filter_size=None,
                 verbose=False,
                 smooth='gaussian', s_para={'s_sigma':5},
                 noise_mode=None, noise_date='auto'):
        """
        fname_part: str; one of chunk file name; '*.hdf5' or '*.fits'
        start: int; chunk file start
        stop: int; chunk file stop
        frange: frequency range
        dfactor: 'W' or int; if 'W', downsampling 'N' or 'F' to 'W' band; if int, downsampling the spectra with a factor of dfactor
        """
        # generate fits filename list
        if fname_part.endswith('.fits'):
            ftype = 'fits'
        elif fname_part.endswith('.hdf5'):
            ftype = 'hdf5'
        else:
            raise(ValueError('File not support'))
        if ftype == 'hdf5' and str(get_str_aft_nB(fname_part, 2)).startswith('specs'):
            filenames = [fname_part,]
            fname_part = fname_part.strip('.hdf5')
        else:
            fname_part = re.sub('[0-9]{4}\.%s\Z' % ftype, '', fname_part)
            if stop is not None:
                if stop < 1:
                    raise(ValueError('input stop < 1'))
            else:
                stop = len(glob(fname_part + f'[0-9][0-9][0-9][0-9].{ftype}'))
            filenames = [f"{fname_part}{i:04d}.{ftype}" for i in range(start,stop+1)]
        if len(filenames) == 0:
            raise(OSError(f"can not find file, please check fname_part:{fname_part}"))
        self.fname_part = fname_part
        super().__init__(filenames, frange=frange, dfactor=dfactor, med_filter_size=med_filter_size, verbose=verbose)

        self.nB= int(re.findall(r'-M[0-1][0-9]', fname_part)[-1][2:])
        self.smooth = smooth
        self.s_para = s_para
        self.noise_mode = noise_mode
        self.noise_date = noise_date

    def _get_freq(self,):
        # only for W,F,N band
        #freq += 0.000476837158203125/2 # using center frequency
        super()._get_freq(center_corr = 0.000476837158203125/2)

    def _get_smoothed(self, power, freq=None, use_ndimage=False, check_nan=False):
        """
        power: ndim 3
        freq: if is None, use self.freq_use
        """
        smooth = self.smooth
        if freq is None:
            freq = self.freq_use
        s_para = self.s_para
        power = power.astype('float64')
        if smooth=='gaussian':
            sigma = s_para['s_sigma']/(np.nanmax(freq) - np.nanmin(freq))*len(freq)
            if use_ndimage:
                from scipy import ndimage
                s_power = ndimage.gaussian_filter1d(power, sigma=sigma, axis=1)
            elif check_nan:
                is_ = np.isfinite(power)
                power[~is_] = 0
                s_power = smooth_axis1_d3(power, method=smooth, sigma=sigma)
                s_power /= smooth_axis1_d3(is_.astype('float64'), method=smooth, sigma=sigma)
                s_power[np.isinf(s_power)] = np.nan
                s_power[~is_] = np.nan
                power[~is_] = np.nan
            else:
                s_power = smooth_axis1_d3(power, method=smooth, sigma=sigma)

        elif smooth=='poly':
            s_power = smooth_axis1_d3(power, method=smooth, x=freq, deg= s_para['s_deg'])
        else:
            raise(ValueError('smooth method not support'))
        return s_power

    def get_Tcal_s(self,):
        """
        load T cal from noise file, smoothed
        """
        if self.noise_date == 'auto':
            _mjd = self.get_field(np.array([0,]), field='UTOBS')[0]
        else:
            _mjd = None
        tc_freq, tc_, self.tcal_file = read_tcal(self.nB, mode= self.noise_mode, date=self.noise_date, mjd=_mjd)
        if self.frange is not None:
            is_use = (tc_freq <= self.frange[1] ) & (tc_freq >= self.frange[0])
        else:
            is_use = np.full(len(tc_freq), True, dtype=bool)
        tc_= tc_[:, is_use]
        tc_freq= tc_freq[is_use]
        tc= self._get_smoothed(tc_.T[None,:,:], freq=tc_freq)
        tc_inter = interp.interp1d(tc_freq, tc, axis=1, kind='linear', fill_value ='extrapolate')(self.freq_use)
        return tc_inter

    def gen_out_name_base(self, outdir):
        fname_add = os.path.basename(os.path.dirname(os.path.abspath(self.fname_part)))
        self.out_name_base = os.path.join(outdir, f"{os.path.basename(self.fname_part)[:-1]}-{fname_add}")


class CalOnOff(FastRawSpec):
    """
    """
    def __init__(self, fname_part, n_delay, n_on, n_off, med_filter_size_cal=5, p_cal_fname=None, **kwargs):
        """

        Parameters
        ----------
        fname_part : str
            file name with out chunk postfix, for example:
            '/data/inspur_disk06/fast_data/3047/M31_Halo_Drift/20200106/M31_Halo_Drift_1_arcdrift-M02_F_'
        start : int
            first chunk number; Tcal on should be in the beginning of this file.
        stop : int
            last chunk number
        n_on : int
            Tcal_on Time divided by Sampling Time
        n_off : int
            Tcal_off Time divided by Sampling Time
        frange : list; [min,max]
            range of freq to use

        """
        super().__init__(fname_part, **kwargs)
        self.inds = np.arange(np.sum(self.lens))
        self.n_on, self.n_off= n_on, n_off
        self.n_delay = n_delay
        if n_delay >= 1 and self.n_on >=2:
            self.tcal_offset = True
        else:
            self.tcal_offset = False
        self.med_filter_size_cal = med_filter_size_cal
        if p_cal_fname is not None:
            self.p_cal_f = h5py.File(p_cal_fname, 'r')
        else:
            self.p_cal_f = None
        self.sep_on_off_inds()
        # set init plot as False
        # self.plot = False

    def get_extra(self,):
        """
        get extra info
        """
        extra = {}

        close_cal= np.full(self.inds.shape, False, dtype=bool)
        ind_tmp= self.inds_ton[:,0]-1
        if ind_tmp[0] ==-1:
            close_cal[ind_tmp[1:]] = True
        else:
            close_cal[ind_tmp] = True
        ind_tmp= self.inds_ton[:,-1]+1
        if ind_tmp[-1] > len(close_cal)-1:
            close_cal[ind_tmp[:-1]] = True
        else:
            close_cal[ind_tmp] = True
        extra['next_to_cal'] = close_cal

        is_on = np.full(self.inds.shape, False, dtype=bool)
        is_on[self.inds_on] = True
        extra['is_on'] = is_on

        is_delay = np.full(self.inds.shape, False, dtype=bool)
        is_delay[self.inds_delay] = True
        extra['is_delay'] = is_delay

        return extra

    def sep_on_off_inds(self):
        """
        separate on and off spec inds
        """
        n_delay, n_on, n_off = self.n_delay, self.n_on, self.n_off
        inds = self.inds
        nth = (inds - n_delay) % (n_on+n_off)
        self.inds_on = self.inds[(nth<n_on) & (inds>=n_delay)]
        self.inds_off = self.inds[(nth>=n_on) | (inds<n_delay)] # delay as off
        self.inds_delay = self.inds[inds<n_delay]
        #
        self._get_tcal_inds()

    def _get_tcal_inds(self):
        """
        find index of the samples used to calculate power of Noise diode
        """
        # drop incomplete tcal on
        n_remainder = len(self.inds_on)%(self.n_on)
        # Tcal_on index
        inds_ton= self.inds_on[:(len(self.inds_on)-n_remainder)].reshape((-1,self.n_on))
        # Tcal_off index
        nbef = self.n_on//2
        inds_toff_bef = inds_ton[:,:nbef] - nbef
        inds_toff_aft = inds_ton[:,nbef:] + self.n_on - nbef
        # fix first
        is_exceed = inds_toff_bef[0] < 0
        n_exceed = np.sum(is_exceed)
        if 0 < n_exceed < len(is_exceed):
            inds_toff_bef[0][is_exceed] = np.full(np.sum(is_exceed), inds_toff_bef[0][~is_exceed][0])
        if 0 < n_exceed and n_exceed == len(is_exceed):
            if self.tcal_offset:
                # drop first
                inds_ton = inds_ton[1:]
                inds_toff_bef = inds_toff_bef[1:]
                inds_toff_aft = inds_toff_aft[1:]
            else:
                inds_toff_bef[0][is_exceed] = inds_ton[0][-1] + 1
        # fix last
        is_exceed = inds_toff_aft[-1] > self.inds[-1]
        n_exceed = np.sum(is_exceed)
        if 0 < n_exceed < len(is_exceed):
            inds_toff_aft[-1][is_exceed] = np.full(np.sum(is_exceed), inds_toff_aft[-1][~is_exceed][-1])
        if  0 < n_exceed and n_exceed == len(is_exceed):
            if self.tcal_offset:
                # drop last
                inds_ton = inds_ton[:-1]
                inds_toff_bef = inds_toff_bef[:-1]
                inds_toff_aft = inds_toff_aft[:-1]
            else:
                inds_toff_aft[-1][is_exceed] = inds_ton[-1][0] - 1

        self.inds_ton = inds_ton
        self.inds_toff_bef = inds_toff_bef
        self.inds_toff_aft = inds_toff_aft

    def _smooth_cal(self, cal):
        if self.med_filter_size_cal is not None and self.med_filter_size_cal !=0:
            cal = median_filter_axis1_d3(cal, self.med_filter_size_cal)
        cal = self._get_smoothed(cal)
        return cal

    def _get_cal_power_s(self, inds_ton, inds_toff_bef, inds_toff_aft):
        """
        use self._get_cal_power to get the power of cal and smooth it
        """
        p_cal = self._get_cal_power(inds_ton, inds_toff_bef, inds_toff_aft)
        return self._smooth_cal(p_cal)

    def _get_cal_power(self, inds_ton, inds_toff_bef, inds_toff_aft):
        """
        get the power of Noise diode
        """
        p_ton = self.get_field(inds_ton,field='DATA')
        p_toff_aft = self.get_field(inds_toff_aft, field='DATA')
        p_toff_bef = self.get_field(inds_toff_bef, field='DATA')
        try:
            if getattr(self, 'plot', False):
                new_shape= (-1,)+p_ton.shape[-2:]
                figname = self.out_name_base + "-cal.pdf"
                plot_sep(inds_ton.flatten(), np.hstack([inds_toff_bef.reshape(-1),inds_toff_aft.reshape(-1)]),
                         p_ton.reshape(new_shape),
                         np.vstack([p_toff_bef.reshape(new_shape),p_toff_aft.reshape(new_shape)]),
                         figname=figname, n_max=1000, re_tick=True)
        except:
            pass

        if self.tcal_offset and inds_ton.shape[1] <4:
            n = p_ton.shape[-2]
            p_ton = p_ton[range(len(p_ton)), np.argmax(np.mean(p_ton[:,:,n//20:n-n//20,:], axis=(2,3), dtype='float64'), axis=1)]
            p_toff = np.concatenate([p_toff_bef, p_toff_aft],axis=1) # meger p_toff; not np.hstack or np.vstack or np.stack
            n = p_toff.shape[-2]
            p_toff = p_toff[range(len(p_toff)), np.argmin(np.mean(p_toff[:,:,n//20:n-n//20,:], axis=(2,3), dtype='float64'), axis=1)]
        if self.tcal_offset and inds_ton.shape[1] >=4:
            p_ton = np.mean(p_ton[:,1:-1,:,:], axis= 1, dtype='float64')
            p_toff = np.concatenate([p_toff_bef[:,:-1,:,:], p_toff_aft[:,1:,:,:]],axis=1) # meger p_toff
            p_toff = np.mean(p_toff, axis= 1, dtype='float64')
        if not self.tcal_offset:
            p_ton = np.mean(p_ton, axis= 1, dtype='float64')
            p_toff = np.concatenate([p_toff_bef, p_toff_aft],axis=1) # meger p_toff
            p_toff = np.mean(p_toff, axis= 1, dtype='float64')
        p_cal = p_ton.astype('float64') - p_toff.astype('float64')
        return p_cal

    def get_count_tcal(self, inds_on, inds_off):
        """
        get count of tcal
        """
        if self.p_cal_f is not None:
            inds_ton = self.p_cal_f['inds_ton'][:]
        else:
            inds_ton, inds_toff_bef, inds_toff_aft = self.inds_ton, self.inds_toff_bef, self.inds_toff_aft
        # determine which cal power sample to calibrate
        tcal_c= inds_ton[:, inds_ton.shape[1]//2]
        inds_in_tcal_on= np.argmin(abs(inds_on[:,None]- tcal_c[None,:]), axis=1)
        inds_in_tcal_off= np.argmin(abs(inds_off[:,None]- tcal_c[None,:]), axis=1)
        uni= np.unique(np.hstack([inds_in_tcal_on,inds_in_tcal_off]))
        # load cal power
        if self.p_cal_f is not None:
            p_cal_s = self.p_cal_f['p_cal_s'][uni]
        else:
            p_cal_s = self._get_cal_power_s(inds_ton[uni], inds_toff_bef[uni], inds_toff_aft[uni])
        # load on and off power
        p_on = self.get_field(inds_on, 'DATA',)
        p_off = self.get_field(inds_off, 'DATA', close_file=True)
        #print(p_on.dtype, p_off.dtype, p_cal_s.dtype)
        try:
            if getattr(self, 'plot', False):
                figname = self.out_name_base + "-sep.pdf"
                plot_sep(inds_on, inds_off, p_on, p_off, figname=figname)
        except:
            pass

        count_on = p_on.astype('float64') / p_cal_s[np.where(inds_in_tcal_on[:,None] - uni[None,:] ==0, )[1]] - 1 # have subtracted cal
        count_off = p_off.astype('float64') / p_cal_s[np.where(inds_in_tcal_off[:,None] - uni[None,:] ==0, )[1]]
        return count_on, count_off, p_cal_s, inds_ton[uni]

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
        outfield = 'Ta' if cali else 'Power'
        self.gen_out_name_base(outdir)
        n_step= self.lens[0]*step if step is not None else len(self.inds)
        inds_range= np.append(np.arange(0, self.inds[-1], n_step), self.inds[-1]+1)
        # inds_splited= np.array_split(self.inds, len(self.inds)//n_step)

        Ts=[]
        mjds = []
        tc_inter = self.get_Tcal_s() if cali else 'none'
        extra = self.get_extra()
        p_cal_s_list = []
        inds_ton_list = []

        self.plot = True
        for i in range(len(inds_range)-1):
            print('part', i)
            b, e = inds_range[i:i+2]
            inds_on = self.inds_on[(self.inds_on >=b) & (self.inds_on < e)]
            inds_off = self.inds_off[((self.inds_off >=b) & (self.inds_off < e))]
            inds = np.hstack([inds_on, inds_off])
            sort = np.argsort(inds)
            mjd = self.get_field(inds[sort], field='UTOBS') # have sorted

            if cali:
                count_tcal_res = self.get_count_tcal(inds_on, inds_off)
            else:
                count_tcal_res = [self.get_field(inds_on, 'DATA',), self.get_field(inds_off, 'DATA', close_file=True)]
            T = np.vstack(count_tcal_res[:2])
            if save_pcals and cali:
                p_cal_s_list += [count_tcal_res[2]]
                inds_ton_list += [count_tcal_res[3]]
            else:
                del count_tcal_res
            if cali: T = T*tc_inter
            T = T[sort] #sort T
            T = T.astype('float32') # finally convert to float32

            # update header
            if i==0:
                try:
                    header = {} if header is None else header.copy()
                    header.update({'tcal_file': self.tcal_file})
                except:
                    pass
            if sep_save:
                res={}
                for key in extra.keys():
                    res[key] = extra[key][inds[sort]]
                res['mjd'] = mjd
                res['freq'] = self.freq_use
                res[outfield] = MjdChanPolar_to_PolarMjdChan(T)
                res['Tcal'] = tc_inter
                #print(res)
                outname= self.out_name_base + f"-specs_T_{i:04d}_{i+1:04d}.hdf5"
                res['Header'] = header
                save_specs_hdf5(outname, res, wcs_data_name=outfield)
                print(f"Saved to {outname}")
                del res
            else:
                Ts += [T,]
                mjds += [mjd]
            self.plot = False #only plot for the first loop
        if not sep_save:
            res={}
            res.update(extra)
            res['mjd'] = np.hstack(mjds)
            res['freq'] = self.freq_use
            res[outfield] = MjdChanPolar_to_PolarMjdChan(np.vstack(Ts))
            res['Tcal'] = tc_inter
            outname= self.out_name_base + f"-specs_T.hdf5"
            res['Header'] = header
            save_specs_hdf5(outname, res, wcs_data_name=outfield)
            print(f"Saved to {outname}")
        if save_pcals and cali:
            p_cal_s_res = {}
            inds_ton = np.vstack(inds_ton_list)
            _, ind_uni = np.unique(inds_ton[:,0],return_index=True)
            p_cal_s_res['inds_ton'] = inds_ton[ind_uni]
            p_cal_s_res['p_cal_s'] = np.vstack(p_cal_s_list)[ind_uni]
            outname = self.out_name_base +'_p_cal_s.hdf5'
            print(f"p_cal in {outname}")
            p_cal_s_res['Header'] = header
            save_specs_hdf5(outname, p_cal_s_res)


def plot_sep(inds_on, inds_off, val_on, val_off, axs=None, figname=None, n_max=1500, re_tick=None, vlines_sep=True):
    if axs is None:
        fig, axs= plt.subplots(3,1,figsize=(15,4*3))
    len_freq= val_off.shape[1]
    wbin= np.min([200,len_freq//4])
    # selecting points in "n_max"
    inds_range= np.sort(np.hstack([inds_on, inds_off]))[:n_max][[0,-1]]
    is_use_on = (inds_on >= inds_range[0]) & (inds_on <= inds_range[1])
    is_use_off = (inds_off >= inds_range[0]) & (inds_off <= inds_range[1])
    x1 = inds_on[is_use_on]
    x2 = inds_off[is_use_off]
    # stacking inds as x
    x = np.hstack([x1,x2])
    ind_sort_x = np.argsort(x)
    x_sort = x[ind_sort_x]
    if re_tick is not None:
        is_1_sort = np.hstack([np.full(len(x1),True),np.full(len(x2),False)])[ind_sort_x]
        x_r = np.arange(len(x))
    for ind_s, ax in zip((len_freq-wbin)//4*np.array([1,2,3]), axs):
        y1 = np.median(val_on[is_use_on, ind_s:ind_s+wbin, 0], axis=1)
        y2 = np.median(val_off[is_use_off, ind_s:ind_s+wbin, 0], axis=1)
        y= np.hstack([y1,y2])
        if re_tick is not None:
            y_sort = y[ind_sort_x]
            ax.scatter(x_r[is_1_sort], y_sort[is_1_sort], color= 'r', label='on', marker='.', s=4)
            ax.scatter(x_r[~is_1_sort], y_sort[~is_1_sort], color= 'b', label='off', marker='.', s=4)
            #ax.plot(x_r, y_sort, 'k', lw=1)
            ax.minorticks_on()
            ax.tick_params(axis='x', which='minor', bottom=False)
            if re_tick == 'full':
                tick = x_r
                ax.set_xticks(tick)
                ax.set_xticklabels(x_sort[x_r], rotation='vertical')
            if vlines_sep:
                lines = x_r
                lines = lines[np.where(np.diff(x_sort) !=1)[0]]
                ax.vlines(lines+0.5, *ax.get_ylim(), linewidth=1, alpha=0.7)

        else:
            ax.scatter(x1, y1, color= 'r',label='on', marker='.', s=4)
            ax.scatter(x2, y2, color= 'b',label='off', marker='.', s=4)
            ax.plot(x[ind_sort_x],y[ind_sort_x],'k',lw=1)
        ax.set_ylabel('xx')
        ax.grid()
    ax.legend(frameon = True)
    ax.set_xlabel('Index')
    if figname is not None:
        fig.savefig(figname, bbox_inches='tight')


def mean_a(arr):
    """
    calculate the mean value of 5th to 95th peercent of the input arr
    """
    b, e = np.percentile(arr, [5,95],)
    return np.mean(arr[(arr <=e ) & (arr >= b)], dtype='float64')


class CalOnOffA(CalOnOff):
    """
    """
    def __init__(self, *args, pcal_vary_frac=0.01, **kwargs):
        self.pcal_vary_frac = pcal_vary_frac
        super().__init__(*args, **kwargs)

    def squeeze_freq(self, arr, axis=-1, method='mean'):
        """
        """
        if method == 'mean':
            return np.mean(arr, axis=axis, dtype='float64')
        elif method == 'mean_a':
            return np.apply_along_axis(mean_a, axis, arr)
        elif method == 'median':
            return np.median(arr.astype('float64'), axis=axis)

    def _check_diff(self, inds):
        """
        check if max diff is not larger than self.pcal_vary_frac
        """
        p = self.get_field(inds)
        ps = self.squeeze_freq(p, axis=2, method='mean_a')
        ps_max, ps_min = np.max(ps, axis=1), np.min(ps, axis=1)
        df = (ps_max - ps_min)/ps_min
        is_use = df < self.pcal_vary_frac
        return is_use[:,0] & is_use[:,1]

    def check_cal(self, inds_ton):
        """
        check cal
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

        # check off
        is_use = self._check_diff(inds_coff)
        # check on
        if inds_ton.shape[1] >= 2:
            is_use &= self._check_diff(inds_ton)
        if self.verbose:
            aba_inds_ton = inds_ton[~is_use]
            if len(aba_inds_ton) > 0:
                print("abandoned Pcal index:")
                print(*aba_inds_ton[:, 0], sep=',')
        return is_use

    def do_check_cal(self, step=None):
        is_use = self.check_cal(self.inds_ton)
        self.inds_ton = self.inds_ton[is_use]
        self.inds_toff_bef = self.inds_toff_bef[is_use]
        self.inds_toff_aft = self.inds_toff_aft[is_use]

    def get_count_tcal(self, inds_on, inds_off):
        """
        get count of tcal
        """
        if self.p_cal_f is not None:
            inds_ton = self.p_cal_f['inds_ton'][:]
        else:
            inds_ton, inds_toff_bef, inds_toff_aft = self.inds_ton, self.inds_toff_bef, self.inds_toff_aft
        # determine which cal power sample to calibrate
        tcal_c= inds_ton[:, inds_ton.shape[1]//2]
        inds_in_tcal_on= np.argmin(abs(inds_on[:,None]- tcal_c[None,:]), axis=1)
        inds_in_tcal_off= np.argmin(abs(inds_off[:,None]- tcal_c[None,:]), axis=1)
        uni= np.unique(np.hstack([inds_in_tcal_on,inds_in_tcal_off]))
        # check cal
        if self.p_cal_f is None:
            is_use = self.check_cal(inds_ton[uni])
            inds_tmp = np.where(is_use)[0]
            if len(inds_tmp) > 0:
                n_c_aft = len(is_use) - inds_tmp[-1] - 1
            else:
                n_c_aft = len(is_use)
            inds_ton_aft = inds_ton[uni[-1]+1 : uni[-1]+1+n_c_aft]
            # check more cal
#             print(inds_ton_aft)
            if len(inds_ton_aft) > 0:
                is_use_aft = self.check_cal(inds_ton_aft)
                inds_ton_use = np.vstack([inds_ton[uni][is_use], inds_ton[uni[-1]+1 : uni[-1]+1+n_c_aft][is_use_aft]])
                inds_toff_bef_use = np.vstack([inds_toff_bef[uni][is_use], inds_toff_bef[uni[-1]+1 : uni[-1]+1+n_c_aft][is_use_aft]])
                inds_toff_aft_use = np.vstack([inds_toff_aft[uni][is_use], inds_toff_aft[uni[-1]+1 : uni[-1]+1+n_c_aft][is_use_aft]])
            else:
                inds_ton_use = inds_ton[uni][is_use]
                inds_toff_bef_use = inds_toff_bef[uni][is_use]
                inds_toff_aft_use = inds_toff_aft[uni][is_use]
#             print(inds_ton_use)
#             print(inds_toff_bef_use)
#             print(inds_toff_aft_use)
        else:
            inds_ton_use = inds_ton[uni]
        # load cal power
        if self.p_cal_f is not None:
            p_cal_s = self.p_cal_f['p_cal_s'][uni]
        else:
            p_cal_s = self._get_cal_power_s(inds_ton_use, inds_toff_bef_use, inds_toff_aft_use)
        # load on and off power
        p_on = self.get_field(inds_on, 'DATA',)
        p_off = self.get_field(inds_off, 'DATA', close_file=True)
        #print(p_on.dtype, p_off.dtype, p_cal_s.dtype)
        try:
            if self.plot:
                figname = self.out_name_base + "-sep.pdf"
                plot_sep(inds_on, inds_off, p_on, p_off, figname=figname)
        except AttributeError:
            pass

        # redetermine which cal to use
        tcal_c = inds_ton_use[:, inds_ton_use.shape[1]//2]
        inds_in_tcal_on = np.argmin(abs(inds_on[:,None]- tcal_c[None,:]), axis=1)
        inds_in_tcal_off = np.argmin(abs(inds_off[:,None]- tcal_c[None,:]), axis=1)
        count_on = p_on.astype('float64') / p_cal_s[inds_in_tcal_on] - 1 # have subtracted cal
        count_off = p_off.astype('float64') / p_cal_s[inds_in_tcal_off]
        return count_on, count_off, p_cal_s, inds_ton_use


class PositionSwitch(CalOnOff):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # freq: self.freq_use
        self.mjd_a = self.get_mjds()
        # noise-diode on
        self.mjd_on = self.mjd_a[self.inds_on]
        # noise-diode off
        self.mjd_off = self.mjd_a[self.inds_off]

        # power of noise-diode on
        self.p_on = self.get_field(self.inds_on, field='DATA')  # power, shape (mjd,freq,polarization)
        # power of noise-diode off
        self.p_off = self.get_field(self.inds_off, field='DATA')  # power, shape (mjd,freq,polarization)

        self.Tcal_s = self.get_Tcal_s()


    @staticmethod
    def in_range(arr, starts, length):
        """
        arr: shape (m,)
        starts: shape (n,)
        length: scalar
        """
        is_in =  (arr[:,None] > starts[None,:]) & (arr[:,None] < (starts+length)[None,:])
        return np.logical_or.reduce(is_in, axis=1)

    def sep(self, t_src, t_ref, n_repeat, t_change=30):
        """
        second
        """
        #######################
        # on: noise-diode on
        # off: noise-diode off
        # src: on-source
        # ref: off-source
        #######################

        # second to day
        t_src = t_src/60/60/24
        t_ref = t_ref/60/60/24
        n_repeat = n_repeat
        t_change = t_change/60/60/24

        mjd_src_start = self.mjd_a[0] + np.arange(n_repeat)*(t_src + t_change + t_ref + t_change)
        mjd_ref_start = mjd_src_start + (t_src + t_change)

        is_ = self.in_range(self.mjd_on, mjd_src_start, t_src)
        self.p_on_src = self.p_on[is_]
        self.inds_on_src = self.inds_on[is_]

        is_ = self.in_range(self.mjd_off, mjd_src_start, t_src)
        self.p_off_src = self.p_off[is_]
        self.inds_off_src = self.inds_off[is_]

        is_ = self.in_range(self.mjd_on, mjd_ref_start, t_ref)
        self.p_on_ref = self.p_on[is_]
        self.inds_on_ref = self.inds_on[is_]

        is_ = self.in_range(self.mjd_off, mjd_ref_start, t_ref)
        self.p_off_ref = self.p_off[is_]
        self.inds_off_ref = self.inds_off[is_]

        self.mjd_src_start = mjd_src_start
        self.mjd_ref_start = mjd_ref_start

    def plot_sep(self,):
        try:
            if getattr(self, 'plot', False):
                figname = self.out_name_base + "-sep.pdf"
                plot_sep(self.inds_on, self.inds_off, self.p_on, self.p_off, figname=figname)
        except:
            pass

    def gen_Ta(self, only_off=False):
        """
        using

        self.p_on_src
        self.p_off_src
        self.p_on_ref
        self.p_off_ref

        to gen

        self.Ta as shape (1,chan,polar)

        """

        p_cal_src = np.mean(self.p_on_src, axis=0, dtype='float64', keepdims=True) - \
                    np.mean(self.p_off_src, axis=0, dtype='float64', keepdims=True)
        p_cal_ref = np.mean(self.p_on_ref, axis=0, dtype='float64', keepdims=True) - \
                    np.mean(self.p_off_ref, axis=0, dtype='float64', keepdims=True)
        # smooth power of cal
        p_cal_src_s = self._get_smoothed(p_cal_src)
        p_cal_ref_s = self._get_smoothed(p_cal_ref)

        Ta_on_src = (np.mean(self.p_on_src, axis=0, dtype='float64', keepdims=True)/p_cal_src_s - 1)*self.Tcal_s
        Ta_off_src = (np.mean(self.p_off_src, axis=0, dtype='float64', keepdims=True)/p_cal_src_s)*self.Tcal_s
        Ta_on_ref = (np.mean(self.p_on_ref, axis=0, dtype='float64', keepdims=True)/p_cal_ref_s - 1)*self.Tcal_s
        Ta_off_ref = (np.mean(self.p_off_ref, axis=0, dtype='float64', keepdims=True)/p_cal_ref_s)*self.Tcal_s

        # use interaged time as weight
        if not only_off:
            Ta = np.average(np.r_[Ta_on_src, Ta_off_src], axis=0, weights=[len(self.inds_on_src),len(self.inds_off_src)]) - \
                 np.average(np.r_[Ta_on_ref, Ta_off_ref], axis=0, weights=[len(self.inds_on_ref),len(self.inds_off_ref)])
        else:
            Ta = np.average(Ta_off_src, axis=0) - \
                 np.average(Ta_off_ref, axis=0)


        Ta = Ta[None,]

        self.Ta = Ta

        # also change mjd
        inds_tmp = self.inds_off_src[:1]
        self.mjd = self.mjd_a[inds_tmp]

    def gen_radec(self,):
        """
        try to calculation for each record: self.ra_a, self.dec_a
        use the position of first src record as self.ra, self.dec
        also set self.mjd and self.mjd_a

        """

        from .radec import get_radec
        nB = self.nB
        guess_str = self.fname_part

        radec = get_radec(self.mjd_a, guess_str=guess_str, nBs=[nB,])
        self.ra_a, self.dec_a = radec[f'ra{nB}'], radec[f'dec{nB}']

        inds_tmp = self.inds_off_src[:1]
        self.ra, self.dec = self.ra_a[inds_tmp], self.dec_a[inds_tmp]
        self.mjd = self.mjd_a[inds_tmp]

    def plot_radec(self, figsize=[10,10], outname=None):
        """
        plot the ra dec to check the separation of src and ref
        """

        from matplotlib import pyplot as plt
        plt.figure(figsize=figsize)
        plt.scatter(self.ra_a[self.inds_on_src], self.dec_a[self.inds_on_src], s=1, color='r')
        plt.scatter(self.ra_a[self.inds_off_src], self.dec_a[self.inds_off_src], s=1, color='r', label='src')
        plt.scatter(self.ra_a[self.inds_on_ref], self.dec_a[self.inds_on_ref], s=1, color='b')
        plt.scatter(self.ra_a[self.inds_off_ref], self.dec_a[self.inds_off_ref], s=1, color='b', label='ref')
        plt.xlabel('ra')
        plt.ylabel('dec')
        plt.grid()
        plt.minorticks_on()
        plt.legend()
        if outname is not None:
            print(f'Saving ra dec plot to {outname}')
            plt.savefig(outname)
