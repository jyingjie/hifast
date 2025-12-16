

__all__ = ['get_ratio', 'ATemp_correction', 'download_flux_gain_Liu2024', 'get_flux_gain_Liu2024', 'Calibrator',
           'FluxCali']



import os
import re
from glob import glob
import warnings

import h5py
import numpy as np
import scipy.interpolate as interp
import json

from . import conf
from astropy.time import Time
from astropy import units

from .gain import Get_gain, gain_diff_from_ZA, Get_ZA
from ..utils.downloader import get_file


def get_ratio(nB, freq=None):
    """get gain ratio to beam 1 in arXiv:2002.01786"""
    
    with open(os.path.dirname(__file__) + '/data/beam_ratios.json', 'r') as f:
        ratios_para = json.load(f)

    # Assumes JSON structure: "M01": {"1050": val, "1050_err": err, ...} (ordered)
    # Original: freq_key = np.array(ratios_para.index[::2], dtype=float)
    # ratios = ratios_para[f'M{nB:02d}'][::2].values
    
    beam_data = ratios_para[f'M{nB:02d}']
    # Keys/Values are ordered in Python 3.7+
    keys = list(beam_data.keys())
    values = list(beam_data.values())
    
    freq_key = np.array(keys[::2], dtype=float)
    ratios = np.array(values[::2])
    
    if freq is not None:
        ratios= interp.interp1d(freq_key,ratios, kind='quadratic', fill_value= "extrapolate")(freq)
        freq_key=freq
    return ratios, freq_key


## flux gain from Liu et al. 2024
def ATemp_correction(ATemp_pre, ATemp_target, freq):
    """
    ATemp: Ambient temperature (in degrees Celsius) of the pre-measured flux gains
    ATemp_target: Ambient temperature (in degrees Celsius) of the target observation.
    freq: frequency (in MHz) of the pre-measured flux gains
    """
    s_factor = 4.7e-5*freq - 1.4e-3 # equation from Liu et al. 2024
    K_Jy_diff = s_factor * (ATemp_target - ATemp_pre)
    return K_Jy_diff

def download_flux_gain_Liu2024():
    url = 'https://download.scidb.cn/download?fileId=5a612e328587f62c27a3435aac1af0fa'
    fpath = os.path.expanduser('~/.cache/hifast/data/flux_gain_Liu2024.hdf5')
    return get_file(url, fpath, max_retries=60, retry_delay=10)


def get_flux_gain_Liu2024(nB, ZA_target, ATemp_target=None, freq=None):
    """
    nB: scalar
        beam number
    ZA_target: array_like
        Zenith angle (in degrees) of the target observation
    ATemp_target: scalar
        Ambient temperature (in degrees Celsius) of the target observation
    freq: array_like
        frequency (in MHz) of the pre-measured flux gains

    """
    fpath = download_flux_gain_Liu2024()
    if fpath is None:
        raise(ValueError('can not obtain flux gain from Liu et al. 2024'))
    # load
    print(f'loading {fpath}')
    with h5py.File(fpath, 'r') as fs:
        freq_pre = fs['freq'][()]
        ATemp_pre = fs['Ambient_Temperature'][()]
        ZA_pre = fs['ZA'][()]
        K_Jy = fs[f'K_Jy{nB}'][()] # shape (1,n,1)

    if freq is not None:
        K_Jy = interp.interp1d(freq_pre, K_Jy, kind='linear', fill_value= "extrapolate", axis=1)(freq)
        freq_sample = freq
    else:
        freq_sample = freq_pre

    if ATemp_target is not None:
        print(f'Ambient temperature correction from {ATemp_pre} degree C to the target {ATemp_target} degree C')
        K_Jy_diff = ATemp_correction(ATemp_pre, ATemp_target, freq_sample)
        K_Jy += K_Jy_diff[None, :, None]
    else:
        warnings.warn(f'no ambient temperature correction')
    # Zenith angle correction
    if ZA_target is not None:
        K_Jy_diff = gain_diff_from_ZA(ZA_pre, ZA_target, nB, freq_sample)
        K_Jy = K_Jy + K_Jy_diff[..., None]
        print(f'Zenith angle correction from {ZA_pre} degree to the target {ZA_target} degree')
    else:
        warnings.warn(f'no Zenith angle correction')

    return K_Jy, freq_sample


class Calibrator:

    def __init__(self, mjd, nB, cbr_store, cbr_name='*', use_counts=True, only_use_19beams=False):
        """
        mjd: used to find nearest cbr
        nB: beam number
        cbr_store: str
           where quasar calibrator stored in, a file path or a directory including several files
           If None, use the gain depended on Zenith angle (arxiv:2002.01786) and need input ra, dec and mjd.
        cbr_name: str, e.g. 3C48
           If not None and cbr_store is a directoy, searching the specified calibrator.
        use_counts: bool
           If True, return gain as counts/Jy, If False, K/Jy
        """
        self.cbr_store = cbr_store
        self.cbr_name = cbr_name
        self.mjd = mjd
        self.nB = nB
        self.use_counts = use_counts
        self.only_use_19beams = only_use_19beams
        self.parser_cbr()

    @staticmethod
    def get_date_nearest(fpath_list, mjd):
        dates_have = [os.path.basename(fpath).rsplit('-FluxGain', 1)[0].rsplit('-', 1)[-1] for fpath in fpath_list]
        dates_have = [f'{s[:4]}-{s[4:6]}-{s[6:8]} {s[8:10]}:{s[10:12]}:00.000' for s in dates_have]
        # date of file name is UTC
        mjds_h = Time(dates_have, format='iso', scale='utc').mjd
## use mjd
#         mjds_h = np.empty(len(fpath_list), dtype='float64')
#         for i, fpath in enumerate(fpath_list):
#             with h5py.File(fpath) as fs:
#                 mjds_h[i] = fs['mjd'][0]
        ind = np.argmin(abs(mjds_h - mjd))
        return fpath_list[ind]

    def _calc_ratio(self, fpath):
        """calc gain ratio to beam 1"""
        with h5py.File(fpath, 'r') as fs:
            freq = fs['freq'][()]
            K_Jy = fs[f'K_Jy{self.nB}'][0] if not self.use_counts else fs[f'count_Jy{self.nB}'][0]
            K_Jy_nB1 = fs[f'K_Jy1'][0] if not self.use_counts else fs[f'count_Jy1'][0]
#             print(K_Jy,K_Jy_nB1)
            ratio = K_Jy/K_Jy_nB1
            ZA_nB = self.get_ZA_cbr(fs, self.nB)
        return freq, ratio, ZA_nB

    def get_ratio(self):
        # cbr_store is directory
        fpaths_nB_all = glob(os.path.join(self.cbr_store, f"**/{self.cbr_name}-20*-FluxGain-All.hdf5"), recursive=True)
        if len(fpaths_nB_all) == 0:
            raise(ValueError('can not find 19-beams calibration file to calculate relative gain'))
        else:
            fpath = self.get_date_nearest(fpaths_nB_all, self.mjd)
            print(f"using {fpath} to calculate relative gain to M01")
            freq, ratio, ZA_nB = self._calc_ratio(fpath)
        return freq, ratio, ZA_nB

    def parser_cbr(self,):
        """assign or serach fpath from self.cbr_store"""
        if not os.path.exists(self.cbr_store):
            raise(ValueError(f'can not find the calibrator: {self.cbr_store}'))
        if os.path.isfile(self.cbr_store):
            self.cbr_fpath = self.cbr_store
        else:
            if self.only_use_19beams:
                fpaths_cand = glob(os.path.join(self.cbr_store, f"**/{self.cbr_name}-20*-FluxGain-All.hdf5"), recursive=True)
            else:
                fpaths_cand = glob(os.path.join(self.cbr_store, f"**/{self.cbr_name}-20*-FluxGain*.hdf5"), recursive=True)
            if len(fpaths_cand) == 0:
                raise(ValueError(f"can not find calibrator file"))
            self.cbr_fpath = self.get_date_nearest(fpaths_cand, self.mjd)

    def get_ZA_cbr(self, fs, nB):
        if f'ZA_cbr{nB}' in fs.keys():
            ZA = fs[f'ZD_cbr{nB}'][()]
        else:
            ZAs = fs[f'ZD{nB}'][()]
            ZA = np.mean(ZAs)
        # need deg
        ZA = np.rad2deg(ZA)
        return ZA

    def __call__(self, freq_new=None, ZAs_obs=None, tcal_spec=None):
        """
        freq_new: if not None, will interplate
        ZAs_obs: ZA of obs spec. If not None, will fix the gain diff from ZA diff
        """
        print(f'loading gain from {self.cbr_fpath}')
        with h5py.File(self.cbr_fpath) as fs:
            freq = fs['freq'][()]
            key_ = f'K_Jy{self.nB}' if not self.use_counts else f'count_Jy{self.nB}'
            if freq_new is None:
                freq_new = freq
            if key_ in fs.keys():
                gain = fs[key_][()] # K/Jy, shape:(1, chan, 2)
                gain = interp.CubicSpline(freq, gain, axis=1, extrapolate=True)(freq_new)
            else:
                # use M01 or relative gain  ratio
                gain_nB1 = fs[f'K_Jy1'][()] if not self.use_counts else fs[f'count_Jy1'][()]
                gain_nB1 = interp.CubicSpline(freq, gain_nB1, axis=1, extrapolate=True)(freq_new)

                freq_ratio, ratio, ZA_cbr_nB_ratio = self.get_ratio()
                ratio = interp.CubicSpline(freq_ratio, ratio, extrapolate=True)(freq_new)
                gain = gain_nB1 * ratio[None,:]

            if self.use_counts:
                # multiply count with spectra's tcal
                gain *= tcal_spec # to K/Jy

            if ZAs_obs is not None:
                # fix ZA diff
                ZA_cbr = self.get_ZA_cbr(fs, self.nB) if key_ in fs.keys() else ZA_cbr_nB_ratio
                diff = gain_diff_from_ZA(ZA_cbr, ZAs_obs, self.nB, freq_new)
                gain = gain + diff[..., None]
                print(f'fixing ZA diff from ZA_cbr:{ZA_cbr} to ZA_obs:{ZAs_obs}')
        return freq_new, gain



class FluxCali:

    def __init__(self, nB, freq, cbr_store=None, cbr_name='*', tcal_spec=None,  only_use_19beams=False,
                           pre_measured='Jiang2020', ATemp=None,
                           mjd=None, ra=None, dec=None, ZAs=None, fix_diff_ZA=False):
        """

        flux calibration using Calibator or Table Gain
        -----------------------
        T: array_like
           Temperature of the spectra. Shape is (m,n) or (m,n,2) i.e. (Mjd, channel) or (Mjd, channel, Polarization)
        nB: int
           Beam numbe
        freq: array_like, shape (n,)
        cbr_store: str
           where quasar calibrator stored in, a file path or a directory including several files
           If None, use the gain depended on Zenith angle (arxiv:2002.01786) and need input ra, dec and mjd.
        cbr_name: str, e.g. 3C48
           If not None and cbr_store is a directoy, searching the specified calibrator.
        tcal_spec: shape: (1,n,2); T.shape need be (m,n,2); if spectra and Calibator temperature used different Tcal, input this to fix it.
        ra, dec, mjd: None or array_like, shape (m,)

        pre_measured: str
           'Jiang2020' or 'Liu2024'
        Atemp: array_like
           Ambient temperature (in degrees Celsius) of the target observation

        """

        if ra is not None and dec is not None:
            if np.isscalar(mjd):
                mjd = np.array([mjd])
            if np.isscalar(ra):
                ra = np.array([ra])
            if np.isscalar(dec):
                dec = np.array([dec])
            assert ra.shape[0] == dec.shape[0] == mjd.shape[0]

        self.nB = nB
        self.freq = freq
        self.cbr_store = cbr_store
        self.cbr_name = cbr_name
        self.tcal_spec = tcal_spec
        self.only_use_19beams = only_use_19beams
        self.mjd = mjd

        self.ra = ra
        self.dec = dec
        self.ZAs = ZAs
        self.fix_diff_ZA = fix_diff_ZA

        self.ATemp = ATemp

        if cbr_store is None or cbr_store=='none':
            if pre_measured == 'Jiang2020':
                print('Using pre-measured flux gain from Jiang2020')
                self.K_Jy = self.get_gain_tabled()
            elif pre_measured == 'Liu2024':
                print('Using pre-measured flux gain from Liu2024')
                self.K_Jy = self.get_gain_Liu2024()
            else:
                raise ValueError(f'pre_measured {pre_measured} not supported')
        else:
            print(f'Using Calibrator from {cbr_store}')
            self.K_Jy = self.get_gain_cbr()

    def get_gain_cbr(self,):

        use_counts = True if self.tcal_spec is not None else False

        CBR = Calibrator(np.mean(self.mjd), self.nB, self.cbr_store,
                         use_counts=use_counts,
                         only_use_19beams=self.only_use_19beams)
        if self.fix_diff_ZA:
            if self.ZAs is not None:
                ZAs = self.ZAs
            else:
                ZAs = Get_ZA(self.ra, self.dec, self.mjd)
        else:
            ZAs = None
        _, gain = CBR(self.freq, ZAs, self.tcal_spec)
        self.cbr_fpath = CBR.cbr_fpath
        K_Jy = gain
        return K_Jy


    def get_gain_tabled(self,):
        """from Jiang2020"""
        K_Jy = Get_gain(self.ra, self.dec, self.mjd, self.nB, self.freq)[0]  # K/Jy

        return K_Jy[:, :, None]

    def get_gain_Liu2024(self,):
        """from Liu2024"""
        if self.ZAs is not None:
            ZA_target = self.ZAs
        else:
            ZA_target = Get_ZA(self.ra, self.dec, self.mjd)
        K_Jy, _ = get_flux_gain_Liu2024(self.nB, ZA_target=ZA_target, ATemp_target=self.ATemp, freq=self.freq)
        return K_Jy

    def __call__(self, T):
        """
        T: array; shape: (Mjd,Chan,Polar) or (Mjd,Chan)
        """
        K_Jy = self.K_Jy
        # cali
        if T.ndim == 2:
            return T/np.mean(K_Jy, axis=2)
        elif T.ndim == 3:
            if T.shape[-1] == 2:
                return T/K_Jy
            elif T.shape[-1] == 1:
                return T/np.mean(K_Jy, axis=2, keepdims=True)
            else:
                raise(ValueError('shape of T'))
        else:
            raise(ValueError('shape of T'))
