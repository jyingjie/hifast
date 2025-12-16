

__all__ = ['robust_smooth', 'smooth_by_baseline', 'CbrInfo', 'CbrProfileFile', 'fpath', 'CbrProfiles', 'calnames',
           'gen_args', 'parse_obs_log', 'Calibrator', 'MBC', 'OnOff', 'Mapping', 'MappingType2', 'main']


# Author: Nekomata  zmliu@nao.cas.cn
# YJ Jing

from traceback import print_exc

from ..utils.io import ArgumentParser, bool_fun, formatter_class, replace_nB, sub_patten
from ..utils import obs_log
from ..utils.io import save_dict_hdf5
from ..core.radec import get_radec
from ..core.cal import CalOnOff
from ..core.cal2 import CalOnOffM


from astropy.coordinates import SkyCoord
from astropy import units as u
import os
import sys
import re
import warnings

import numpy as np
import scipy.interpolate as interp

from astropy.time import Time


from ..utils.misc import smooth1d
from ..utils.misc import apply_along_axis
from scipy import interpolate as interp
@apply_along_axis
def robust_smooth(T, freq, sigma, sigma_bad=0, dlimit=0.025, return_outlier=False):
    """
    Interpolate NaN values, mitigate outliers, and smooth a 1D array

    freq: frequency in MHz
    sigma: sigma for smoothing
    sigma_bad: sigma for bad data, default 0
    dlimit: limit of the derivative of the smoothed data, default 0.025
    """

    # from Mhz to number of channels
    nch_1MHz = (len(freq)-1) / (freq[-1] - freq[0])
    sigma = int(sigma * nch_1MHz)
    sigma_bad = int(sigma_bad * nch_1MHz)

    # Interpolate over NaN values first
    nans = np.isnan(T)
    if np.any(nans):
        good_data_indices = ~nans
        T_interpolator = interp.interp1d(freq[good_data_indices], T[good_data_indices], kind='linear', fill_value='extrapolate')
        T = T_interpolator(freq)

    if sigma_bad > 0:
        T_smoothed = smooth1d(T, 'gaussian_fft', sigma_bad, axis=0)
        T_diff = np.abs(np.diff(T_smoothed))
        exceed_diff = T_diff / ((T_smoothed[:-1] + T_smoothed[1:]) / 2) * nch_1MHz >= dlimit
        k = 2*sigma_bad

        exceed_indices = np.where(exceed_diff)[0]
        mask_indices = (exceed_indices[:, None] + np.arange(-k, k + 1)).flatten()
        mask_indices = np.clip(mask_indices, 0, len(T) - 1)
        is_outlier = np.full(len(T), False)

        np.put(is_outlier, mask_indices, True)

        T_non_outlier = T[~is_outlier]
        interpolator = interp.interp1d(freq[~is_outlier], T_non_outlier, kind='linear', fill_value='extrapolate')
        T = interpolator(freq)
    else:
        is_outlier = None

    Ts = smooth1d(T, 'gaussian_fft', sigma, axis=0)
    if return_outlier:
        return Ts, is_outlier
    else:
        return Ts


def smooth_by_baseline(T, freq, sigma, axis=1):
    """
    Smooth the input array `T` along the specified `axis` using a baseline fitting method.

    Parameters:
    - T: numpy.ndarray
        The input array to be smoothed.
    - freq: numpy.ndarray
        The frequency values corresponding to the data points in `T`.
    - sigma: float
        The standard deviation of the Gaussian kernel used for smoothing. The value is specified in MHz.
    - axis: int, optional
        The axis along which the smoothing is performed. Default is 1.

    Returns:
    - ys: numpy.ndarray
        The smoothed array by fitting baseline.
    """
    from ..core.baseline import get_baseline

    # Convert frequency from MHz to number of channels
    nch_1MHz = (len(freq) - 1) / (freq[-1] - freq[0])
    sigma = int(sigma * nch_1MHz)

    bl_para = {'lam': sigma, 'ratio': 0.01}
    method = 'Gauss-sym1'

    ys = get_baseline(freq, T, axis=axis, method=method, bl_para=bl_para, check=True, interp_nan=True)

    return ys


import json
class CbrInfo():
    def __init__(self, ra, dec, coeff=None, coeff_err=None):
        self.ra = ra
        self.dec = dec
        self.coeff = coeff
        self.coeff_err = coeff_err

    def get_profile(self, freq: 'GHz') -> 'Jy':
        freq_log10 = np.log10(freq)
        flux = 10**np.polyval(self.coeff[::-1], freq_log10)
        return flux

class CbrProfileFile():
    def __init__(self, fname):
        self.fname = fname
        with open(fname) as f:
            self.data = json.load(f)
        self.names = list(self.data.keys())

    def __call__(self, name):
        if name not in self.data:
            return None
        return CbrInfo(self.data[name]['ra'], self.data[name]['dec'], self.data[name]['coeff'], self.data[name]['coeff_err'])


fpath = os.path.dirname(__file__) + '/data/PerleyButler2017.0.json'
CbrProfiles = CbrProfileFile(fpath)
calnames = CbrProfiles.names


def gen_args(argv=None):

    sep_line = '##'+'#'*70+'##'
    parser = ArgumentParser(prog=f"python -m hifast.{os.path.basename(sys.argv[0])[:-3]}", formatter_class=formatter_class,
                            allow_abbrev=False,
                            description='Processing the calibrator data', )
    parser.add_argument('fname',
                        help="Path to the raw data file. This should be used in conjunction with the `--nBs` parameter. Example: '*M01*-0001.fits'.")
    parser.add_argument('--nBs', type=int, nargs='*',
                        help='Beam numbers to process. Input integers separated by space, such as `1 3 11`. \
                            If no numbers are provided, all beams will be processed when the `--obsmode` is `MultiBeamCalibration` or `MultiBeamOTF`. \
                                Otherwise, only beam 1 will be processed.', default=None)
    parser.add_argument('--outdir', type=str, required=True,
                    help='Path to the directory where output files will be saved.',default=None)
    parser.add_argument('-f', action='store_true', dest='force',
                        help='if set, overwriting file if output file exists')

    parser.add_argument('--nproc', '--np', type=int, default=3,
                   help='number of process used')
    parser.add_argument('--cache', type=bool_fun, choices=[True, False], default='False',
                       help='If set to True and not using `-f`, will cache/load the results of all beams except for the initial one.')
    parser.add_argument('--obs_log',
                       help='If set, will use the obs_log to get `-d`, `-m`, `-n`, `--noise_mode`.\
                        It requires a file path ending in `.log` or `.txt`, or a directory path where the log file can be searched for.')

    group = parser.add_argument_group(f'*Calibrator setting\n{sep_line}')
    group.add_argument('--obsmode',type=str, choices=['MultiBeamCalibration', 'OnOff', 'MultiBeamOTF','DriftWithAngle','DecDriftWithAngle','Drift'],
                        help='Observation mode. This can be `MultiBeamCalibration`, `OnOff`, or some special scanning modes.',)
    group.add_argument('--fit_sep', type=bool_fun, choices=[True, False], default='False',
                        help='fit Cal on and Cal off separately. Used in scanning mode.')
    group.add_argument('--cbrname', '--calname', type=str, required=True,
                    help='Name of the calibrator, for example, "3C48".',)
    group.add_argument('--crd', type=float, nargs=2,
                    help='Coordinates of the calibrator in degrees. Input two float numbers for RA and DEC. \
                        If not provided, the program will use Simbad.query_object to query by `--calname`.', default=None)


    group.add_argument('--fluxProfilePara', type=str, nargs='*',
                    help=(f'Input the flux density profile parameters for calibrator. If the calibrator is not included in the list of {calnames}, '
                           'you may enter several floating-point numbers: ``a0, a1, a2, a3, ...``. '
                           'The function will then be represented by the equation ``a0 + a1*log10(f) + a2*log10(f)**2 + a3*log10(f)**3 + ...``. '
                           'Alternatively, you can provide a path (ends with `.hdf5`) of a hdf5 file that contains fields "freq" (frequency in MHz) and "flux" (flux density in Jy).'
                          ),
                           default=None)
    group.add_argument('--s_method', choices=['A', 'B'],  default='B',
                    help='A and B both utilize Gaussian smoothing but with different strategies to suppress RFI. Default is B.')
    group.add_argument('--s_sigma_cbr', type=float,
                    help='Smooth sigma (in MHz) of Gaussian smooth along frequency for calibrator spectra.', default=10)

    group = parser.add_argument_group(f'*Parameters with hifast.radec \n{sep_line}')
    group.add_argument('--tol', type=float, default=1,
                       help='max allowed extrapolate time; unit: second')
    group.add_argument('--ky_files', nargs='*',
                        help='KY files, if not given, try to search in .')


    group = parser.add_argument_group(f'*Parameters same with hifast.sep\n{sep_line}')
    group.add_argument('-d','--d', '--n_delay', type=int,
                    help='Time of delay divided by sampling time')
    group.add_argument('-m','--m', '--n_on', type=int,
                    help='Time of Tcal_on divided by sampling time')
    group.add_argument('-n','--n', '--n_off', type=int,
                        help='Time of Tcal_off divided by sampling time')
    group.add_argument('--Cal_fillratio', type=float, default=1,
                    help='Cal on fill ratio. when the time of Calon is less than one sampling time, {Calon_time}/{sampling_time}. Default is 1.')
    group.add_argument('--frange', type=float, nargs=2,
                    help='Frequency range. Input two float numbers.',default= [1000,1500] )
    group.add_argument('--noise_mode', default='high', choices=['high','low'],
                        help='Noise diode mode. This can be "high" or "low".')
    group.add_argument('--noise_date', default='auto',type=str,
                        help='Date of noise diode observation. Default is "auto".')
    group.add_argument('--s_sigma', '--smt_sigma', type=float,
                    help='Smooth sigma (in MHz) of Gaussian smooth along frequency for power and temperature of noise diode.', default=1)

    group = parser.add_argument_group(f'*Parameters in MultiBeamCalibration or OnOff\n{sep_line}')
    group.add_argument('--t_src', type=int,
                    help='Tracking time of On-Source in seconds. This is used in conjunction with Off-Source and change time to determine the On/Off source spectra.\
                        In the `OnOff` mode, the Off-Source time used is identical to the On-Source time, with a change time of 30 seconds. \
                    In the `MultiBeamCalibration` mode, the change time is 40 seconds.)', default=60)
    group.add_argument('--n_cir', '--n_repeat', type=int, default=1,
                    help='The number of on-source off-source cycles in `OnOff`. It is 1 in `MultiBeamCalibration`.')
    group.add_argument('--ref_tlim_nsrc', type=float, default=6,
                    help=('Assumed the starting time to track a beam is `Ts`, the spectra with time in \n'
                          '[Ts - t_change - ref_tlim_nsrc*t_src), Ts - t_change]\n'
                          'or\n[Ts + t_src+t_change, Ts + t_src + t_change + ref_tlim_nsrc*t_src)]\n'
                          'will be used as Off-Source spectra.'))
    group.add_argument('--src_drange', nargs=2, type=float, default=[0., 30.],
                    help='Limits the distance of On-Source spectra to the Calibrator coordinates in this range. The unit is arcsecond.')
    group.add_argument('--ref_drange', nargs=2, type=float, default=[648., 3600.],
                    help='Limits the distance of Off-Source spectra to the Calibrator coordinates in this range. The unit is arcsecond.')

    group = parser.add_argument_group(f'*Parameters in Mapping mode\n{sep_line}')
    group.add_argument('--dra_range_top', nargs=2, type=float, default=[0., 0.025],
                    help='The sampling points with absoulute RA separation in this range are used to fit the maximum temperature/power. The unit is deg.')
    group.add_argument('--dra_range_bottom', nargs=2, type=float, default=[0.18, 0.5],
                    help='The sampling points with absoulute RA separation in this range are used to fit the system temperature/power. The unit is deg.')
    group.add_argument('--ddec_max', type=float, default=0.025,
                    help='absoulute DEC separation limit for both. The unit is deg.')
    group.add_argument('--top_fit_method', type=str, choices=['gauss', 'skewgauss'], default='gauss',
                    help='Method used to fit the maximum temperature/power. This can be `gauss` or `skewgauss`.')


    group = parser.add_argument_group(f'*Others\n{sep_line}')
    group.add_argument('--saveT', type=bool_fun, choices=[True, False], default='False',
                    help='Choose whether to save T files or not.')

    group = parser.add_argument_group(f'*Plot\n{sep_line}')
    group.add_argument('--plot_near_radec', type=float, default=0.182,
                    help='Plot the sampling points near to the Calibrator in this range. The unit is deg.')
    group.add_argument('--plot', action='store_true',
                       help='if set, plot FluxGain with frequency into a png.')
    return parser.parse_args(argv)


def parse_obs_log(args):
    args = args

    if args.obs_log is not None:
        if not os.path.exists(args.obs_log):
            raise ValueError(f'`--obs_log` not exists: `{args.obs_log}`.')
        if os.path.isfile(args.obs_log):
            obs_log_path = args.obs_log
        elif os.path.isdir(args.obs_log):
            print(f'try to search for log file in the directory: `{args.obs_log}`.')
            obs_log_path = obs_log.search_log_fpath(args.fname, args.obs_log)
        else:
            raise ValueError(f'`--obs_log` should be a file path ending in `.log` or `.txt`, or a directory path where the log file can be searched for.\
            \nBut got `{args.obs_log}`.')
        print(f'Use obs_log file: `{obs_log_path}`.')
        paras = obs_log.get_sep_para_from_log(obs_log_path)[0]
        print(f'got paras: {paras}')
        args.d = paras['n_delay']
        args.m = paras['n_on']
        args.n = paras['n_off']
        args.noise_mode = paras['noise_mode']

    elif args.d is None or args.m is None or args.n is None or args.noise_mode is None:
        raise ValueError(f'`--obs_log` is not set, but `--d`, `--m`, `--n`, `--noise_mode` are not set either.')


class Calibrator():
    Cal_class = CalOnOff
    def __init__(self, args):
        self.args = args

        if args.nBs is None:
            if args.obsmode in ['MultiBeamCalibration','MultiBeamOTF',]:
                nBs= np.arange(1,20)
            else:
                nBs= [1,]
        else:
            nBs= args.nBs
        self.nBs_list = nBs
        self.nB_init = self.nBs_list[0]
        self._gen_sep_para()

        # init first beam to get mjd used in all beams and in ra-dec calculation
        self.init_sep(self.nB_init)
        self.gen_outname()
        self._check_fout()
        self.gen_mjd()


        self.freq_use = self.raw_sep.freq_use
        self.gen_crd()
        self.gen_flux_profile(self.freq_use)

    def _gen_sep_para(self):
        """
        """
        args = self.args
        para = {}
        para['n_delay'] = args.d
        para['n_on'] = args.m
        para['n_off'] = args.n
        para['start'] = 1
        para['stop'] = None
        para['frange'] = args.frange
        para['verbose'] = True
        para['smooth'] ='gaussian'
        para['s_para'] = {'s_sigma':args.s_sigma}
        para['dfactor'] = None
        para['med_filter_size'] = None
        para['noise_mode'] = args.noise_mode
        para['noise_date'] = args.noise_date
        self.sep_para = para

    def init_sep(self, nB):
        # load the raw data of the calibrator
        fpath = replace_nB(self.args.fname, nB)
        self.raw_sep = self.Cal_class(fpath, **self.sep_para)

    def gen_mjd(self):
        # get the mjd of the calibrator
        if not hasattr(self,'mjd'):
            self.mjd = self.raw_sep.get_mjds()

    def gen_outname(self):
        try:
            # DATE    = '2021-09-08T18:38:00Z', UTC
            obsdate = cbr.raw_sep.hd0s[0]['DATE']
        except:
            self.gen_mjd()
            obsdate= Time(self.mjd[0], format='mjd').to_value('isot')
        obsdate= re.sub('\D','',obsdate)[:12]
        outdir = sub_patten(self.args.outdir, date=obsdate+'UTC', project=self.args.cbrname)
        outdir = os.path.expanduser(outdir)
        if outdir!='' and (not os.path.exists(outdir)):
            print(f'outdir {outdir} not exists. Create it now')
            os.makedirs(outdir, exist_ok=True)
        self.args.outdir = outdir # update outdir
        self.outname= f"{outdir}/{self.args.cbrname}-{self.args.obsmode}-{obsdate}UTC"

    def _check_fout(self,):
        """
        check whether self.fout exists
        """
        args = self.args
        fpath_out = self.outname + '-FluxGain.hdf5'
        # or
        fpath_out2 = self.outname + '-FluxGain-All.hdf5'

        if os.path.exists(fpath_out) or os.path.exists(fpath_out2):
            if args.force:
                print(f"will overwrite the existing output file {fpath_out} or {fpath_out2}")
            else:
                print(f"File exists {fpath_out} or {fpath_out2}")
                print("exit... Use ' -f ' to overwrite it.")
                sys.exit(0)
        self._import()

    @staticmethod
    def _import():
        # run after _check_fout()
        from matplotlib import pyplot as plt
        import h5py
        # register them to global
        globals().update(locals())

    def gen_crd(self):
        """
        load the coordinates of the calibrator
        """
        args = self.args
        obj = args.cbrname
        crds = args.crd

        if crds is None:
            print('try to get coordinates of the calibrator using astroquery.simbad')
            try:
                from astroquery.simbad import Simbad
            except ImportError:
                print('Can not import astroquery.simbad, please input coordinates using ``--crd``.')
                sys.exit(1)
            try:
                sbo = Simbad.query_object(obj)
                crd = SkyCoord((sbo['RA'][0]+sbo['DEC'][0]),unit=(u.hourangle,u.deg))
            except:
                print('Calibrator can not find in Simbad.query_object, please input coordinates using ``--crd``.')
                sys.exit(1)
        else:
            print(f'Use input coordinates of the calibrator: {crds}')
            crd = SkyCoord(*crds,unit=(u.deg))
        print(f'Got coordinates of the calibrator: {crd}')
        self.crd = crd

    def plot_cbr_flux_profile(self, freq, flux):
        plt.plot(freq, flux)
        plt.grid()
        plt.xlabel('Frequency [MHz]')
        plt.ylabel('Flux [Jy]')
        plt.savefig(self.outname+'-cbr_FluxProfile.png')


    @staticmethod
    def get_cbr_flux_from_file(fpath, freq_sample):
        # freq_sample should in MHz
        with h5py.File(fpath, 'r') as f:
            freq = f['freq'][()] # in MHz
            flux = f['flux'][()] # in Jy
        # interpolate
        flux_sample = interp.interp1d(freq, flux)(freq_sample)
        return flux_sample

    def gen_flux_profile(self, freq_key):
        """
        freq_key should in MHz
        calname should be string
        """
        calname = self.args.cbrname

        fpparas = self.args.fluxProfilePara

        freq_key= freq_key/1000 # to GHz
        if fpparas is not None:
            # if input hdf5 file
            if fpparas[0].endswith('.hdf5'):
                print(f'Loading flux profile of the calibrator from {fpparas[0]}')
                flux = self.get_cbr_flux_from_file(fpparas[0], freq_key*1e3)
            # if polynomial cofficients (float nubmers) are input
            else:
                print(f'Get flux profile of the calibrator from polynomial cofficients (log): {fpparas}')
                fpparas = list(map(float, fpparas))
                c = CbrInfo(ra=None, dec=None, fpparas=fpparas)
                flux = c.get_profile(freq_key)
        else:
            if calname in CbrProfiles.names:
                c = CbrProfiles(calname)
                flux = c.get_profile(freq_key)
            else:
                raise ValueError('Can not find calibrator name in build-in flux profile. Please check calname or input flux profile parameters.')
                sys.exit(1)
        print('plot the flux profile of the calibrator')
        self.plot_cbr_flux_profile(freq_key*1e3, flux)
        self.flux_profile = flux

    def gen_T(self, nB):
        # process the calibrator data to get their temperature
        # self._init_sep(nB)
        pass

    def gen_ra_dec(self, nB):
        """calculate the ra dec of each spetrum for the calibrator"""
        if hasattr(self,'radec'):
            radec = self.radec
        else:
            radec = {}
        if nB == 1:
            print(f'Calculating RA-DEC for Beam {list(self.nBs_list)}')
            radec = get_radec(self.mjd, guess_str=self.args.fname,nBs=list(self.nBs_list),nproc=self.args.nproc,
                             tol=self.args.tol,
                             ky_files=self.args.ky_files)
            save_dict_hdf5(self.outname+'-radec.hdf5',radec)
            print('Saved RA-DEC to '+self.outname+'-radec.hdf5')
        else:
            try:
                _,_ = radec[f'ra{nB}'], radec[f'dec{nB}']
            except:
                print(f'Calculating RA-DEC for Beam {nB:02d}')
                radec = get_radec(self.mjd, guess_str=self.args.fname,tol=1,nBs=[nB,])
        self.ra, self.dec = radec[f'ra{nB}'], radec[f'dec{nB}']
        self.ZD = radec[f'ZD{nB}']
        self.radec = radec

        self._plot_radec_near_cbr(nB)

    def _plot_radec_near_cbr(self, nB, r_max=0.081):
        """plot the RA-DEC near the CBR"""
        from matplotlib import pyplot as plt

        if self.args.plot_near_radec is not None:
            r_max = self.args.plot_near_radec
        ra = self.ra
        dec = self.dec
        ra_cbr = self.crd.ra.deg
        dec_cbr = self.crd.dec.deg
        figname = self.outname + f'-{nB:02d}-radec.png'

        is_use = (ra < ra_cbr + r_max/np.cos(np.deg2rad(dec_cbr)) ) & (ra > ra_cbr-r_max/np.cos(np.deg2rad(dec_cbr)))
        is_use &= (dec < dec_cbr+r_max ) & (dec > dec_cbr-r_max)

        plt.figure(figsize=(5,5))
        plt.scatter(ra[is_use], dec[is_use], marker='.', s=1, color='k')
        plt.scatter(ra_cbr, dec_cbr, marker='.', s=7, color='r')
        plt.xlabel('RA [deg]')
        plt.ylabel('DEC [deg]')
        plt.grid()
        plt.minorticks_on()
        plt.suptitle(os.path.basename(figname))
        plt.tight_layout()
        plt.savefig(figname)

    def gen_flux_gain(self, nB):

        self.K_Jy = self.Ta_s/self.flux_profile[None,:, None]
        self.count_Jy = self.count_s/self.flux_profile[None,:, None]

    def process_beam(self, nB, check_cache=True):

        if check_cache:
            cache_name = self.outname + f'-M{nB:02d}-cache.h5'
            if os.path.exists(cache_name):
                import h5py
                print(f'Load FluxGain of beam {nB:02d} from cached file: {cache_name}')
                outdata = {k:v[()] for k,v in h5py.File(cache_name,'r').items()}
                for key in outdata.keys():
                    if isinstance(outdata[key], bytes):
                        outdata[key] = outdata[key].decode('utf-8')
                return outdata

        outdata = {}
        try:
            self.init_sep(nB)
            self.gen_ra_dec(nB)
            self.gen_T(nB)
            self.gen_flux_gain(nB)
        except Exception:
            print(f'Fail to process Beam {nB:02d}')
            print('Error:')
            print_exc()
            return None

        outdata[f'K_Jy{nB}'] = self.K_Jy
        outdata[f'count_Jy{nB}'] = self.count_Jy
        outdata[f'Ta{nB}'] = self.Ta
        outdata[f'Ta_s{nB}'] = self.Ta_s
        outdata[f'count{nB}'] = self.count
        outdata[f'count_s{nB}'] = self.count_s

        outdata[f'Tcal{nB}'] = self.Tcal_s
        outdata[f'ra{nB}']= self.ra
        outdata[f'dec{nB}']= self.dec
        outdata[f'ZD{nB}']= self.ZD
        outdata[f'ZD_cbr{nB}']= self.ZD_cbr
        # outdata[f'ONp{nB}']= self.ONp
        # outdata[f'OFFp{nB}']= self.OFFp
        # outdata[f'Tmax{nB}']= self.TMAX
        # {cbr.fit_k*100:.3f}%
        outdata[f'frac_fail{nB}'] = 1-self.frac_src_used
        print('Used Tcal of ' + self.tcal_file)
        print(f'Calibration of M{nB:02d} complete. {(1-self.frac_src_used)*100:.3f}% of input data ignored or failed to fit.')
        print('#'*50)

        # common
        outdata['Tcal_file']= self.tcal_file
        outdata['freq'] = self.freq_use
        outdata['flux_profile'] = self.flux_profile
        outdata['mjd'] = self.mjd
        #
        outdata['outname'] = self.outname
        if check_cache:
            print(f'Save FluxGain of beam {nB:02d} to cache file {cache_name}')
            save_dict_hdf5(cache_name, outdata)
        return outdata

    def save(self, outdata, nBs_have):
        if nBs_have == 19:
            print("19 beams have been successfully processed and saving to file with name containing 'All'.")
            outname = outdata['outname'] + '-FluxGain-All.hdf5'
        else:
            print("Not all 19 beams were processed successfully and saved to a file that doesn't have 'All' in its name.")
            outname= outdata['outname'] + '-FluxGain.hdf5'
        save_dict_hdf5(outname, outdata)
        print('Saved to '+ outname)

        if self.args.plot:
            from .plot import plot
            fig = plot(outdata, suptitle=outname)
            fig_path = outname.rsplit('.', 1)[0] + '.png'
            print(f'plot fulx gain to {fig_path}')
            fig.savefig(fig_path)

    def process_beam_mp(self, q, sema, *args, **kwargs):
        try:
            res = self.process_beam(*args, **kwargs)
        except SystemExit as e:
            res = None
        q.put(res)
        sema.release()

    def __call__(self):
        args = self.args
        nBs = self.nBs_list
        # process the first one, M01 in general, and the radec are generated
        print(f'Process Beam {self.nBs_list[0]:02d}')
        outdata = self.process_beam(nBs[0], False)
        if outdata is None:
            print('Fail to process the first beam in ``--nBs``. Abort.')
            sys.exit(1)
        nBs_have = 1
        if len(nBs) > 1:
            print(f'Process rest Beams with {args.nproc} CPUs')

            from functools import partial
            from multiprocessing import Process, Queue, Semaphore
            ps = []
            qs = []
            sema = Semaphore(args.nproc)

            use_cache = args.cache if not args.force else False

            for nB in nBs[1:]:
                sema.acquire()
                q = Queue()
                p = Process(target=self.process_beam_mp, args=(q, sema, nB, use_cache))
                ps += [p]
                qs += [q]

                p.start()
            result = [q.get() for q in qs]

            for j, d in enumerate(result):
                if d is not None:
                    outdata.update(d)
                    nBs_have += 1
                else:
                    print(f'Fail to process Beam {nBs[1:][j]:02d}')
        # save
        self.outdata = outdata
        self.nBs_have = nBs_have
        self.save(outdata, nBs_have)


class MBC(Calibrator):

    def _limit_r(self, inds, drange, nB, verbose=True):
        num_in = len(inds)
        crds = SkyCoord(self.ra[inds], self.dec[inds], unit=u.deg)
        r = crds.separation(self.crd).to_value(u.arcsec)
        is_ = (r>=drange[0])&(r<=drange[1])
        inds = inds[is_]
        num_use = len(inds)
        if num_use == 0:
            print(f'Beam {nB:02d}: Input {num_in} spectra, no spectra inside {drange} arcsec , please check input parameters.')
            sys.exit(1)
        else:
            if verbose:
                print(f'Beam {nB:02d}: Input {num_in} spectra, {num_use} spectra is in {drange} arcsec')
        return inds

    def limit_r_src_ref(self, inds_src, inds_ref, nB, verbose=True):
        if verbose:
            print(f'Beam {nB:02d}: Further restrict the on-source (src) spectra within {self.args.src_drange} arcsec')
        inds_src = self._limit_r(inds_src, self.args.src_drange, nB, verbose=verbose)
        if verbose:
            print(f'Beam {nB:02d}: Further restrict the off-source (ref) spectra within {self.args.ref_drange} arcsec')
        inds_ref = self._limit_r(inds_ref, self.args.ref_drange, nB, verbose=verbose)
        return inds_src, inds_ref

    def _sep_src_ref(self, inds_use, nB):
        """
        """
        t_src = self.args.t_src/60/60/24 # to day
        mjd = self.mjd[inds_use]
        ref_tlim_nsrc = self.args.ref_tlim_nsrc

        # obs order
        # 波束 1 重复跟踪一次。即多波束观测顺序：1-1-2-3-4-5-6-7-19-8-9-10-11-12-13-14-15-16-17-18。
        # 单次切换时间为 40 秒。
        t_change = 40/60/60/24 # fixed
        t_cir = t_src + t_change
        mjd_splits = self.mjd[0] + np.arange(20)*t_cir
        obslist = np.array([1,2,3,4,5,6,7,19,8,9,10,11,12,13,14,15,16,17,18])
        mjd_start = mjd_splits[np.where(obslist==nB)[0][0]+1]
        is_src = (mjd_start < mjd) & (mjd_start + t_src> mjd)
        if nB == 1:
            is_src |= (mjd<(self.mjd[0] + t_src))

        if nB == 1:
            ref_tlim_nsrc *= 2
        is_ref = (mjd > mjd_start + t_cir) & (mjd < mjd_start + t_cir + ref_tlim_nsrc*t_src)
        if nB != 1:
            is_ref |= (mjd > mjd_start - t_change - ref_tlim_nsrc*t_src) & (mjd < mjd_start - t_change)
        inds_src, inds_ref = inds_use[is_src], inds_use[is_ref]
        return inds_src, inds_ref

    def plot_with_t(self, ax, **kwargs):
        """
        ax: ax to plot, if is None, create a new one
        kwargs: label=[t, power]

        """
        spec = self.raw_sep
        # if [1395, 1405] in self.spec.freq_use, use it. else use from middle-5 to middle+5
        frange_cand = [1395, 1405]
        if frange_cand[0] >= spec.freq_use[0] and frange_cand[1] <= spec.freq_use[-1]:
            frange = frange_cand
        else:
            f_middle = np.median(spec.freq_use)
            frange = [f_middle-5, f_middle+5]
        is_ = (spec.freq_use>=frange[0]) & (spec.freq_use<=frange[1])

        if ax is None:
            fig, ax = plt.subplots()
        for i, label in enumerate(kwargs.keys()):
            t, p = kwargs[label]
            p_m = np.nanmedian(p[:,is_], axis=1)
            ax.plot(t, p_m[:,0], color='C'+str(i), linestyle='none', marker='o', ms=3,label=label+' XX',)
            ax.plot(t, p_m[:,1], color='C'+str(i), linestyle='-', ms=3,label=label+' YY',)
        ax.legend()
        ax.set_ylabel(f'Median Power in {frange} MHz')
        return ax


    def gen_T(self, nB):

        spec = self.raw_sep
        spec.gen_out_name_base(self.args.outdir)
        spec.plot = True

        # plot power with t
        fig, ax = plt.subplots(figsize=(10, 5))
        time_sec = (self.mjd - self.mjd[0])*24*60*60

        inds_off_src, inds_off_ref = self._sep_src_ref(spec.inds_off, nB)
        inds_off_src_bak = inds_off_src
        inds_off_src, inds_off_ref = self.limit_r_src_ref(inds_off_src, inds_off_ref, nB, verbose=True)
        self.frac_src_used = len(inds_off_src)/len(inds_off_src_bak)
        # ZD
        self.ZD_cbr = np.nanmedian(self.ZD[inds_off_src])

        p_off_src = spec.get_field(inds_off_src, 'DATA', close_file=False)
        p_off_ref = spec.get_field(inds_off_ref, 'DATA', close_file=False)

        # plot power with t
        self.plot_with_t(ax,
                         cal_off_src=[time_sec[inds_off_src], p_off_src],
                         cal_off_ref=[time_sec[inds_off_ref], p_off_ref],
                        )
        ax.set_xlabel('Time from beginning [Second]')
        figname = self.outname + f'-M{nB:02d}-power_time.png'
        ax.set_title(f'M{nB:02d}, ZD={np.rad2deg(self.ZD_cbr):.3f} deg')
        fig.tight_layout()
        fig.savefig(figname)

        inds_cal = spec.inds_ton[:, spec.inds_ton.shape[1]//2]
        inds_cal_src, inds_cal_ref = self._sep_src_ref(inds_cal, nB)
        inds_cal_src, inds_cal_ref = self.limit_r_src_ref(inds_cal_src, inds_cal_ref, nB, verbose=False)

        # only use pcal in ref
        ind_tmp = np.where(np.isin(inds_cal, inds_cal_ref))[0]
        pcals = spec._get_cal_power(spec.inds_ton[ind_tmp], spec.inds_toff_bef[ind_tmp], spec.inds_toff_aft[ind_tmp])
        #self.pcals = pcals
        # to do: mitigate rfi effect
        self.pcals_m = pcals_m = np.nanmedian(pcals, axis=0, keepdims=True)
        if self.args.s_method == 'A':
            self.pcals_m_s = pcals_m_s = robust_smooth(pcals_m, 1, self.freq_use, self.args.s_sigma, 1, dlimit=0.016, return_outlier=False)
        elif self.args.s_method == 'B':
            self.pcals_m_s = pcals_m_s = smooth_by_baseline(pcals_m, self.freq_use, self.args.s_sigma, axis=1)
        else:
            raise ValueError(f'`--s_method` should be `A` or `B`, but got `{self.args.s_method}`.')
        # only use noise off spectra
        self.Tcal_s = Tcal_s = spec.get_Tcal_s()
        self.tcal_file = spec.tcal_file
        T_off_src_m = np.nanmedian(p_off_src/pcals_m_s*Tcal_s, axis=0, keepdims=True)
        T_off_ref_m = np.nanmedian(p_off_ref/pcals_m_s*Tcal_s, axis=0, keepdims=True)
        self.Ta = Ta = T_off_src_m - T_off_ref_m

        count_off_src_m = np.nanmedian(p_off_src/pcals_m_s, axis=0, keepdims=True)
        count_off_ref_m = np.nanmedian(p_off_ref/pcals_m_s, axis=0, keepdims=True)
        self.count = count = count_off_src_m - count_off_ref_m

        if self.args.s_method == 'A':
            self.Ta_s = robust_smooth(Ta, 1, self.freq_use, self.args.s_sigma_cbr, 0.5, dlimit=0.015, return_outlier=False)
            self.count_s = robust_smooth(count, 1, self.freq_use, self.args.s_sigma_cbr, 0.5, dlimit=0.015, return_outlier=False)
        elif self.args.s_method == 'B':
            self.Ta_s = smooth_by_baseline(Ta, self.freq_use, self.args.s_sigma_cbr, axis=1)
            self.count_s = smooth_by_baseline(count, self.freq_use, self.args.s_sigma_cbr, axis=1)
        else:
            raise ValueError(f'`--s_method` should be `A` or `B`, but got `{self.args.s_method}`.')


class OnOff(MBC):
    def _sep_src_ref(self, inds_use, nB):
        """
        """
        mjd = self.mjd[inds_use]
        t_src = self.args.t_src/60/60/24 # to day
        t_change = 30/60/60/24 # fixed to 30s for now
        # also assuming t_ref equal to t_src
        t_cir = t_src + t_change
        mjd_splits = self.mjd[0] + np.arange(self.args.n_cir)*t_cir

        is_src,is_ref = np.full(len(mjd), False), np.full(len(mjd), False)
        for st in mjd_splits:
            is_src |= ((mjd > st) & (mjd < st+t_src))
            is_ref |= ((mjd > st+t_src+t_change) & (mjd < st+2*t_src+t_change))

        inds_src, inds_ref = inds_use[is_src], inds_use[is_ref]
        if nB != 1:
            # swap
            inds_src, inds_ref = inds_ref, inds_src
        return inds_src, inds_ref


class Mapping(Calibrator):
    Cal_class = CalOnOffM

    def fit_along_ra(self, T, dra, figname=None, ):
        """
        fitting method of Scan modes
        T: temperature
        dra: ra separation
        top_range: abs(dra) range of fitting the top of curve: peak of the cbr
        bottom_range: abs(dra) range of fitting the bottom of curve: system temperature


        """
        from scipy.optimize import curve_fit
        dra_range_top = self.args.dra_range_top
        dra_range_bottom = self.args.dra_range_bottom
        top_fit_method = self.args.top_fit_method

        dra_abs = np.abs(dra)
        # fit bottom curve
        # limit ra for bottom
        rcoff = (dra_abs >= dra_range_bottom[0]) & (dra_abs <= dra_range_bottom[1])

        tmax = np.max(T)
        tsys = np.mean(T[rcoff])

        offbound= ([-1,tsys-0.1*np.abs(tsys)], [1,tsys+0.1*np.abs(tsys)])
        offp, offcon = curve_fit(Mapping.fit_scoff, dra[rcoff], T[rcoff], bounds=offbound)

        # fit top curve
        # limit ra for top
        rcon = (dra_abs >= dra_range_top[0]) & (dra_abs <= dra_range_top[1])
        # rcon = T >= 0.5*(tmax + offp[1]) # limit dra in half height of T
        # rcon &= dra_abs < 0.1


        if top_fit_method == 'gauss':
            fit_scon = Mapping.fit_scon_gauss
        elif top_fit_method == 'skewgauss':
            fit_scon = Mapping.fit_scon_skewgauss
        else:
            raise ValueError(f'Unknown top_fit_method: {top_fit_method}')
        onbound = ([tmax-0.05*np.abs(tmax)-tsys,-0.01,0,tsys-0.05*np.abs(tsys),],[tmax+0.05*np.abs(tmax)-tsys,0.01,1,tsys+0.05*np.abs(tsys),])
        if top_fit_method == 'skewgauss':
            onbound[0].append(0)
            onbound[1].append(5)
        onp, oncon = curve_fit(fit_scon, dra[rcon], T[rcon], bounds=onbound)

        #ton= onp[0]+onp[3]
        #toff= offp[1]
        #return np.array([ton,toff,tmax,onp[1]])

        # plot fit result
        if figname is not None:
            from matplotlib import pyplot as plt
            plt.figure(figsize=(8,4))
            plt.scatter(dra[rcon], T[rcon], marker='.',color='r')
            plt.scatter(dra[rcoff],T[rcoff], marker='.',color='b')
            dra_plot = np.arange(-0.1, 0.1, 1/3600)
            plt.plot(dra_plot, fit_scon(dra_plot, *onp), color='k')
            plt.plot(dra, Mapping.fit_scoff(dra, *offp), color='k')
            plt.xlabel('RA separation [deg]')
            plt.ylabel('T')
            plt.grid()
            plt.minorticks_on()
            plt.suptitle(os.path.basename(figname))
            plt.tight_layout()
            plt.savefig(figname)

        off_value = offp[1]
        if top_fit_method == 'gauss':
            on_value = onp[0] + onp[3]
        elif top_fit_method == 'skewgauss':
            _x = np.arange(-1/60, 1/60, 0.1/3600)
            _y = fit_scon(_x, *onp)
            _ind = np.argmax(_y)
            on_value = _y[_ind]
            if np.abs(_y[_ind] -_y[_ind+1]) > 1e-3*on_value:
                on_value = None

        return tmax, on_value, off_value, onp, offp

    @staticmethod
    def fit_scon_gauss(x,a,b,c,d):
        return a*np.exp(-(x-b)**2/(2*c**2))+d

    @staticmethod
    def fit_scon_skewgauss(x,a,b,c,d, skew):
        from scipy.stats import skewnorm
        return a*skewnorm.pdf(x, skew, loc=b, scale=c)*(2*np.pi)**0.5*c + d

    @staticmethod
    def fit_scoff(x,k,b):
        return k*x+b

    def _gen_sep_para(self):
        super()._gen_sep_para()
        self.sep_para['freq_step_c'] = 5 #self.freq_use[-1] - self.freq_use[0]
        self.sep_para['pcal_vary_lim_bin'] = 0.02
        self.sep_para['pcal_bad_lim_freq'] = 0.5

    def fit_scan(self, data, dra, nB, figname_add=''):
        """
        fit the scanning data along ra for each frequency bins
        """

        Tuse = data
        rause= dra
        ONv = []
        OFFv = []
        ONp = []
        OFFp = []

        TMAX = []
        freq = self.freq_use
        wbin = 2
        freq_key = fbins = np.arange(freq[0]+wbin/2, freq[-1], wbin/2)
        fit_k = 0

        for i in range(len(freq_key)):
            f_cut = (freq > freq_key[i]-wbin/2) & (freq < freq_key[i]+wbin/2)
            T_fit = np.nanmedian(Tuse[:,f_cut,:], axis=1)
            try:
                if i % (len(freq_key)//3) == 0:
                    figname_XX = self.outname + f'-M{nB:02d}-freq_{freq_key[i]:.2f}-XX-{figname_add}.pdf'
                    figname_YY = self.outname + f'-M{nB:02d}-freq_{freq_key[i]:.2f}-YY-{figname_add}.pdf'
                else:
                    figname_XX = None
                    figname_YY = None
                res_XX  = self.fit_along_ra(T_fit[:,0], rause, figname=figname_XX)
                res_YY  = self.fit_along_ra(T_fit[:,1], rause, figname=figname_YY)
                TMAX.append([res_XX[0], res_YY[0]])
                ONv.append([res_XX[1], res_YY[1]])
                OFFv.append([res_XX[2], res_YY[2]])
                ONp.append([res_XX[3], res_YY[3]])
                OFFp.append([res_XX[4], res_YY[4]])
            except Exception as e:
                print('Fit failed in '+str(int(freq_key[i])) + f'MHz in beam {nB}, please check.')
                print('Error:')
                print_exc()
                TMAX.append([None, None])
                ONv.append([None, None])
                OFFv.append([None, None])
                ONp.append([None, None])
                OFFp.append([None, None])
                fit_k += 1

        if fit_k >= int(0.1*len(freq_key)):
            warnings.warn(f'Too many freq channels failed to fit in beam {nB}. Please check.')
        self.frac_src_used = 1 - fit_k/len(freq_key) # frac of freq channels used

        ONv = np.array(ONv, dtype='float64')
        OFFv = np.array(OFFv, dtype='float64')
        TMAX = np.array(TMAX, dtype='float64')

        Tsrc = interp.interp1d(freq_key, ONv, kind='linear', axis=0, fill_value ='extrapolate')(self.freq_use)
        Tref = interp.interp1d(freq_key, OFFv, kind='linear', axis=0, fill_value ='extrapolate')(self.freq_use)
        TMAX = interp.interp1d(freq_key, TMAX, kind='linear', axis=0, fill_value ='extrapolate')(self.freq_use)

        fit_k = fit_k/len(freq_key)

        return Tsrc, Tref, TMAX, ONp, OFFp, fit_k

    def get_is_fit(self, inds):
        """
        spectra use to fit the Tsys and T_peak
        """
        dra = (self.ra[inds] - self.crd.ra.deg) * np.cos(self.crd.dec.rad)
        ddec = self.dec[inds] - self.crd.dec.deg
        # limit ra and dec
        is_fit = (np.abs(ddec) < self.args.ddec_max) & (np.abs(dra) < self.args.dra_range_bottom[1])
        return is_fit, dra, ddec

    def gen_T(self, nB):
        spec = self.raw_sep
        spec.gen_out_name_base(self.args.outdir)
        spec.plot = True
        spec.plot_pcals = True
        spec.set_para_pcals(calc_diff_method='div',
                            squeeze_diff_freq='median',
                            squeeze_diff_freq_bad_lim=0.5,
                            squeeze_diff_freq_frange=None,
                            method_merge='median',
                            merge_cal_pre_process='scale',
                            method_interp='poly1d',
                            method_interp_edges='none',
                            method_interp_sigma_t=0,
                            )
        spec.prepare_pcals()
        is_fit, dra, _, = self.get_is_fit(spec.inds_off)
        self.ZD_cbr = np.nanmedian(self.ZD[spec.inds_off[is_fit]])
        _, count_fit = spec.get_count_tcal(spec.inds_on[0:0], spec.inds_off[is_fit])
        self.res_count_fit = res_count_fit = self.fit_scan(count_fit, dra[is_fit], nB)

        self.count = count = (res_count_fit[0] - res_count_fit[1])[None,]
        self.Tcal_s = Tcal_s = spec.get_Tcal_s()
        self.tcal_file = spec.tcal_file
        self.Ta = Ta = self.count * Tcal_s
        if self.args.s_method == 'A':
            self.Ta_s = robust_smooth(Ta, 1, self.freq_use, self.args.s_sigma_cbr, 0.5, dlimit=0.015, return_outlier=False)
            self.count_s = robust_smooth(count, 1, self.freq_use, self.args.s_sigma_cbr, 0.5, dlimit=0.015, return_outlier=False)
        elif self.args.s_method == 'B':
            self.Ta_s = smooth_by_baseline(Ta, self.freq_use, self.args.s_sigma_cbr, axis=1)
            self.count_s = smooth_by_baseline(count, self.freq_use, self.args.s_sigma_cbr, axis=1)
        else:
            raise ValueError(f'`--s_method` should be `A` or `B`, but got `{self.args.s_method}`.')


class MappingType2(Mapping):

    """
    fit Cal On and Call off separately
    """
    Cal_class = CalOnOff


    def _gen_sep_para(self):
        super(Mapping, self)._gen_sep_para()


    def gen_T(self, nB):
        args = self.args

        spec = self.raw_sep
        spec.gen_out_name_base(self.args.outdir)
        spec.plot = True

        # fit power of Cal off
        is_fit_off, dra, _, = self.get_is_fit(spec.inds_off)
        self.ZD_cbr = np.nanmedian(self.ZD[spec.inds_off[is_fit_off]])
        p_fit_off = spec.get_field(spec.inds_off[is_fit_off], 'DATA', close_file=False)
        self.res_p_fit_off = res_p_fit_off = self.fit_scan(p_fit_off, dra[is_fit_off], nB, figname_add='-CalOff')

        # fit power of Cal on
        is_fit_on, dra, _, = self.get_is_fit(spec.inds_on)
        #self.ZD_cbr = np.nanmedian(self.ZD[spec.inds_on[is_fit_on]])
        p_fit_on = spec.get_field(spec.inds_on[is_fit_on], 'DATA', close_file=False)
        self.res_p_fit_on = res_p_fit_on = self.fit_scan(p_fit_on, dra[is_fit_on], nB, figname_add='-CalOn')

        # pcal
        #pcal_src = res_p_fit_on[0] - res_p_fit_off[0]
        # use pcal of ref bcz of pcal of src has large noise
        # Cal on is not filled in one sampling time.
        pcal_ref = res_p_fit_on[1] - res_p_fit_off[1]
        self.pcal = pcal_ref
        # if Cal on is not filled in one sampling time
        if hasattr(args, 'Cal_fillratio') and args.Cal_fillratio !=1:
            print(f'The input ``Cal_fillratio`` {args.Cal_fillratio} is used to fix that Cal on is not filled in one sampling time. ')
            self.pcal =  self.pcal/args.Cal_fillratio


        # power of target
        #p_target = (res_p_fit_off[0] + res_p_fit_on[0] - self.pcal_ori)/2.
        p_target = res_p_fit_off[0] - res_p_fit_off[1]
        self.count = count = (p_target / self.pcal)[None,:]
        # Tcal
        self.Tcal_s = Tcal_s = spec.get_Tcal_s()
        self.tcal_file = spec.tcal_file

        self.Ta = Ta = self.count * Tcal_s

        if self.args.s_method == 'A':
            self.Ta_s = robust_smooth(Ta, 1, self.freq_use, self.args.s_sigma_cbr, 0.5, dlimit=0.015, return_outlier=False)
            self.count_s = robust_smooth(count, 1, self.freq_use, self.args.s_sigma_cbr, 0.5, dlimit=0.015, return_outlier=False)
        elif self.args.s_method == 'B':
            self.Ta_s = smooth_by_baseline(Ta, self.freq_use, self.args.s_sigma_cbr, axis=1)
            self.count_s = smooth_by_baseline(count, self.freq_use, self.args.s_sigma_cbr, axis=1)
        else:
            raise ValueError(f'`--s_method` should be `A` or `B`, but got `{self.args.s_method}`.')


def main():
    args = gen_args()
    try:
        parse_obs_log(args)
    except Exception:
        print('Error:')
        print_exc()
        sys.exit(1)
    if args.obsmode == 'MultiBeamCalibration':
        cbr = MBC(args)
    elif args.obsmode == 'OnOff':
        cbr = OnOff(args)
    elif args.obsmode in ['MultiBeamOTF','DriftWithAngle','DecDriftWithAngle','Drift']:
        if args.fit_sep:
            cbr = MappingType2(args)
        elif args.m != args.n:
            cbr = Mapping(args)
        else:
            print('n_on should not equal to n_off in scan modes and `--fit_sep False`, please check.')
            sys.exit(1)
    else:
        print('Unknow obsmode, please check.')
        sys.exit(1)
    cbr()
