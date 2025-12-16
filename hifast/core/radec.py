

__all__ = ['ky_dir_default', 'plot_radec', 'guess_ky', 'process_ky', 'kydata2radec', 'radec_interp', 'get_radec']


import os
import warnings
import re
import sys
from glob import glob
import json

import numpy as np
import scipy.interpolate as interp

import h5py

from astropy import units as u
from astropy.time import Time
from astropy.utils import iers
#iers.Conf.iers_auto_url.set("https://datacenter.iers.org/data/9/finals2000A.all")
from . import kypara2radec
from ..utils.io import rec_his, save_specs_hdf5
import erfa


ky_dir_default = [os.path.expanduser("~")+'/KY/',
                    '/data/inspur_disk06/fast_data/KY/',
                    '/data31/KY/',]


def _tight_ra(ra):
    if len(ra)<=1:
        return ra
    ra_s= np.sort(ra)
    diff_s= np.diff(ra_s)
    ind_max= np.argmax(diff_s)
    if diff_s[ind_max] < ra_s[0]+360-ra_s[-1]:
        return ra
    else:
        ra= np.copy(ra)
        is_c= ra>=ra_s[ind_max+1]
        ra[is_c]= ra[is_c]-360
        return ra


def plot_radec(radec,outname=None,ax=None):
    if ax is None:
        from matplotlib import pyplot as plt
        if outname is not None:
            plt.switch_backend('agg')
        fig, ax = plt.subplots(1,1,figsize=(8,8))
    ind_sort = np.argsort(radec['mjd'])
    for i in range(1,20):
        ax.plot(_tight_ra(radec[f'ra{i}'][ind_sort]), radec[f'dec{i}'][ind_sort])

    ax.set_xlabel('RA')
    ax.set_ylabel('DEC')
    ax.minorticks_on()
    if outname is not None:
        ax.set_title(outname)
        fig.tight_layout()
        fig.savefig(outname)


def guess_ky(fname, obs_mjd, ky_dir=None):
    """
    fname: string
    obs_mjd: float; mjd
    """
    obs_data= (Time(obs_mjd,format='mjd')+8*u.hour).to_datetime().strftime('%Y_%m_%d')

    if ky_dir is None:
        for ky_dir in ky_dir_default:
            print(f'Check if KY directory {ky_dir} exists...')
            if os.path.exists(ky_dir):
                print(f'Search KY .xlsx file in {ky_dir}')
                break
            else:
                print(f'KY directory {ky_dir} does not exist.')
    if ky_dir is None:
        raise(ValueError('No KY directory found'))
    keys= re.match('^.*M[0-1][0-9]_[W,N,F]',os.path.basename(fname)).group(0)[:-6].split('_')
    ky_files=[] #set variable first, in case len(keys)==1
    for i in range(1,len(keys)):
        ky_files= glob(f"{ky_dir}/**/{'_'.join(keys[:-i])}*{obs_data}*.xlsx", recursive=True)
        if len(ky_files)>=1:
            break
    return ky_files


def process_ky(ky_files, mjds_src=None, tol=0, ky_fixed=False):
    '''
    processing KY files, return ky_data, mjd, ky_file

    ky_files: list; path of feed (馈源舱) recording file (.xlsx)
    mjds_src: array; if not None, will check if the feed file cover the mjds_src
    tol: float; tolerable time, unit is second
    '''
    tol = float(tol)/24/3600 # second to day
    for ky_file in ky_files:
        print(f"check {ky_file}")
        this = True
        utcoffset = 8*u.hour
        #read ky data
        sheet_name = '整控-馈源舱数据'
        import openpyxl
        wb = openpyxl.load_workbook(ky_file, data_only=True, read_only=True)
        # Select sheet
        if sheet_name in wb.sheetnames:
             ws = wb[sheet_name]
        else:
             print("can't find '整控-馈源舱数据', try to use the first sheet")
             ws = wb.active

        # Read headers
        rows = ws.rows
        headers = [cell.value for cell in next(rows)]
        # Filter None headers
        headers = [h if h else f"Unnamed:{i}" for i, h in enumerate(headers)] # Handle empty headers if any? 
        # But ky_data usually has nice headers.

        data_dict = {h: [] for h in headers}
        
        count = 0
        for row in rows:
             count += 1
             for i, cell in enumerate(row):
                 if i < len(headers):
                     data_dict[headers[i]].append(cell.value)
                     
        if count == 0:
            print('empty sheet')
            continue
            
        ky_data = data_dict
        systime = np.array(ky_data['SysTime'], dtype=str)
        systime = Time(systime, format='iso', scale='utc') - utcoffset # covert local time to UTC
        mjd = systime.mjd
        ky_data['systime_utc'] = systime
        if mjds_src is None:
            # if only deal with ky file
            return ky_data, mjd, ky_file

        diff = mjds_src.min()-mjd.min()
        if -tol < diff < 0:
            if not ky_fixed:
                print(f"The earliest input mjd is earlier by {abs(diff)*24*3600:.3f}s than the mjd in KY. Which is in the set tolerable time, RA and Dec will be extrapolated.")
            else:
                print(f"The earliest input mjd is earlier by {abs(diff)*24*3600:.3f}s than the mjd in KY.")
        if -tol > diff:
            this = False
            print(f"The earliest input mjd is earlier by {abs(diff)*24*3600:.3f}s than the mjd in KY.")
            #raise(ValueError(f"input mjd exceed {abs(diff)*24*3600}s to the mjd in KY,"))
        if ky_fixed:
            diff = mjd.max()-mjds_src.min()
            if -tol < diff< 0:
                print(f"the begin of input mjd is later {abs(diff)*24*3600:.3f}s than the end of mjd in KY. ")
            if -tol > diff:
                this = False
                print(f"the begin of input mjd is later {abs(diff)*24*3600:.3f}s than the end of mjd in KY. Abort.")
                #raise(ValueError(f"input mjd exceed {abs(diff)*24*3600}s to the mjd in KY,"))

        diff = mjd.max() - mjds_src.max()
        if -tol < diff< 0:
            if not ky_fixed:
                print(f"The last input mjd is later by {abs(diff)*24*3600:.3f}s than the mjd in KY. Which is in the set tolerable time, RA and Dec will be extrapolated.")
            else:
                print(f"The last input mjd is later by {abs(diff)*24*3600:.3f}s than the mjd in KY.")
        if -tol > diff:
            #raise(ValueError(f"input mjd exceed {abs(diff)*24*3600}s to the mjd in KY,"))
            if not ky_fixed:
                this = False
                print(f"The last input mjd is later by {abs(diff)*24*3600:.3f}s than the mjd in KY. Abort.")

        if this:
            break

    if not this:
        raise(ValueError(f"The input mjd are not in the mjd range of those KY files"))
    else:
        print(f"Using {ky_file}")
        return ky_data, mjd, ky_file


def kydata2radec(ky_data, mjd, nBs='All', ky_fixed=False, nproc=1, backend='astropy'):
    """
    calculating ra dec of ky from ky data
    ky_data: data in feed file (.xlsx)
    mjd: mjd of ky_data if ky_fixed is False, else mjd of the source
    nBs: list or str:'All'
    """
    if nBs== 'All':
        nBs = range(1,20)

    multibeamAngle= np.asarray(ky_data['SDP_AngleM'], dtype='float64')

    #实测中心波束相对中心的全局坐标
    globalCenterX = np.asarray(ky_data['SDP_PhaPos_X'], dtype='float64')
    globalCenterY = np.asarray(ky_data['SDP_PhaPos_Y'], dtype='float64')
    globalCenterZ = np.asarray(ky_data['SDP_PhaPos_Z'], dtype='float64')
    #实测下平台的全局姿态角
    globalYaw = np.asarray(ky_data['SDP_SwtDPose_Y'], dtype='float64')
    globalPitch = np.asarray(ky_data['SDP_SwtDPose_P'], dtype='float64')
    globalRoll = np.asarray(ky_data['SDP_SwtDPose_R'], dtype='float64')

    if ky_fixed:
        warnings.warn("By using the ky_fixed parameter, the feed is assumed to be stationary during observation.")
        std= np.std(multibeamAngle[-10:])/np.pi*180
        if std > 3:
            raise(ValueError(f'multibeamAngle variation at last 10 points is {std} degree.'))
        full_fun= lambda x:np.full(len(mjd), np.mean(x[-5:]))
        multibeamAngle= full_fun(multibeamAngle)
        globalCenterX = full_fun(globalCenterX)
        globalCenterY = full_fun(globalCenterY)
        globalCenterZ = full_fun(globalCenterZ)
        globalYaw =     full_fun(globalYaw)
        globalPitch =   full_fun(globalPitch)
        globalRoll =    full_fun(globalRoll)

    radec= {}
    radec['mjd']= mjd
    radec['angle'] = multibeamAngle
    #ymdhms = ky_data['systime_utc'].ymdhms
    #utc1, utc2 = erfa.dtf2d(b"UTC", ymdhms['year'], ymdhms['month'], ymdhms['day'], ymdhms['hour'], ymdhms['minute'], ymdhms['second'])
    if ky_fixed:
        obstime = Time(mjd, format='mjd', scale='utc')
    else:
        obstime = ky_data['systime_utc']

    if nproc > 1:
        import multiprocessing
        global fun_mp
        def fun_mp(para):
            nB, verbose = para
            print(f'beam {nB}')
            return kypara2radec.kypara2radec(obstime, multibeamAngle, nB,
                                              globalCenterX,  globalCenterY, globalCenterZ, globalYaw, globalPitch, globalRoll,
                                               backend=backend, verbose=verbose)
        with multiprocessing.Pool(processes=nproc) as pool:
            res = pool.map(fun_mp, zip(nBs, [True,] + [False, ]*(len(nBs)-1)))
        for i, nB in enumerate(nBs):
            radec['ra'+str(nB)]= res[i][0]
            radec['dec'+str(nB)]= res[i][1]
            radec['Az'+str(nB)]= res[i][2]
            radec['ZD'+str(nB)]= res[i][3]
    else:
        for ii, nB in enumerate(nBs):
            verbose = True if ii == 0 else False
            ra_, dec_, Az_, ZD_ = kypara2radec.kypara2radec(obstime, multibeamAngle, nB,
                                          globalCenterX,  globalCenterY, globalCenterZ, globalYaw, globalPitch, globalRoll, backend=backend, verbose=verbose)
            radec['ra'+str(nB)]= ra_
            radec['dec'+str(nB)]= dec_
            radec['Az'+str(nB)]= Az_
            radec['ZD'+str(nB)]= ZD_
            print(f'beam {nB}')
    return radec


def radec_interp(radec_ky, mjds_src):
    """
    get ra dec of source through interpolating that of ky
    """
    mjds_ky= radec_ky['mjd']
    radec_src={}
    radec_src['mjd']= mjds_src
    print('interpolating...')
    radec_src['is_extrapo']= (mjds_src < np.nanmin(mjds_ky)) | (mjds_src > np.nanmax(mjds_ky))
    for key in radec_ky:
        if 'ra' in key or 'dec' in key or 'angle' in key or 'Az' in key or 'ZD' in key:
            value = radec_ky[key]
            if 'ra' in key:
                value = _tight_ra(value)
            radec_src[key]= interp.interp1d(mjds_ky, value, kind='linear', fill_value ='extrapolate')(mjds_src)
            if 'ra' in key:
                radec_src[key][radec_src[key]<0] += 360 # PyAstronomy.pyasl don't support negative ra
    return radec_src


def get_radec(parse_mjd, guess_str=None, ky_files=None, tol=0.3, ky_fixed=False, use_cache=True, nproc=1,
              backend='astropy', env_para=None, dUT1=None, nBs='All', cache_kyradec=True):
    """
    parse_mjd: str or array
               str: *.hdf5 or *.xlsx
               array: shape is (m,); need input guess_str or ky_files
    guess_str: str; used to guess feed file name, for example, "M31_Drift_v3_6_arcdrift-M19_W"
               if None, use parse_mjd instead and parse_mjd must be str.
    ky_files: list of str; specify feed files; if None, will guess from guess_str or parse_mjd (str).
    tol: float; unit is second
    ky_fixed: bool; if True, assuming the feed is stationary during observation.
    use_cache: bool; use cached radec of feed
    nproc: int; parallel cpu number
    backend: str; erfa or astropy
    env_para: dict; e.g.
            env_para = {}
            env_para['phpa'] = 925.
            env_para['temperature'] = 15
            env_para['humidity'] = 0.8
    dUT1: UT1 - UTC
    """
    if nBs== 'All':
        nBs = range(1,20)
    elif isinstance(nBs, list):
        use_cache = False
        cache_kyradec = False
    else:
        raise(ValueError("nBs need be 'All' or a list with numbers"))

    ident = {}
    if env_para is not None:
        kypara2radec.set_obs_env_para(env_para['phpa'], env_para['temperature'], env_para['humidity'])
        ident.update(env_para)
    if dUT1 is not None:
        kypara2radec.dUT1 = dUT1
        if backend == 'erfa':
            ident['dUT1'] = dUT1
    import json
    from hashlib import md5
    ident_str = md5(json.dumps(ident).encode('utf-8')).hexdigest()

    # if input feed file, return ra dec in it
    if isinstance(parse_mjd, str):
        fname = parse_mjd
        if 'xls' in fname.split('.')[-1]:
            ky_data, mjd, ky_file = process_ky([fname,])
            return kydata2radec(ky_data, mjd, nBs=nBs, ky_fixed=False, backend=backend, nproc=nproc)
        elif 'hdf5' in fname.split('.')[-1]:
            with h5py.File(fname,'r') as f:
                    mjds_src = f['S']['mjd'][:]
        else:
            raise(ValueError('fname needs end with xlsx or hdf5'))
    else:
        mjds_src = parse_mjd
        if guess_str is None and ky_files is None:
            raise(ValueError('need input `guess_str` or `ky_files`'))
    if guess_str is None and isinstance(parse_mjd, str):
        guess_str = parse_mjd


    # obtain cache_file and cache_key
    cache_dir= os.path.expanduser('~/.cache/fast_python/kyradec/')
    os.makedirs(cache_dir, exist_ok=True)
    cache_file = cache_dir + 'known_kyfile.json'
    obs_data = (Time(mjds_src[0], format='mjd')+8*u.hour).to_datetime().strftime('%Y_%m_%d')
    try:
        cache_key = re.match('^.*M[0-1][0-9]_[W,N,F]', os.path.basename(guess_str)).group(0)[:-6] + '-' + obs_data
    except:
        cache_key = None
        use_cache = False

    # find feed file corresponding to guess_str
    ky_file = None
    if ky_files is None:
        # try to find feed file corresponding to guess_str
        if use_cache:
            #if cached file is not exits, init
            if not os.path.exists(cache_file):
                    with open(cache_file, "w") as f:
                        json.dump({}, f)
            try:
                with open(cache_file,'r') as f:
                    ky_file= json.load(f)[cache_key]
                print('get kyfile name from ', cache_file)
                #print('using ', ky_file)
                ky_data = None
            except KeyError:
                pass
    # if can't get ky_file from cache
    if ky_file is None:
        if ky_files is None:
            obs_mjd = mjds_src[0]
            ky_files = guess_ky(guess_str, obs_mjd)
            if len(ky_files) == 0:
                raise(ValueError("can not find KY file similar to input file, please appoint KY files" ))
    else:
        ky_files = [ky_file,]
    # test kyfiles and load data
    res = process_ky(ky_files, mjds_src, tol, ky_fixed=ky_fixed)
    if res is None:
        return
    else:
        ky_data, mjds_ky, ky_file = res
    # cache ky_file
    if cache_key is not None and cache_kyradec:
        import tempfile
        file_tmp = tempfile.mktemp(prefix=os.path.basename(cache_file)+'.tmp.', dir=os.path.dirname(cache_file))
        with open(cache_file, "r") as f:
            cache = json.load(f)
        cache.update({cache_key:ky_file})
        with open(file_tmp,'w') as f:
            json.dump(cache, f, indent=4)
            f.flush()
            os.fsync(f.fileno())
        os.rename(file_tmp, cache_file)
    # obtain radec of ky
    cache_fname_radec = os.path.join(cache_dir,'.'.join(os.path.basename(ky_file).split('.')[:-1]) + f'_{backend}_cache_{ident_str}.npz')
    cache_fname_radec_b = os.path.join(cache_dir,'.'.join(os.path.basename(ky_file).split('.')[:-1]) + f'_{backend}_cache_{ident_str}.npy')
    radec_ky = None
    if use_cache and not ky_fixed:
        try:
            radec_ky = dict(np.load(cache_fname_radec).items())
            print('using cached ra dec of KY: ', cache_fname_radec)
        except:
            try:
                radec_ky = np.load(cache_fname_radec_b, allow_pickle=True)[()]
                print('using cached ra dec of KY: ', cache_fname_radec_b)
            except:
                pass
    if radec_ky is None:
        if ky_data is None:
            # ky_file is obtained from cache_file
            ky_data, mjds_ky, ky_file= process_ky([ky_file,], mjds_src, tol, ky_fixed=ky_fixed)
        if ky_fixed:
            return kydata2radec(ky_data, mjds_src, ky_fixed=ky_fixed, nproc=nproc, backend=backend, nBs=nBs)
        else:
            radec_ky= kydata2radec(ky_data, mjds_ky, ky_fixed=ky_fixed, nproc=nproc, backend=backend, nBs=nBs)
            if cache_kyradec:
                print('saving radec of KY to cached file')
                np.savez(cache_fname_radec, **radec_ky)
    return radec_interp(radec_ky, mjds_src)
