## cal
import numpy as np

def read_tcal_sav(nB, s_type='w', tcal_dir=None, mode='high', date='20190115'):
    """
    nB: int
      beam number
    """
    from scipy.io.idl import readsav
    if tcal_dir is None:
        tcal_dir= os.path.expanduser("~")+'/Tcal/'
    fname = tcal_dir+f'{date}/median_{date}.Tcal-results.HI_{s_type}.{mode}.sav'
    tc_info= readsav(fname)[f'{mode}_{s_type}'][0]
    tc_freq= tc_info['freq']
    return tc_freq, tc_info['M%02d_TC'%nB], fname

def read_tcal_fits(nB, s_type='w', tcal_dir=None, mode='high', date=''):
    """
    nB: int
      beam number
    """
    from astropy.io import fits
    import os
    if tcal_dir is None:
        tcal_dir= os.path.expanduser("~")+'/Tcal/'
    fname = tcal_dir + f'{date}/CAL.{date}.{mode}.{s_type.upper()}.fits'
    f = fits.open(fname)
    tc_freq = f[1].data['FREQ'][0]
    tc_T = f[1].data['TCAL'][0,nB-1].T
    return tc_freq, tc_T, fname

def read_tcal(nB, s_type='w', tcal_dir=None, mode='high', date='auto', mjd=None):
    """
    nB: int
       beam number
    s_type: str
       type, w or n
    tcal_dir: str
       if not set, use ~/Tcal/
    mode: str
       high or low
    date: str
       example: 20190115 or 20200531, default is None
    mjd:
       if date is 'auto', using the nearest date of tcal
    """
    from glob import glob
    import os
    if tcal_dir is None:
        tcal_dir= os.path.expanduser("~")+'/Tcal/'
    dates_have = [os.path.basename(i) for i in glob(tcal_dir + '/20[0-9][0-9][0-9][0-9][0-9][0-9]')]
    if len(dates_have) == 0:
        raise(ValueError('can not find tcal file'))
    if date == 'auto':
        from astropy.time import Time
        if mjd is None:
            raise(ValueError('need input mjd if date is auto'))
        mjds_h = Time([f'{s[:4]}-{s[4:6]}-{s[6:8]} 00:00:00.000' for s in dates_have], format='iso').mjd
        #print(mjds_h)
        date = dates_have[np.argmin(abs(mjds_h - mjd))]
    else:
        if date not in dates_have:
            raise(ValueError(f'can not find tcal file in {date}'))
    ## read tcal
    if date == '20190115':
        return read_tcal_sav(nB, s_type, tcal_dir, mode, date='20190115')
    else:
        return read_tcal_fits(nB, s_type, tcal_dir, mode, date=date)
