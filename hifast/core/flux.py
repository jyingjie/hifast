import h5py
import numpy as np
import scipy.interpolate as interp
import pandas as pd
import os
from os.path import basename
import re

from .gain import Get_gain

def get_ratio(nB, freq=None):
    ratios_para= pd.read_json(os.path.dirname(__file__)+'/data/beam_ratios.json')
    #ratios_para= pd.read_json('beam_ratios.json')
    freq_key= np.array(ratios_para.index[::2],dtype=float)
    ratios= ratios_para[f'M{nB:02d}'][::2].values
    if freq is not None:
        #ratios= interp.interp1d(freq_key,ratios, kind='nearest', fill_value= "extrapolate")(freq)
        ratios= interp.interp1d(freq_key,ratios, kind='quadratic', fill_value= "extrapolate")(freq)
        freq_key=freq
    return ratios, freq_key

def get_K_Jy_cali(cali_fname, nB, freq, ra=None, dec=None, mjd=None, tcal_spec=None):
    """
    K/Jy from Calibator
    ----------------------
    cali_fname: quasar calibration file name
    nB: source beam number
    freq: 
    ra, dec, mjd: properties of spectra, needed if only Beam 1 in cali_fname
    tcal_spec: 
    """
    
    with h5py.File(cali_fname,'r') as fs:
        freq_c= fs['freq'][()]
        if f'M{nB:02d}' in fs.keys():
            K_Jy = fs[f'M{nB:02d}'][0] # K/Jy
            if tcal_spec is not None:
                K_Jy /= fs[f'Tcal{nB}'][0] # to count of tcal
            need_ratio = False
        else:
            K_Jy = fs[f'M01'][()][0] # K/Jy
            if tcal_spec is not None:
                K_Jy /= fs[f'Tcal1'][0] # to count of tcal
            need_ratio = True 
    K_Jy = interp.interp1d(freq_c, K_Jy, kind='quadratic', fill_value= "extrapolate", axis=0)(freq)
    if tcal_spec is not None:
        K_Jy *= tcal_spec[0]
    K_Jy = K_Jy[None,...] # mjd axis
    if K_Jy.ndim == 2:
        K_Jy = K_Jy[..., None] # polar axis
    if need_ratio:
        K_Jy = K_Jy*get_ratio(nB, freq)[0][None,:,None]
    return K_Jy

def cali_src(T, nB, freq, cali_fname=None, ra=None, dec=None, mjd=None, tcal_spec=None):
    """
    flux calibration using fixed factor or Calibator 
    -----------------------
    T: array_like
       Temperature of the spectra. Shape is (m,n) or (m,n,2) i.e. (Mjd, channel) or (Mjd, channel, Polarization)
    nB: int
       Beam numbe
    freq: array_like, shape (n,)
    cali_fname: str
       quasar calibration file name; hdf5 file
       If None, use the gain depended on Zenith angle (arxiv:2002.01786) and need input ra, dec and mjd.
    ra, dec, mjd: None or array_like, shape (m,)
    tcal_spec: shape: (1,n,2); T.shape need be (m,n,2); if spec and Calibator used different Tcal, input this to fix it.
    """
    
    if cali_fname is None or cali_fname=='none':
        #check
        if np.isscalar(mjd):
            mjd = np.array([mjd])
        if np.isscalar(ra):
            ra = np.array([ra])
        if np.isscalar(dec):
            dec = np.array([dec])
        assert ra.shape[0] == dec.shape[0] == mjd.shape[0] == T.shape[0]   
        
        K_Jy = Get_gain(ra, dec, mjd, nB, freq)[0] * 25.6   #K/Jy
        K_Jy= K_Jy[:,:,None]
        
    else:
        if tcal_spec is not None:
            if T.ndim != 3 or T.shape[-1] !=2:
                raise(ValueError('tcal_spec is not None, T.shape need be (m,n,2)'))
            if tcal_spec.ndim != 3 or tcal_spec.shape[0] !=1 or tcal_spec.shape[-1] !=2:
                raise(ValueError('tcal_spec.shape need be (1, n, 2)'))
        K_Jy = get_K_Jy_cali(cali_fname, nB, freq, tcal_spec=tcal_spec)
    
    if T.ndim == 2:
        return T/np.mean(K_Jy, axis=2)
    elif T.ndim == 3:
        if T.shape[-1] == 2:
            return T/K_Jy
        elif T.shape[-1] == 1:
            return T/np.mean(K_Jy, axis=2, keepdims=True)
    else:
        raise(ValueError('shape of T'))

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(allow_abbrev=False)
    parser.add_argument('fname',
                        help='file name')
    parser.add_argument('-c', '--cali_fname',
                       help='quasar calibration file name')
    parser.add_argument('--outdir',
                       help='default is same with the input file')

    args = parser.parse_args()
    fname=args.fname
    outdir=args.outdir
    
    cali_fname=args.cali_fname
    nB= int(re.findall(r'-M[0-1][0-9]',basename(fname))[0][2:])

    f= h5py.File(fname,'r')
    T= f['Ta'][()]
    freq= f['freq'][()]
    mjd= f['mjd'][()]
    ra= f['ra'][()]
    dec= f['dec'][()]
    if 'Header' in f.keys():
        from collections import OrderedDict
        header_in= OrderedDict(f['Header'].attrs.items())
    else:
        header_in=None
    f.close()
    flux= cali_src(T, nB, freq, cali_fname, ra=ra,dec=dec,mjd=mjd)
    
    #save
    if outdir is None:
        outdir= os.path.dirname(fname)
    fileout= os.path.join(outdir, '.'.join(os.path.basename(fname).split('.')[:-1]) +'-flux.hdf5' )

    print(f"saving to {fileout}")
    dict_out= {}
    dict_out['freq']=freq
    dict_out['ra']=ra
    dict_out['dec']=dec
    dict_out['mjd']=mjd
    dict_out['flux']=flux.astype('float32')
    #save file
    from ..utils.io import rec_his,save_dict_hdf5
    header=rec_his(args=args)
    if header_in is not None: header.update(header_in)
    save_dict_hdf5(fileout, dict_out, header=header)
