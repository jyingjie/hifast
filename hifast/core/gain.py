#!/usr/bin/env python
# coding: utf-8
# Author: Ziming liu zmliu@nao.cas.cn
# Yingjie Jing
import numpy as np
import scipy.interpolate as interp

import pandas as pd

from astropy.coordinates import SkyCoord, EarthLocation
from astropy import coordinates as coord
from astropy import units as u
from astropy.time import Time

import os

from . import conf

def Gain_para():
    gain_para= pd.read_csv(os.path.dirname(__file__)+'/data/FAST_gain_curve.txt', header=None, sep='\s+')
    gain_para.set_index([0,1],inplace=True,)
    gain_para= gain_para[gain_para.columns[::2]]
    gain_para.columns= list(range(1050,1500,50))
    #gain_para.loc['M01','a']
    return gain_para 

def Get_ZA(ra, dec, mjd):
    """
    ra: deg
    dec: deg
    mjd: day
    -------------------
    return ZA: deg
    """
    obs_location = EarthLocation.from_geodetic(lat=conf.lat*u.rad, lon=conf.long*u.rad, height=conf.height*u.m)
    aa_frame =  coord.AltAz(obstime = Time(mjd,format='mjd'), location=obs_location)

    crd= SkyCoord(ra, dec,unit='deg',)
    crd_aa = crd.transform_to(aa_frame)
    za =90 - crd_aa.alt.deg
    return za

def ZA2gain(ZA, nB):
    """
    nB: Beam number
    ZA: zenith angle; deg
    """
    # get eta
    gain_para= Gain_para()
    freq_key= np.array(gain_para.columns)
    
    gain= np.zeros((len(ZA),len(freq_key)))
    is_use= ZA > 26.4
    for i,fre in enumerate(freq_key):
        a,b,c= gain_para.loc[f'M{nB:02d}'][fre][['a','b','c']]
        gain[is_use,i]= c * ZA[is_use] + b + 26.4*(a - c)
        gain[~is_use,i]= a * ZA[~is_use] + b
    # 25.6*eta
    gain *= 25.6
    return gain, freq_key

def Get_gain(ra, dec, mjd, nB, freq=None):
    """
    ra: deg
    dec: deg
    mjd: day
    nB: int
    freq: None or array
    ------------------------
    return gain
           freq_key
    """
    ZA= Get_ZA(ra, dec, mjd)
    gain,freq_key=ZA2gain(ZA, nB)
    if freq is not None:
        gain= interp.interp1d(freq_key,gain,kind='quadratic',fill_value='extrapolate')(freq)
        freq_key=freq
    return gain, freq_key

