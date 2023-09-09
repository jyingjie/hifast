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

import json

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


def gain_diff_from_ZA(ZA_1: float, ZA_2: np.array, nB: int, freq: np.array) -> np.array:
    """
    Calculate the gain difference between two different ZA's (zenith angle)

    :param ZA_1: scalar value for first ZA
    :param ZA_2: 1d array values for second ZA
    :param nB: scalar integer 
    :param freq: 1d array values for frequency
   
    :return: array of gain differences for each frequency. (ZA_2 - ZA_1)
    """
    # eta diff
    
    fpath = os.path.dirname(__file__) + '/data/gain_ZA_fit.json'
    
    # Open and load the json data into a dictionary
    with open(fpath,'r') as f:
        coeffs_dict = json.load(f)

    fun_a = np.poly1d(coeffs_dict['a'][f'M{nB:02}'])
    fun_c = np.poly1d(coeffs_dict['c'][f'M{nB:02}'])

    # Calculate the delta between two ZA's
    delta = (ZA_2 - ZA_1)[:, None]
    diffs = np.zeros((len(delta), len(freq)), dtype=float)
    
    # Check for ZA_2 angles below 26.4
    is_s = ZA_2 <= 26.4

    # Compute diffs based on the whether the ZA's are below or above 26.4
    if ZA_1 <= 26.4:
        diffs[is_s] = delta[is_s] * fun_a(freq)
        diffs[~is_s] = ((ZA_2[~is_s] - 26.4)[:, None] * fun_c(freq) +
                         (26.4 - ZA_1) * fun_a(freq))
    elif ZA_1 > 26.4:
        diffs[~is_s] = delta[~is_s] * fun_c(freq)
        diffs[is_s] = ((ZA_2[is_s] - 26.4)[:, None] * fun_a(freq) +
                         (26.4 - ZA_1) * fun_c(freq))
    #  25.6*eta
    diffs *= 25.6
    return diffs