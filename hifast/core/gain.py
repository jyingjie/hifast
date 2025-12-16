#!/usr/bin/env python
# coding: utf-8
# Author: Ziming liu zmliu@nao.cas.cn
# Yingjie Jing
import numpy as np
import scipy.interpolate as interp



from astropy.coordinates import SkyCoord, EarthLocation
from astropy import coordinates as coord
from astropy import units as u
from astropy.time import Time

import os

import json

from . import conf

def Gain_para():
    fpath = os.path.dirname(__file__)+'/data/FAST_gain_curve.txt'
    # Manually parse the file to replicate:
    # gain_para.set_index([0,1],inplace=True)
    # gain_para = gain_para[gain_para.columns[::2]] (taking every 2nd column?)
    # original code: gain_para= pd.read_csv(..., header=None, sep='\s+')
    # gain_para.set_index([0,1]) -> first two cols are index: Beam, Param
    # gain_para= gain_para[gain_para.columns[::2]] -> if original cols were 0,1,2,3,4... after set_index, cols are 2,3,4...
    # The ::2 on columns likely means skipping error columns or similar.
    # gain_para.columns = list(range(1050, 1500, 50))
    
    # Structure needed: data[beam][freq][param] or data[beam][param][freq]
    # Existing usage: gain_para.loc[f'M{nB:02d}'][fre][['a','b','c']]
    # This implies: result[beam][freq] -> gives a series/dict where keys are 'a','b','c'
    # Wait, gain_para.loc['M01'] gives a DF indexed by 'a','b','c'.
    # gain_para.loc['M01'][fre] gives the column at 'fre', which is a Series with index 'a','b','c'.
    # So we need: data[beam][freq] = {'a': val, 'b': val, 'c': val}
    
    data = {}
    with open(fpath, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if not parts: continue
            beam = parts[0]
            param = parts[1]
            values = parts[2:]
            
            # Original: gain_para = gain_para[gain_para.columns[::2]]
            # If values indices are 0,1,2,3... ::2 takes 0, 2, 4...
            selected_values = values[::2]
            
            freqs = list(range(1050, 1500, 50))
            if len(selected_values) != len(freqs):
                 # Fallback or error? Assuming file matches expected length
                 pass
            
            if beam not in data: data[beam] = {}
            
            for i, freq in enumerate(freqs):
                if freq not in data[beam]: data[beam][freq] = {}
                data[beam][freq][param] = float(selected_values[i])
                
    # To support syntax: gain_para.loc['M01'][fre]['a']
    # We can wrap it in a class or just return a dict-like object.
    # But existing code: gain_para.loc[f'M{nB:02d}'][fre][['a','b','c']]
    # If data is dict: data['M01'][fre] returns {'a':..., 'b':..., 'c':...}
    # Dicts don't support list indexing [['a','b','c']].
    # We need to change the usage in `ZA2gain` as well.
    return data 

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
    # freq_key= np.array(gain_para.columns) 
    # Use keys from first available beam (e.g. M01)
    # The keys are frequencies.
    # gain_para structure: data[beam][freq][param]
    first_beam = list(gain_para.keys())[0]
    freq_key = np.array(list(gain_para[first_beam].keys()))
    
    gain= np.zeros((len(ZA),len(freq_key)))
    is_use= ZA > 26.4
    for i,fre in enumerate(freq_key):
        # a,b,c= gain_para.loc[f'M{nB:02d}'][fre][['a','b','c']]
        params = gain_para[f'M{nB:02d}'][fre]
        a, b, c = params['a'], params['b'], params['c']
        
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