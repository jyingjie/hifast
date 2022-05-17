#!/usr/bin/env python
# coding: utf-8
# Adapted from FAST_ResultDataInversTool.exe

import numpy as np
import math
import erfa

from . import conf

#观测波长
wl = 300000.0;
gap = 0.27 #多波束中心间距 米

#属性-可更新的参数：
dUT1 =  0.1 # UT1-UTC is tabulated in IERS bulletins.
phpa = 925. #气压
temperature =  15. #气温
humidity =  0.8 #相对湿度

def set_obs_env_para(phpa_new=None, temperature_new=None, humidity_new=None, reset=False):
    
    global phpa #气压
    global temperature #气温
    global humidity #相对湿度
    
    if reset:
        phpa = 925. #气压
        temperature =  15. #气温
        humidity =  0.8 #相对湿度
    else:
        if phpa_new is not None:
            phpa = phpa_new
        if temperature_new is not None:
            temperature = temperature_new
        if humidity_new is not None:
            humidity = humidity_new
    
feedLocalCoord = [  [ 0.0  , 0.0  , 0.0   ],
                    [ -0.27, 0.0, 0.0],
                    [ -0.135, 0.233826859, 0.0],
                    [ 0.135, 0.233826859, 0.0],
                    [ 0.27, 0.0, 0.0],
                    [ 0.135, -0.233826859, 0.0],
                    [ -0.135, -0.233826859, 0.0],
                    [ -0.54, 0.0, 0.0],
                    [ -0.405, 0.233826859, 0.0],
                    [ -0.27, 0.467653718, 0.0],
                    [ 0.0, 0.467653718, 0.0],
                    [ 0.27, 0.467653718, 0.0],
                    [ 0.405, 0.233826859, 0.0],
                    [ 0.54, 0.0, 0.0],
                    [ 0.405, -0.233826859, 0.0],
                    [ 0.27, -0.467653718, 0.0],
                    [ 0.0, -0.467653718, 0.0],
                    [ -0.27, -0.467653718, 0.0],
                    [ -0.405, -0.233826859, 0.0] ]


def kypara2radec(obstime, multibeamAngle, nB, globalCenterX,  globalCenterY, globalCenterZ, globalYaw, globalPitch, globalRoll, backend='erfa', verbose=True):
    """
    array except nB
    obstime: observation time, astropy.time.Time
    """
    
    Az, ZD = kypara2AzZD_ufun(multibeamAngle, nB, globalCenterX,  globalCenterY, globalCenterZ, globalYaw, globalPitch, globalRoll)
    Az = np.array(Az, dtype='float64')
    ZD = np.array(ZD, dtype='float64')
    
    if verbose:
        print(f"Use {backend} to convert Az, ZD to RA, DEC")
        print(f"atmospheric pressure: {phpa} hPa")
        print(f"ground-level temperature: {temperature} deg C")
        print(f"relative humidity: {humidity}")
        if backend == 'erfa':
            print(f"dUT1: {dUT1}")
    ra, dec = AzZD2radec(obstime, Az, ZD, backend=backend)
    return ra, dec, Az, ZD
    
    
    
def kypara2AzZD(multibeamAngle, nB, globalCenterX,  globalCenterY, globalCenterZ, globalYaw, globalPitch, globalRoll):
    """
    return Az: Azimuth in radian, ZD: zenith distance in radian
    """
    #1.多波束转角带来的旋转
    rotationMatrixMultiBeam= CalMultiBeamRotationMatrix(multibeamAngle)

    #2.Stewart下平台带来的旋转、
    rotationMatrixPlatform= CalPlatformRotationMatrix(globalYaw, globalPitch, globalRoll)

    #3.计算旋转后的位置
    useFeed = feedLocalCoord[nB-1]
    posRelative = VectorTransform(rotationMatrixPlatform, VectorTransform(rotationMatrixMultiBeam, useFeed))

    #4.计算全局坐标
    posAbsolute = posRelative + np.array([globalCenterX,  globalCenterY, globalCenterZ])

    #5.计算地平坐标
    #Az
    Az = math.atan2(-posAbsolute[0], -posAbsolute[1])
    #ZD
    ZD = math.atan2(math.sqrt(posAbsolute[0]**2 + posAbsolute[1] **2), -posAbsolute[2])
    
    return Az, ZD

kypara2AzZD_ufun = np.frompyfunc(kypara2AzZD, 8, 2)

def AzZD2radec(obstime, Az, ZD, backend='erfa'):
    """
    obstime: observation time, astropy.time.Time
    Az: array, Azimuth in radian
    ZD: array, zenith distance in radian
    backend: str, 'astropy' or 'erfa'
    """
    if backend == 'erfa':
        #     The polar motion xp,yp can be obtained from IERS bulletins.  The
        #         values are the coordinates (in radians) of the Celestial
        #         Intermediate Pole with respect to the International Terrestrial
        #         Reference System (see IERS Conventions 2003), measured along the
        #         meridians 0 and 90 deg west respectively.  For many
        #         applications, xp and yp can be set to zero.
        
#         If hm, the height above the ellipsoid of the observing station
#         in meters, is not known but phpa, the pressure in hPa (=mB), is
#         available, an adequate estimate of hm can be obtained from the
#         expression

#               hm = -29.3 * tsl * log ( phpa / 1013.25 );

#         where tsl is the approximate sea-level air temperature in K
#         (See Astrophysical Quantities, C.W.Allen, 3rd edition, section
#         52).  Similarly, if the pressure phpa is not known, it can be
#         estimated from the height of the observing station, hm, as
#         follows:

#               phpa = 1013.25 * exp ( -hm / ( 29.3 * tsl ) );

#         ***Note, however, that the refraction is nearly proportional to
#         the pressure and that an accurate phpa value is important for
#         precise work.***
        #from astropy.coordinates.builtin_frames.utils import get_polar_motion
        xp = 0.
        yp = 0.
        
        if dUT1 == 'auto':
            dUT1_this = obstime.delta_ut1_utc
        else:
            dUT1_this = dUT1
            
        ra, dec = erfa.atoc13(b"A", Az, ZD, obstime.jd1, obstime.jd2, dUT1_this, elong, phi, hm, xp, yp, phpa, temperature, humidity, wl)
        ra = np.rad2deg(ra)
        dec = np.rad2deg(dec)
        
    elif backend == 'astropy':
        from astropy.coordinates import SkyCoord, EarthLocation
        from astropy.time import Time
        import astropy.units as u
        
        location = EarthLocation(lat=conf.lat*u.radian, lon=conf.long*u.radian, height=conf.height*u.m)
        
        coords = SkyCoord(frame='altaz', alt=90*u.deg-ZD*u.radian, az=Az*u.radian, obstime=obstime, location=location,
                pressure=phpa*u.hPa, temperature=temperature*u.deg_C, relative_humidity=humidity, obswl=wl*u.um)
        coords_icrs = coords.icrs
        ra, dec = coords_icrs.ra.value, coords_icrs.dec.value
        
    return ra, dec


#other functions
#从多波束转角计算旋转矩阵
def CalMultiBeamRotationMatrix(multibeamAngle):          
    """
    multibeamAngle: scalar  float64  Radian 
    """    
    ca= math.cos(multibeamAngle) 
    sa= math.sin(multibeamAngle) 
    rotationMatrixMultiBeam= np.zeros((3,3), dtype='float64')
    
    rotationMatrixMultiBeam[0,0] = ca 
    rotationMatrixMultiBeam[0,1] = -sa 
    rotationMatrixMultiBeam[0,2] = 0.0 
    rotationMatrixMultiBeam[1,0] = sa 
    rotationMatrixMultiBeam[1,1] = ca 
    rotationMatrixMultiBeam[1,2] = 0.0 
    rotationMatrixMultiBeam[2,0] = 0.0 
    rotationMatrixMultiBeam[2,1] = 0.0 
    rotationMatrixMultiBeam[2,2] = 1.0 
    return rotationMatrixMultiBeam

#从偏航、俯仰、翻滚计算旋转矩阵。
def CalPlatformRotationMatrix(globalYaw, globalPitch, globalRoll):
    """
    globalYaw, globalPitch, globalRoll: scalar  float64
    """
    
    cy = math.cos(globalYaw) 
    sy = math.sin(globalYaw) 
    cp = math.cos(globalPitch) 
    sp = math.sin(globalPitch) 
    cr = math.cos(globalRoll) 
    sr = math.sin(globalRoll) 
    
    rotationMatrixPlatform = np.zeros((3,3), dtype='float64')
    rotationMatrixPlatform[0,0] = cy * cp 
    rotationMatrixPlatform[0,1] = cy * sp * sr - sy * cr 
    rotationMatrixPlatform[0,2] = sy * sr + cy * sp * cr 
    rotationMatrixPlatform[1,0] = sy * cp 
    rotationMatrixPlatform[1,1] = cy * cr + sy * sp * sr 
    rotationMatrixPlatform[1,2] = sy * sp * cr - cy * sr 
    rotationMatrixPlatform[2,0] = -sp 
    rotationMatrixPlatform[2,1] = cp * sr 
    rotationMatrixPlatform[2,2] = cp * cr 
        
    return rotationMatrixPlatform

def VectorTransform(matrixA, vectorB):
    """
    """ 
    
    m11 = matrixA[0, 0] * vectorB[0] + matrixA[0, 1] * vectorB[1] + matrixA[0, 2] * vectorB[2]
    m21 = matrixA[1, 0] * vectorB[0] + matrixA[1, 1] * vectorB[1] + matrixA[1, 2] * vectorB[2]
    m31 = matrixA[2, 0] * vectorB[0] + matrixA[2, 1] * vectorB[1] + matrixA[2, 2] * vectorB[2]

   
    return  np.array([m11, m21, m31])
