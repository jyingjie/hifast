#from __future__ import print_function
import numpy as np
import sys
from PyAstronomy import pyasl

###########################
#### This Python function is used to calculate Doppler veolcity at FAST site
#### Inputs: 
####     ra, dec: The coordinate in deg. List or array.
####     jd:   Julian date. Same lenght as ra and dec
####
#### Calling sequence:
####     from ugdopplerfast import ugdopplerfast
####     vlsrcor = ugdopplerfast(ra,dec,jd)  
####
#### Outputs:
####     vlsrcor: array with LSR correction values
####
#### History:
####   Written by Ningyu Tang, July 4, 2018


# Coordinates of FAST Site
dawodangelong = 106.8566667    # East positive, deg
dawodangnlat  = 25.65294444    # North positive, deg
altitude      = 1110.0288      # Altitude, m 

## Coordinates of GBT Site
## Testing with GBT tool, http://www.gb.nrao.edu/cgi-bin/radvelcalc.py
## Consistent in 1 m/s
#dawodangelong =  280.1603    # East positive, deg
#dawodangnlat  =  38.4331    # North positive, deg
#altitude      =  807.43      # Altitude, m 
## Coordinates 
#rasrc  = 45.0
#decsrc = 20.

## (Mid-)Time of observation
#jd = 2457963.201

def ugdopplerfast(rasrc,decsrc,jd,frame="LSRK"): 
    """
    frame: LSRK or HELIOCEN
    """
    rasrc  = np.array(rasrc)
    decsrc = np.array(decsrc) 
    jd     = np.array(jd) 

    if rasrc.shape[0] == decsrc.shape[0] == jd.shape[0]:
        nin =  rasrc.shape[0]
        porbh   = np.zeros(nin, dtype= np.double)
        porbhjd = np.zeros(nin, dtype= np.double)

        for nr in range(nin):
            
           ########################Orbital Section#######################
           
           # Heli-centric and Barycentric velocity 
           # NONONONo correction of  Earth rotation 
           #vh1, vb1 = pyasl.baryCorr(jd,rasrc,decsrc,deq=2000.0)
           
           # Heli-centric  velocity and jd for the the middle of exposure
           # With correction of Earth rotation!!!
           # Usage as follows:
           # corr, hjd = pyasl.helcorr(longitude, latitude, altitude, 
           #             ra2000, dec2000, jd, debug=False)
           
            vh, hjd = pyasl.helcorr(dawodangelong,dawodangnlat,altitude,\
                                   rasrc[nr],decsrc[nr],jd[nr])
            porbh[nr]   = vh
            porbhjd[nr] = hjd       
        if frame == "HELIOCEN":
            return porbh
           #print("Heliocentric velocity [km/s]: ", vh)
        elif frame=='LSRK':
            ########################LSR Section#######################

            rasrc_rad  = np.deg2rad(rasrc)
            decsrc_rad = np.deg2rad(decsrc)

            xxsource      = np.zeros((3,nin), dtype= np.double) 
            xxsource[0,:] = np.cos(decsrc_rad) * np.cos(rasrc_rad)
            xxsource[1,:] = np.cos(decsrc_rad) * np.sin(rasrc_rad)
            xxsource[2,:] = np.sin(decsrc_rad)
            pvlsr         = np.zeros(nin,dtype= np.double)


            ##----- LSR SECTION----------
            # THE STANDARD LSR IS DEFINED AS FOLLOWS: THE SUN MOVES AT 20.0 KM/S
            # TOWARD RA=18H, DEC=30.0 DEG IN 1900 EPOCH COORDS
            # Precessed J2000: 18:03:50.24   30:00:16.8 (18.06395556,30.00466667)
            # using PRECESS, this works out to ra=18.063955 dec=30.004661 in 2000 coords.

            # Definition toward 18 h in 1900 Epoch
            ralsr_rad  = np.deg2rad(18.06395556*15.)    
            declsr_rad = np.deg2rad(30.00466667) 


            #FIND THE COMPONENTS OF THE VELOCITY OF THE SUN WRT THE LSR FRAME 
            xxlsr      =  np.zeros((3,nin),dtype=np.double) 
            xxlsr[0,:] =  np.cos(declsr_rad)* np.cos(ralsr_rad)
            xxlsr[1,:] =  np.cos(declsr_rad)* np.sin(ralsr_rad)
            xxlsr[2,:] =  np.sin(declsr_rad)
            vvlsr      =  20.*xxlsr

            #PROJECTED VELOCITY OF THE SUN WRT LSR TO THE SOURCE

            for nr in range(nin):
                pvlsr[nr] = np.sum(vvlsr[:,nr]*xxsource[:,nr])

            vvvlst = -porbh-pvlsr

            return vvvlst
            #return vvvlst,-porbh
        else:
            raise(ValueError('frame not support'))
    else:  
        sys.exit('Error:the data length of RA,DEC,and JD are not equal.')
      

