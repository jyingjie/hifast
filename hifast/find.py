#!/usr/bin/env python
# coding: utf-8
"""
hifast.find
===========

Search for spectra near a target coordinate.

This module searches for beams within a specified radius of a target RA/DEC in a list of HDF5 files.
"""
import h5py
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
import os
import argparse
from glob import glob
import sys

class Find_point(object):
    def __init__(self, p_str, r, p_u=(u.hourangle,u.deg), r_u=u.deg):
        if ':' in p_str:
            # Sexagesimal format (HH:MM:SS DD:MM:SS) -> Hourangle, Deg
            self.point = SkyCoord(p_str, unit=p_u)
        else:
            try:
                # Try auto-detect units (e.g., "30.3deg 40.44deg", "12h30m 45d")
                self.point = SkyCoord(p_str)
            except Exception:
                # Fallback to decimal degrees (e.g., "30.3 40.44")
                self.point = SkyCoord(p_str, unit=(u.deg, u.deg))
        
        self.r = r*u.deg
    @staticmethod
    def get_radec(fname):
        f = h5py.File(fname,'r')
        S = f["S"] if "S" in f.keys() else f
        ra_list = []
        beam_list = []
        for key in S.keys():
            if 'ra' == key[:2]:
                try:
                    beam_list += [key[2:]]
                except:
                    beam_list += ['']
                ra_list += [S[key][()]]

        dec_list = [S[f'dec{beam}'][()] for beam in beam_list]
        radec = {}
        for ra, dec, beam in zip(ra_list, dec_list, beam_list):
            radec[beam] = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
        f.close()
        return radec
    def find_in(self, fname):
        radec = self.get_radec(fname)
        for beam in radec.keys():
            # Handle both scalar and array coordinates
            sep = radec[beam].separation(self.point)
            if sep.isscalar:
                sep = np.atleast_1d(sep)
            
            inds = np.where(sep < self.r)[0]
            if len(inds)>0:
                print(f'Beam {beam:>2} in', fname)
                print(inds)

def create_parser():
    parser = argparse.ArgumentParser(prog=f"python -m hifast.{os.path.basename(sys.argv[0])[:-3]}",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     allow_abbrev=False,
                                     description='Search for spectra within a radius of a target coordinate.')
    
    parser.add_argument('fnames', nargs='+', metavar='FILE',
                        help='Input HDF5 file paths (supports glob patterns like "data/*.hdf5").')
    
    parser.add_argument('-p', required=True, metavar='RA_DEC',
                        help='Target coordinate string. Supports "HH:MM:SS DD:MM:SS" (e.g., "0:48:26.39 +42:34:08"), decimal degrees (e.g., "30.3 40.44"), or explicit units (e.g., "30.3deg 40.44deg").')
    
    parser.add_argument('-r', type=float, default=1, metavar='ARCMIN',
                        help='Search radius [arcmin].')
    return parser

parser = create_parser()

if __name__ == '__main__':
    args = parser.parse_args()
    fnames = args.fnames
    p = args.p
    r = args.r
    files=[]
    for thefname in fnames:
        files+= glob(thefname, recursive=True)
    files= list(set(files))
    if len(files)==0:
        import sys
        print('cannot find file')
        sys.exit()
    files.sort()
    #print(files)
    F = Find_point(p, r/60)
    [F.find_in(fname) for fname in files]
