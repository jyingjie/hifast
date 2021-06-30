#!/usr/bin/env python
# coding: utf-8
import h5py
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u

class Find_point(object):
    def __init__(self, p_str, r, p_u=(u.hourangle,u.deg), r_u=u.deg):
        self.point = SkyCoord(p_str, unit=p_u)
        self.r = r*u.deg
    @staticmethod
    def get_radec(fname):
        f = h5py.File(fname,'r')
        ra_list = []
        beam_list = []
        for key in f.keys():
            if 'ra' == key[:2]:
                try:
                    beam_list += [key[2:]]
                except:
                    beam_list += ['']
                ra_list += [f[key][()]]

        dec_list = [f[f'dec{beam}'][()] for beam in beam_list]
        radec = {}
        for ra, dec, beam in zip(ra_list, dec_list, beam_list):
            radec[beam] = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
        f.close()
        return radec
    def find_in(self, fname):
        radec = self.get_radec(fname)
        for beam in radec.keys():
            inds = np.where(radec[beam].separation(self.point) < self.r)[0]
            if len(inds)>0:
                print(f'Beam {beam:>2} in', fname)
                print(inds)

if __name__ == '__main__':
    import os
    import argparse
    from glob import glob
    
    parser = argparse.ArgumentParser(allow_abbrev=False)
    parser.add_argument('fnames',nargs='+',
                        help='fnames')
    parser.add_argument('-p', required=True,
                       help="point ra dec str, example: \"0:48:26.3991 +42:34:08.808\"")
    parser.add_argument('-r', type=float, default=1, 
                       help='radius, unit: arcmin, default 1')
    
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