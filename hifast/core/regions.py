

__all__ = ['read_regions', 'box', 'line', 'replace_region', 'mask_srcs']


import warnings
import math
import numpy as np


def read_regions(fpath):
    with open(fpath, 'r') as f:
        # skip comment
        line = f.readline().strip()
        while line:
            if line.startswith('#'):
                line = f.readline().strip()
            else:
                break

        if not line.startswith('global'):
            warnings.warn(f"format error in {fpath}, the first line not started with '#' should start with 'global'. Only DS9 format is supported, please check {fpath}")
            return None
        # region unit
        line = f.readline().strip()
        if not line:
            warnings.warn(f"no region found in {fpath}")
            return None
        elif line.strip() not in ['physical', 'image']:
            warnings.warn(f"only 'pixel' type regions are supported, no regions are read from {fpath}, please check {fpath}")
            return None
        # read region
        regions = {}
        line = f.readline().strip()
        while line:
            key, remain = line.split('(', 1)
            reg = tuple(map(float, remain.split(')', 1)[0].split(',')))
            if key in regions:
                regions[key] += [reg, ]
            else:
                regions[key] = [reg, ]
            line = f.readline().strip()
        return regions


def box(center_x, center_y, size_x, size_y, PA=None):
    """
    box: (center_x, center_y, size_x, size_y, PA?)
    """
    bottom_left = center_x - size_x/2, center_y - size_y/2
    top_right = center_x + size_x/2, center_y + size_y/2

    bottom_left = [max(math.ceil(num), 0) for num in bottom_left]
    top_right = [max(math.floor(num), 0) for num in top_right]

    return slice(bottom_left[0], top_right[0]+1), slice(bottom_left[1], top_right[1]+1)


def line(start_x, start_y, end_x, end_y):
    """
    line(start_x, start_y, end_x, end_y)
    """

    dx_abs = abs(start_x - end_x)
    dy_abs = abs(start_y - end_y)

    if dx_abs > dy_abs:
        s = max(math.ceil(min(start_x, end_x)), 0)
        e = max(math.floor(max(start_x, end_x)), 0)
        return slice(s, e+1), slice(None)
    else:
        s = max(math.ceil(min(start_y, end_y)), 0)
        e = max(math.floor(max(start_y, end_y)), 0)
        return slice(None), slice(s, e+1)


def replace_region(arr, regions, fill_value=np.nan, inplace=True):
    """
    arr: array; ndim >=2, slice the first two dim to replace
    regions:
    fill_value:
    inplace:
    """

    if not inplace:
        arr = np.copy(arr)
    for key in regions.keys():
        if key == 'box':
            for reg in regions[key]:
                arr[box(*reg)[::-1]] = fill_value # x axis is freq
        elif key == 'line':
            for reg in regions[key]:
                arr[line(*reg)[::-1]] = fill_value
        else:
            warnings.warn(f'region type {key} not support')
    return arr


from . import corr_vel as cv
from astropy.coordinates import SkyCoord


def mask_srcs(fpath, is_exclued, ra, dec, freq, mjd=None, inplace=True, rest_frame='LSRK'):
    """
    fpath: file inculde src position. example:
            #ra[deg], dec[deg], R[arcmin], freq_min[MHz], freq_max[MHz]
            349.63, 30.96, 20, 1300, 1302
            349.6510, 30.96226, 10, 1310, 1311
    is_exclued:
    ra, dec: in deg
    mjd:
    inplace: True change is_exclude inplace.
    rest_frame: if not 'none', fix rest frame of freq
    """
    if not inplace:
        is_exclued = np.copy(is_exclued)
    srcs = np.loadtxt(fpath, comments='#', delimiter=',', ndmin=2)
    cata = SkyCoord(ra=ra, dec=dec, unit='deg')
    for src in srcs:
        c = SkyCoord(src[0], src[1], unit='deg')
        ind = np.where(cata.separation(c).arcmin < src[2])[0]
        if len(ind) > 0:
            if rest_frame != 'none':
                freq_C = cv.frame_correct_freq(freq, ra[ind[0]], dec[ind[0]], mjd[ind[0]], frame=rest_frame)
            else:
                freq_C = freq
            ind_freq = np.where((freq_C >= src[3]) & (freq_C <= src[4]))[0]
            # freq is in ascending
            if len(ind_freq) > 0:
                is_exclued[ind, ind_freq[0]:ind_freq[-1]+1] = True
    return is_exclued
