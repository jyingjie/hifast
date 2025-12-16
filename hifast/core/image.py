

__all__ = ['SpecFile', 'align_chan', 'Imaging', 'CalcSpecPixel', 'args']


import numpy as np
import h5py
from glob import glob
import copy
import os
import sys
import copy

from astropy.coordinates import SkyCoord
from astropy import units as u

from .wcs import gen_header, gen_grid_radec, header_3to2
from .beam import ConvFun
from .radec import _tight_ra
from . import conf
from .corr_vel import freq2vel, vel2freq
from .parallel_fill_grid import optimized_search_around_sky
from ..utils.io import get_nB
from ..utils.io import replace_nB

from multiprocessing import RawArray
from multiprocessing import Process, Queue

try:
    import bottleneck as bn
    NANMEAN = bn.nanmean
    NANSUM = bn.nansum
except ImportError:
    NANMEAN = np.nanmean
    NANSUM = np.nansum


import re
def groups_print_fnames(fnames):
    from itertools import groupby
    nBs = [get_nB(fname) for fname in fnames]
    replace_nB = lambda path,s:re.sub(r'[0-9][0-1]M-', f"-M{s}"[::-1], path[::-1], count=1)[::-1]
    pairs = [(replace_nB(fname, '{Beam}'), get_nB(fname)) for fname in fnames]
    pairs.sort()
    for key, group in groupby(pairs, lambda x:x[0]):
        print(key)
        print('Beams:', ', '.join([str(pp[1]) for pp in group]))


class SpecFile():
    """
    load spec info from hdf5 file
    """
    def __init__(self, fpath, type3='freq'):

        self.fpath = fpath
        self.type3 = type3
        self.nB = get_nB(fpath)
        self.load_misc()

    def __repr__(self):
        _repr = f"File: {self.fpath}"
        return _repr

    def get_frame(self, f):
        """
        try to read rest frame from header
        """
        try:
            self.frame = f['Header'].attrs['frame']
        except:
            self.frame = ''

    def load_misc(self, ):
        with h5py.File(self.fpath, 'r') as f:
            S = f['S']
            self.ra = S['ra'][()]
            self.dec = S['dec'][()]
            self.mjd = S['mjd'][()]
            self.freq = S['freq'][()]
            self.get_frame(f)
        # not here
        # self.ra = _tight_ra(self.ra)

        self.nchan_ori = len(self.freq)
        self.nspec_ori = len(self.ra)

        if self.type3.upper() == 'FREQ':
            self.arr3 = self.freq
        elif self.type3.upper() == 'VOPT':
            self.arr3 = freq2vel(self.freq, vtype='optical')
        elif self.type3.upper() == 'VRAD':
            self.arr3 = freq2vel(self.freq, vtype='radio')

    def select_spec(self, ra_range, dec_range, r_add):
        """
        select spectra in ra and dec range; generate self.is_use_t; also change ra, dec
        """
        if ra_range is not None:
            ra_tmp = self.ra
            is_use_t = (ra_tmp >= (ra_range[0] - 2*r_add)) & (ra_tmp <= (ra_range[1] + 2*r_add))
            is_use_t |= ((ra_tmp >= (ra_range[0]+360 - 2*r_add)) & (ra_tmp <= (ra_range[1]+360 + 2*r_add)))
            is_use_t |= ((ra_tmp >= (ra_range[0]-360 - 2*r_add)) & (ra_tmp <= (ra_range[1]-360 + 2*r_add)))
        else:
            is_use_t = np.full(len(self.ra), True)
        if dec_range is not None:
            dec_ = self.dec
            is_use_t &= (dec_ >= (dec_range[0] - 2*r_add)) & (dec_ <= (dec_range[1] + 2*r_add))
        if np.all(is_use_t):
            self.is_use_t = None
        else:
            self.is_use_t = is_use_t
            self.ra = self.ra[is_use_t]
            self.dec = self.dec[is_use_t]


def align_chan(arr3s, inds, range3=None):
    """
    vlrs, Ta: list; arr3s[i] in descending order
    """
    ranges = np.array([np.min([i[0] for i in arr3s]), np.max([i[-1] for i in arr3s])])
    delta = (arr3s[0][0] - arr3s[0][1])/2.
    vmax, vmin = ranges.max() + delta, ranges.min()-delta
    if range3 is not None:
        vmax = min(vmax, range3[1])
        vmin = max(vmin, range3[0])
    _arr3s = []
    _inds = []

    for i, j in zip(arr3s, inds):
        is_use = (i >= vmin) & (i <= vmax)
        _arr3s += [i[:][is_use],]
        _inds += [j[:][is_use],]

    arr3s_len = np.array([len(i) for i in _arr3s])
    len_use = arr3s_len.min()

    arr3s = [i[:len_use] for i in _arr3s]
    inds = [i[:len_use] for i in _inds]

    return arr3s, inds


class Imaging():
    """
    """
    def __init__(self, args, inplace_args=False):

        args = args if inplace_args else copy.deepcopy(args)
        if len(args.bwidth) == 1:
            args.bwidth = args.bwidth*2
        self.args = args

        # Set temporary directory for multiprocessing
        self._setup_temp_dir()

        self.check_outname()

    def _setup_temp_dir(self):
        """
        Configure temporary directory for multiprocessing.
        Only enabled when --share_mem is True.
        """
        args = self.args

        # Only set temp_dir when share_mem is enabled
        if not getattr(args, 'share_mem', False):
            return

        # Only apply on Linux
        if sys.platform != 'linux':
            print("Warning: --share_mem is only supported on Linux, ignoring temp_dir setting")
            return

        temp_dir = getattr(args, 'temp_dir', '/dev/shm')

        # Check if temp_dir exists and is writable
        if not os.path.exists(temp_dir):
            print(f"Error: Temporary directory '{temp_dir}' does not exist.")
            print("Please create the directory or specify a different one with --temp_dir")
            sys.exit(1)

        if not os.access(temp_dir, os.W_OK):
            print(f"Error: Temporary directory '{temp_dir}' is not writable.")
            print("Please check permissions or specify a different directory with --temp_dir")
            sys.exit(1)

        import tempfile
        tempfile.tempdir = temp_dir
        print(f"Using temporary directory for multiprocessing: {temp_dir}")
        if temp_dir == '/dev/shm':
            print("  (in-memory, fast but limited by RAM)")
        else:
            print("  (disk-based, slower but more space)")

    def check_outname(self):
        args = self.args
        if os.path.exists(args.outname):
            if args.force:
                print(f"will overwrite the existing output file {args.outname}")
            else:
                print(f"File exists {args.outname}")
                print("exit... Use ' -f ' to overwrite it or change outname.")
                sys.exit(0)


    def gen_fpaths(self, ):
        """fpaths from fpattern"""
        args = self.args

        fpaths = []
        for fpattern in args.fpatterns:
            fpaths += glob(fpattern)
        fpaths = list(set(fpaths))
        fpaths.sort()
        if len(fpaths) == 0:
            raise(ValueError('no file input'))
        try:
            groups_print_fnames(fpaths)
        except:
            print(*fpaths, sep='\n')
        self.fpaths = fpaths

    def gen_specfiles(self, ):
        """
        read basic info from fpath, and align chan axis
        """
        args = self.args

        specfiles = [SpecFile(fpath, args.type3) for fpath in self.fpaths]
        [s.select_spec(args.ra_range, args.dec_range, (args.r_cut+max(args.bwidth))/3600)
         for s in specfiles]
        specfiles = np.asarray(specfiles, dtype='object')

        ## align chan axis
        is_desc = specfiles[0].arr3[1] - specfiles[0].arr3[0] < 0
        if is_desc:
            arr3s = [sf.arr3 for sf in specfiles]
            inds = [np.arange(sf.nchan_ori) for sf in specfiles]
        else:
            arr3s = [sf.arr3[::-1] for sf in specfiles]
            inds = [np.arange(sf.nchan_ori)[::-1] for sf in specfiles]
        arr3s, inds = align_chan(arr3s, inds, args.range3)
        for i, (ind_use, arr3) in enumerate(zip(inds, arr3s)):
            if not is_desc:
                ind_use = ind_use[::-1]
                arr3 = arr3[::-1]
            if len(ind_use) != specfiles[i].nchan_ori:
                specfiles[i].ind_use_chan = ind_use
            else:
                specfiles[i].ind_use_chan = None
            specfiles[i].arr3 = arr3
        ## align end
        self.specfiles = specfiles

    def check_spec_key(self,):
        args = self.args
        with h5py.File(self.specfiles[0].fpath, 'r') as f:
            if args.key is None:
                if 'flux' in f['S'].keys():
                    args.key = 'flux'
                elif 'Ta' in f['S'].keys():
                    args.key = 'Ta'
                elif 'Power' in f['S'].keys():
                    args.key = 'Power'
                else:
                    raise(ValueError('no flux or temperature information'))
            print(f"Use field `{args.key}`")

    def check_axis3(self,):
        """"""
        args = self.args
        #check channel align
        arr3_stack = np.vstack([sf.arr3 for sf in self.specfiles])
        std_ = np.std(arr3_stack, axis=0, dtype=np.float64) # single precision can be inaccurate
        print(f'{args.type3} dispersion (std) at the same channel is between {np.nanmin(std_)} and {np.nanmax(std_)}.')
        if np.nanmax(std_)>0.1:
            raise(ValueError('vel dispersion (std) at the same channel is two large'))

    def get_history(self,):
        args = self.args
        #record history
        histories = []
        try:
            from .._version import get_versions
            histories = ['version: ' + get_versions()['version']]
            histories += [f"{k}: {v}" if k !='fpatterns' else "args" for k,v in args.__dict__.items()]
        except:
            pass
        histories += self.fpaths
        return histories

    def change_header(self,):
        """
        change header if args.wcs_from
        """
        args = self.args
        if args.wcs_from is not None:
            print(f'use the wcs sky coordinates parameters from {args.wcs_from}')
            from astropy.io import fits
            fa = fits.open(args.wcs_from)
            header2 = fa[0].header
            for key in self.header.keys():
                if key[-1:] in ['1', '2'] or key == 'LONPOLE' or key == 'LATPOLE':
                    print(f'replacing {key}')
                    try:
                        self.header[key] = header2[key]
                    except:
                        print(f'replace {key} fail')

    def gen_header_grid(self, ):
        """

        """
        args = self.args

        self.ra_stack = ra_stack = _tight_ra(np.hstack([sf.ra for sf in self.specfiles]))
        self.dec_stack = dec_stack = np.hstack([sf.dec for sf in self.specfiles])
        if args.ra_range is None:
            args.ra_range = ra_stack.min(), ra_stack.max()
        if args.dec_range is None:
            args.dec_range = dec_stack.min(), dec_stack.max()

        frame = self.specfiles[0].frame

        #
        print('generating WCS header')
        print('ra range:', args.ra_range,' deg')
        print('dec range:', args.dec_range,' deg')
        print('bwidth:', *args.bwidth, ' arcsec')
        self.header= gen_header(args.ra_range, args.dec_range,
                                args.bwidth[0]/3600, args.bwidth[1]/3600,
                                self.specfiles[0].arr3, args.proj, args.type3, frame,
                                histories=self.get_history(),
                                beam_fwhw = args.beam_fwhw,
                                ra_center=args.wcs_ra_center, dec_center=args.wcs_dec_center)
        if args.key == 'flux':
            self.header["BUNIT"] = 'Jy/beam'
        elif args.key == 'Ta':
            self.header["BUNIT"] = 'K'
        # args.wcs_from
        self.change_header()
        # gen grid points
        self.ra_grid, self.dec_grid = gen_grid_radec(self.header)



    def fill_grid(self, ):
        args = self.args

        cata = SkyCoord(self.ra_stack, self.dec_stack, unit=(u.degree, u.degree))
        grid = SkyCoord(self.ra_grid, self.dec_grid, unit=(u.degree, u.degree))

        use_optimized = getattr(args, 'use_optimized_search', True)

        if use_optimized:
            print(f'searching spectra in r_cut for each grid (optimized, n_workers={args.nproc})...')
            # optimized_search_around_sky(coords1, coords2) returns (idx1, idx2, d2d, d3d)
            # cata.search_around_sky(grid) calls search_around_sky(grid, cata) internally
            # So we need to call optimized_search_around_sky(grid, cata) to match
            # compute_3d=False to skip unused 3D distance calculation
            self.ind_g, self.ind_cata, self.d2d, _ = optimized_search_around_sky(
                grid.ravel(), cata.ravel(), args.r_cut*u.arcsec,
                n_workers=args.nproc, compute_3d=False, verbose=False
            )
        else:
            print('searching spectra in r_cut for each grid (original astropy)...')
            # cata.search_around_sky(grid) returns (idx_grid, idx_cata, d2d, d3d)
            self.ind_g, self.ind_cata, self.d2d, _ = cata.ravel().search_around_sky(
                grid.ravel(), args.r_cut*u.arcsec
            )

    def gen_scale_beams(self,):
        args = self.args
        import json
        if args.scale_beams_file is not None:
            beams_c = json.load(open(args.scale_beams_file))
#             # two dim Volume = 2*pi*A*sigma^2
#             self.scale_beams = (beams_c['beam_fwhw_synt'] / np.array(beams_c['beam_fwhw_19']))**2
            self.scale_beams = np.array(beams_c['scale_beams_19'])
            assert len(self.scale_beams) == 19
            args.beam_fwhw = beams_c['beam_fwhw_synt'] # arcmin
        else:
            self.scale_beams = None

    def set_r_cut(self,):
        args = self.args
        if args.r_cut is None:
            cf = ConvFun(args.method, args.beam_fwhw,
                         gaussian_fwhw=args.gaussian_fwhw,
                         bsize=getattr(args, 'bsize', None), gsize=getattr(args, 'gsize', None))
            args.r_cut = getattr(cf, 'dis_cut_suggest')*60 # to arcsec
        print(f'r_cut: {args.r_cut} arcsec')


    def gen_conv_weis(self,):
        args = self.args
        cf = ConvFun(args.method, args.beam_fwhw,
                     gaussian_fwhw=args.gaussian_fwhw,
                     bsize=getattr(args, 'bsize', None), gsize=getattr(args, 'gsize', None))
        self.conv_weis = cf(self.d2d.arcmin)
        self.conv_obj = cf


    def save_cubes(self):
        """
        Save data cubes to FITS files sequentially to minimize memory usage.

        This method creates cube views on-demand and immediately deletes the underlying
        pixel arrays after each save, avoiding memory spikes from keeping both in memory.

        Memory behavior notes:
        - Data is computed in float64 for accuracy during accumulation/normalization
        - Saved as float32 to reduce file size (standard for FITS cubes)
        - astype('float32') creates a temporary copy in RAM, regardless of whether
          the source data is in RAM or disk-based shared memory (RawArray with temp_dir)
        - Sequential save-and-delete approach minimizes peak RAM usage by only having
          one float32 conversion in memory at a time
        - When share_mem=True: source data may be on disk (temp_dir), but astype()
          still allocates the float32 copy in RAM. This is unavoidable but brief.
        """
        args = self.args
        from astropy.io import fits
        import gc

        # Calculate shape for reshaping pixel arrays into cubes
        _shape = self.ra_grid.shape + self.pixel_data.shape[-1:]

        # Apply beam correction to pixel_data before reshaping
        # self.data_cube is a view, so modifying it modifies self.pixel_data
        self.data_cube = self.pixel_data.reshape(_shape).transpose(2,0,1)
        # Apply beam correction (modifies self.header and self.data_cube)
        # Must be done before creating _header (header copy) and saving other cubes
        # Note: Only pixel_data needs correction; pixel_finite_nums_chan and pixel_weis_sum_chan do not.
        self.apply_beam_correction()

        _header = copy.deepcopy(self.header)
        _header['BUNIT'] = ''

        # Save nums_chan_cube first and delete pixel_finite_nums_chan immediately
        outname_ = '.'.join(args.outname.split('.')[:-1]) + '-count.fits'
        nums_chan_cube = self.pixel_finite_nums_chan.reshape(_shape).transpose(2,0,1)
        hdu = fits.PrimaryHDU(nums_chan_cube, header=_header)
        print(f'Saving to {outname_}.')
        hdu.writeto(outname_, overwrite=True)
        print(f'Saved count cube.')
        del nums_chan_cube
        del self.pixel_finite_nums_chan  # Free underlying memory
        gc.collect()  # Force garbage collection to ensure memory is freed


        # Save main data cube and delete pixel_data immediately
        hdu = fits.PrimaryHDU(self.data_cube.astype('float32'), header=self.header)
        print(f'Saving to {args.outname}.')
        hdu.writeto(args.outname, overwrite=args.force)
        print(f'Saved data cube.')
        del self.data_cube
        del self.pixel_data  # Free underlying memory
        gc.collect()  # Force garbage collection to ensure memory is freed

        # Save weights cube and delete pixel_weis_sum_chan immediately
        outname_ = '.'.join(args.outname.split('.')[:-1]) + '-weights.fits'
        weis_chan_cube = self.pixel_weis_sum_chan.reshape(_shape).transpose(2,0,1)
        hdu = fits.PrimaryHDU(weis_chan_cube.astype('float32'), header=_header)
        print(f'Saving to {outname_}.')
        hdu.writeto(outname_, overwrite=True)
        print(f'Saved weights cube.')
        del weis_chan_cube
        del self.pixel_weis_sum_chan  # Free underlying memory
        del self.pixel_specs_nums  # Clean up remaining arrays
        gc.collect()  # Force garbage collection to ensure memory is freed


    def init_out(self, DataType='float64'):
        args = self.args
        npixel = np.prod(self.ra_grid.shape)
        # init
        shape = (npixel, len(self.specfiles[0].arr3))
        if args.share_mem:
            # ctype: https://docs.python.org/3/library/array.html#module-array
            if DataType == 'float64':
                shm_dtype = 'd'
            elif DataType == 'float32':
                shm_dtype = 'f'
            else:
                raise(ValueError('DataType not support'))
            self.DataType = DataType

            shm = RawArray(shm_dtype, int(np.prod(shape)))
            self.pixel_data = np.frombuffer(shm, dtype=DataType).reshape(shape)
            self.pixel_data[:] = 0
            #
            shm = RawArray(shm_dtype, int(np.prod(shape)))
            self.pixel_weis_sum_chan = np.frombuffer(shm, dtype=DataType).reshape(shape)
            self.pixel_weis_sum_chan[:] = 0
            # deal with nan value
            shm = RawArray('i', int(np.prod(shape)))
            self.pixel_finite_nums_chan = np.frombuffer(shm, dtype='int32').reshape(shape)
            self.pixel_finite_nums_chan[:] = 0

        else:
            self.pixel_data = np.zeros(shape, dtype=DataType)
            #
            self.pixel_weis_sum_chan = np.zeros(shape, dtype=DataType)
            #
            self.pixel_finite_nums_chan = np.zeros(shape, dtype='int32')
        self.pixel_specs_nums = np.zeros(npixel, dtype='int32')

    def apply_beam_correction(self, ):
        args = self.args
        new_beam_fwhw = self.conv_obj.new_beam_fwhw
        if args.apply_beam_correction:
            print(f"Modifying BMAJ and BMIN in header from {self.header['BMAJ']*60, self.header['BMIN']*60} arcmin to {new_beam_fwhw, new_beam_fwhw} arcmin, ",
                 "to reflect the spatial resolution change due to the convolution.")
            self.header['BMAJ'] = new_beam_fwhw/60 # deg
            self.header['BMIN'] = new_beam_fwhw/60 # deg

            # multiple the ratio factor only when BUNIT is Jy/beam
            print("If the BUNIT is Jy/beam, a ratio factor is needed to correct the data.")
            if self.header["BUNIT"].lower() == 'jy/beam':
                print(f"The BUNIT is Jy/beam, multiplying the beam ratio factor {self.conv_obj.beam_factor}")
                self.data_cube *= self.conv_obj.beam_factor
            else:
                print("The BUNIT is not Jy/beam, no need to multiply the beam ratio factor")
                pass
        else:
            self.header['BMAJ_new'] = new_beam_fwhw/60
            self.header['BMIN_new'] = new_beam_fwhw/60


    def __call__(self, ):
        """
        """
        # prepare
        args = self.args
        step = args.step

        self.gen_fpaths()
        self.gen_scale_beams() # also change args.beam_fwhw
        self.set_r_cut()
        self.gen_specfiles()
        self.check_spec_key()
        self.check_axis3()
        self.gen_header_grid()
        self.fill_grid()
        self.gen_conv_weis()

        ## todo: consider other weights
        weis = self.conv_weis

        self.init_out()
        # loop step
        nspecs = np.array([len(sf.ra) for sf in self.specfiles])
        starts = np.arange(0, len(self.specfiles), step)
        for i, start in enumerate(starts):
            print(f"Part [{i+1}/{len(starts)}]")
            sfs = self.specfiles[start:start+step]
            # index in current specs
            ind_cata_the = self.ind_cata - np.sum(nspecs[:start])
            # select current ind_g, ind_cata, wei,
            is_ = (ind_cata_the < np.sum(nspecs[start:start+step])) & (ind_cata_the >=0)
            ind_cata_the_use = ind_cata_the[is_]
            # if empty
            if len(ind_cata_the_use) == 0:
                # notice: make sure the skip not affect the following steps
                continue
            ind_g_use = self.ind_g[is_]
            weis_use = weis[is_]
            # check
            assert (self.dec_stack[self.ind_cata[is_]] == np.hstack([sf.dec for sf in sfs])[ind_cata_the_use]).all()

            # Note: CalcSpecPixel's share_mem is temporarily disabled (set to False)
            # to avoid potential issues with shared memory in multiprocessing.
            # The Imaging class's share_mem (in init_out) is still functional.
            csp = CalcSpecPixel(sfs, key=args.key, polar=args.polar, scale_beams=self.scale_beams, share_mem=False)
            res = csp(ind_g_use, ind_cata_the_use, weis_use, n_worker=args.nproc)
            ii = res[0]
            self.pixel_data[ii] += res[1]
            self.pixel_weis_sum_chan[ii] += res[2]
            self.pixel_specs_nums[ii] += res[3]
            self.pixel_finite_nums_chan[ii] += res[4]
        # normalization
        old_set = np.seterr()
        np.seterr(invalid='ignore')
        self.pixel_data /= self.pixel_weis_sum_chan
        self.pixel_data[self.pixel_weis_sum_chan==0] = np.nan
        if args.frac_finite_min > 0:
            self.pixel_data[self.pixel_finite_nums_chan < (self.pixel_specs_nums*args.frac_finite_min)[:, None]] = np.nan
        np.seterr(**old_set)

        # Save each cube directly and delete underlying pixel arrays immediately
        # to minimize memory usage (especially important with share_mem)
        self.save_cubes()


def _get_start_stop(arr,arr_in):
    """
    arr_in: array
       sorted
    """
    if not np.all(arr_in[:-1] <= arr_in[1:]):
        raise(ValueError('arr_in should be sorted'))
    inds_left= np.searchsorted(arr_in, arr, side='left',)
    inds_right=np.searchsorted(arr_in, arr, side='right',)
    return inds_left, inds_right


class CalcSpecPixel():
    def __init__(self, sfs, *, key='flux', polar='M', scale_beams=None, share_mem=False):
        """
        key:
        polar: 'M','XX','YY'
        """
        # init shared specs
        self.polar = polar
        self.scale_beams = scale_beams
        self.share_mem = share_mem
        self.ns = ns = [len(sf.ra) for sf in sfs]

        DataType = 'float32' # use same with that in file, i.e. float32, will convert to float64 in self.calc
        specs_shape = (np.sum(ns), len(sfs[0].arr3))
        if self.share_mem:
            # ctype: https://docs.python.org/3/library/array.html#module-array
            if DataType == 'float64':
                shm_dtype = 'd'
            elif DataType == 'float32':
                shm_dtype = 'f'
            else:
                raise(ValueError('DataType not support'))
            self.shm_specs = RawArray(shm_dtype, int(np.prod(specs_shape)))
            self.specs = np.frombuffer(self.shm_specs, dtype=DataType).reshape(specs_shape)
        else:
            self.shm_specs = self.specs = np.zeros(specs_shape, dtype=DataType)

        self.DataType = DataType
        self.key = key
        self.load_data(sfs)

    def load_data(self, sfs):

        for sf, s in zip(sfs, np.cumsum([0, *self.ns[:-1]])):
            f = h5py.File(sf.fpath)
            print(f"loading: {sf.fpath}", end='\033[K\r')
            f_spec = f['S'][self.key]
            if self.polar == 'XX':
                sli = 0
            elif self.polar == 'YY':
                if f_spec.shape[0] < 2:
                    raise(ValueError(f'polar method is YY, but spectra have only one polar.'))
                sli = 1
            else:
                sli = slice(None)
            if sf.is_use_t is not None:
                ind_use_t = np.where(sf.is_use_t)[0]
                if sf.ind_use_chan is not None:
                    if (np.diff(sf.ind_use_chan) == 1).all():
                        _specs = f_spec[sli, ind_use_t, sf.ind_use_chan[0]:sf.ind_use_chan[-1]+1]
                    else:
                        raise(ValueError('select chans fail'))
                else:
                    _specs = f_spec[sli, ind_use_t]
            else:
                if sf.ind_use_chan is not None:
                    if (np.diff(sf.ind_use_chan) == 1).all():
                        _specs = f_spec[sli, :, sf.ind_use_chan[0]:sf.ind_use_chan[-1]+1]
                    else:
                        raise(ValueError('select chans fail'))
                else:
                    _specs = f_spec[sli]
            if self.polar == 'M':
                _specs = np.mean(_specs, axis=0) # only two float32 to average
            elif self.polar in ['XX', 'YY']:
                ...
            elif self.polar == 'B':
                ...
            else:
                raise(ValueError(f'polar method not support {self.polar}'))
            if self.scale_beams is not None:
                _specs *= self.scale_beams[sf.nB-1]
                print(self.scale_beams[sf.nB-1])
            self.specs[s:s+len(_specs)] = _specs.astype(self.DataType)
            f.close()
        print('')

    @staticmethod
    def calc(q, s, e,
             inds, inds_dtype, inds_shape,
             specs, specs_dtype, specs_shape,
             weis, weis_dtype, weis_shape):

        if not isinstance(inds, np.ndarray):
            inds = np.frombuffer(inds, dtype=inds_dtype).reshape(inds_shape)
        if not isinstance(specs, np.ndarray):
            specs = np.frombuffer(specs, dtype=specs_dtype).reshape(specs_shape)
        if not isinstance(weis, np.ndarray):
            weis = np.frombuffer(weis, dtype=weis_dtype).reshape(weis_shape)

        _SUM = NANSUM
        val_c = []
        num = []
        num_finite_chan = []
        wei_sum_chan = []
        for s_, e_ in zip(s, e):
            val = specs[inds[s_:e_]]
            wei = weis[s_:e_] # wei is align with inds
            wei = wei.reshape((-1,)+(1,)*(val.ndim-1))
            # deal with nan value
            is_finite = np.isfinite(val)
            val_c += [_SUM(val.astype('float64')*wei, axis=0),]
            num += [len(wei),]
            num_finite_chan += [np.sum(is_finite, axis=0),]
            wei_sum_chan += [_SUM(wei*is_finite, axis=0),]
        res = np.vstack(val_c), np.vstack(wei_sum_chan), np.hstack(num), np.vstack(num_finite_chan)

#         import time
#         print('ssssss')
#         time.sleep(15)
#         res = [None, None, None]
        if q is None:
            # not multiprocessing
            return res
        else:
            q.put(res)

    def __call__(self, ind_g_use, ind_cata_the_use, weis_use, n_worker=3):

        if self.share_mem:
            # ctype: https://docs.python.org/3/library/array.html#module-array
            shm1 = RawArray('q', int(np.prod(ind_cata_the_use.shape)))
            shm1_arr = np.frombuffer(shm1, dtype='int64').reshape(ind_cata_the_use.shape)
            shm1_arr[:] = ind_cata_the_use
            shm2 = RawArray('d', int(np.prod(weis_use.shape)))
            shm2_arr = np.frombuffer(shm2, dtype='float64').reshape(weis_use.shape)
            shm2_arr[:] = weis_use
        else:
            shm1 = shm1_arr = ind_cata_the_use
            shm2 = shm2_arr = weis_use

        print('spread these spectra into grids')
        ind_g_use_uni = np.unique(ind_g_use)
        start, stop = _get_start_stop(ind_g_use_uni, ind_g_use)

#         print('dist')
#         import time
#         time.sleep(10)
        n_worker = min(len(start), n_worker)
        if n_worker > 1:
            print(f"assign gird to {n_worker} processes")
            ps = []
            qs = []
            start_list = np.array_split(start, n_worker)
            stop_list = np.array_split(stop, n_worker)
            for start_, stop_ in zip(start_list, stop_list):
                q = Queue()
                p = Process(target=self.calc, args=(q, start_, stop_,
                                                   shm1, shm1_arr.dtype, shm1_arr.shape,
                                                   self.shm_specs, self.specs.dtype, self.specs.shape,
                                                   shm2, shm2_arr.dtype, shm2_arr.shape),)
                ps += [p]
                qs += [q]
            for p in ps:
                p.start()
            res = [q.get() for q in qs]
            for p in ps:
                p.join()  # need after q.get()
            print('gathering results')
            pixel_data = np.vstack([r[0] for r in res])
            pixel_weis_sum_chan = np.vstack([r[1] for r in res])
            pixel_specs_nums = np.hstack([r[2] for r in res])
            pixel_finite_nums_chan = np.vstack([r[3] for r in res])
        else:
            pixel_data, pixel_weis_sum_chan, pixel_specs_nums, pixel_finite_nums_chan = self.calc(None, start, stop,
                                                   shm1, shm1_arr.dtype, shm1_arr.shape,
                                                   self.shm_specs, self.specs.dtype, self.specs.shape,
                                                   shm2, shm2_arr.dtype, shm2_arr.shape)

        return ind_g_use_uni, pixel_data, pixel_weis_sum_chan, pixel_specs_nums, pixel_finite_nums_chan


from argparse import Namespace

## exampler parameters
args = Namespace()

args.ra_range = [20.254, 27.25]
args.dec_range = None

args.wcs_ra_center = None
args.wcs_dec_center = None

args.range3 = [-2000, 1500]
args.type3 = 'vrad'
args.bwidth = [60,]
args.proj = 'SIN'
args.r_cut = 90 # arcsec
##
args.method = 'bessel_gaussian'
args.beam_fwhw = 2.9 # arcmin; beam Full width at half maximum
args.gaussian_fwhw = None
args.bsize = None # for Bessel in Bessel*gauss or sin in sin*gauss
args.gsize = None # for gauss in Bessel*gauss or sin*gauss

args.fpatterns = ['/home/jyj/jingyj/pipeline_test/M31_v4_highz/S2/med/*-bld-fc-ds.hdf5', ]

args.outname = 'test.fits'
args.force = True

args.nproc = 5
args.step = 19
args.share_mem = False

args.key = None
args.wcs_from = '/home/jyj/jingyj/M33/M33_local_areciob.fits'

args.scale_beams_file = None

args.polar = 'M'

args.frac_finite_min = 0.01
