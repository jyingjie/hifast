

__all__ = ['Imaging2', 'args']


import numpy as np
from matplotlib import pyplot as plt
import h5py

from multiprocessing import Process, Queue

from .image import Imaging, _get_start_stop

try:
    import bottleneck as bn
    MEDIAN = bn.median
    NANMEDIAN = bn.nanmedian
    MEAN = bn.nanmean
    NANMEAN = bn.nanmean
    SUM = bn.nansum
    NANSUM = bn.nansum
    STD = bn.nanstd
except ImportError:
    MEDIAN = np.median
    NANMEDIAN = np.nanmedian
    MEAN = np.mean
    NANMEAN = np.nanmean
    SUM = np.sum
    NANSUM = np.nansum
    STD = np.std


class Imaging2(Imaging):

    def fill_grid(self):
        super().fill_grid()
        # sort by self.ind_g then by self.ind_cata
        inds = np.lexsort((self.ind_cata, self.ind_g))
        self.ind_cata = self.ind_cata[inds]
        self.ind_g = self.ind_g[inds]
        self.d2d = self.d2d[inds]
        # d3d if needed

    def set_r_cut(self,):
        args = self.args
        if args.method in ['gaussian', 'bessel_gaussian', 'sinc_gaussian']:
            super().set_r_cut()
        elif args.method in ['wmean', 'wmedian']:
            # use gaussian weights
            method_bak = args.method
            args.method = 'gaussian'
            super().set_r_cut()
            args.method = method_bak
        elif args.method in ['mean', 'median']:
            # no weights needed
            if args.r_cut is None:
                raise(ValueError(f'The method `{args.method}` requires the user to input r_cut'))
            else:
                print(f'r_cut: {args.r_cut} arcsec')

    def _process_sfs(self,):
        """
        add extra info to specfiles for speeding _load_spec_f1
        """
        for sf in self.specfiles:
            if sf.ind_use_chan is not None:
                if (np.diff(sf.ind_use_chan) == 1).all():
                    sf.chan_use_sli = slice(sf.ind_use_chan[0],sf.ind_use_chan[-1]+1)
                else:
                    raise(ValueError('freq shoud be sliced in continuous way'))
            if sf.is_use_t is not None:
                sf.ind_use_t = np.where(sf.is_use_t)[0]

    @staticmethod
    def _load_spec_f1(sf, inds, key, polar, scale_beams=None):
        """
        sf
        inds: indices after limit ra and dec range
        """
        f = h5py.File(sf.fpath)
        f_spec = f['S'][key]
        if polar == 'XX':
            sli = 0
        elif polar == 'YY':
            if f_spec.shape[0] < 2:
                raise(ValueError(f'polar method is YY, but spectra have only one polar.'))
            sli = 1
        else:
            sli = slice(None)
        if sf.is_use_t is not None:
            ind_use_ori = sf.ind_use_t[inds]
        else:
            ind_use_ori = inds

        if sf.ind_use_chan is not None:
            _specs = f_spec[sli, ind_use_ori, sf.chan_use_sli]
        else:
            _specs = f_spec[sli, ind_use_ori]

        if polar == 'M':
            _specs = np.mean(_specs, axis=0) # only two float32 to average
        elif polar in ['XX', 'YY']:
            ...
        elif polar == 'B':
            ...
        else:
            raise(ValueError(f'polar method not support {polar}'))
        if scale_beams is not None:
            _specs *= scale_beams[sf.nB-1]
            print(self.scale_beams[sf.nB-1])
        return _specs



    def _load_spec(self, s, e):
        """
        load spec for one grid from `self.specfiles`
        """
        args = self.args
        ifile_ = self.ifile[s:e]
        ind_ifile_ = self.ind_ifile[s:e]
        # ifile_ and ind_ifile_ is sorted in ifile_ order and then in ind_ifile_ order after
        # self.ind_g and self.ind_cata (and ...) sorted by self.ind_g then by self.ind_cata
        ifile_uni, ifile_num = np.unique(ifile_, return_counts=True)
        ind_ifile_list = np.split(ind_ifile_, np.cumsum(ifile_num)[:-1])
        specs = np.vstack([self._load_spec_f1(self.specfiles[i], ii, args.key, args.polar)
                          for i, ii in zip(ifile_uni,ind_ifile_list)])

        return specs

    @staticmethod
    def calc(q, obj, start_, stop_, ind_g_uni_, verbose=False):
        args = obj.args

        _MEAN = NANMEAN
        _MEDIAN = NANMEDIAN
        _SUM = NANSUM

        old_set = np.seterr()
        np.seterr(invalid='ignore')

        iter_ = zip(start_, stop_, ind_g_uni_)
        if verbose:
            from tqdm import tqdm
            iter_ = tqdm(iter_, total=len(start_), desc='CPU 0: ', mininterval=2)

        for s, e, ii in iter_:
            specs = obj._load_spec(s, e)
            # deal with nan value
            is_finite = np.isfinite(specs)
            num = len(specs)
            num_finite_chan = np.sum(is_finite, axis=0)

            if args.method in ['gaussian', 'bessel_gaussian', 'sinc_gaussian']:
                # convolve
                weis = obj.conv_weis[s:e]
                weis = weis.reshape((-1,)+(1,)*(specs.ndim-1))
                weis_sum_chan = _SUM(weis*is_finite, axis=0)
                spec_g = _SUM(specs.astype('float64')*weis, axis=0) / weis_sum_chan
                # `weis_sum_chan` only those methods
                obj.pixel_weis_sum_chan[ii] = weis_sum_chan
            elif args.method in ['wmean', 'wmedian']:
                weis = obj.conv_weis[s:e]
                weis = weis.reshape((-1,)+(1,)*(specs.ndim-1))
                if args.method == 'wmean':
                    STAT = _MEAN
                else:
                    STAT = _MEDIAN
                spec_g = STAT(specs.astype('float64')*weis, axis=0) / STAT(weis*is_finite, axis=0)
            elif args.method == 'mean':
                spec_g = _MEAN(specs.astype('float64'), axis=0)
            elif args.method == 'median':
                spec_g = _MEDIAN(specs.astype('float64'), axis=0)
            # limit nan frac
            spec_g[num_finite_chan < num*args.frac_finite_min] = np.nan
            # change spec inplace
            obj.pixel_data[ii] = spec_g
            obj.pixel_specs_nums[ii] = num
            obj.pixel_finite_nums_chan[ii] = num_finite_chan

        np.seterr(**old_set)
        if verbose:
            print('Waiting other CPUs')
#         res =
#         if q is None:
#             # not multiprocessing
#             return res
#         else:
#             q.put(res)


    def __call__(self):
        args = self.args
        self.gen_fpaths()
        self.gen_scale_beams() # also change args.beam_fwhw
        self.set_r_cut()
        self.gen_specfiles()
        self.check_spec_key()
        self.check_axis3()
        self.gen_header_grid()
        self.fill_grid()


        if args.method in ['gaussian', 'bessel_gaussian', 'sinc_gaussian']:
            # convolve
            self.gen_conv_weis()
        elif args.method in ['wmean', 'wmedian']:
            # use gaussian weights
            method_bak = args.method
            args.method = 'gaussian'
            self.gen_conv_weis()
            args.method = method_bak
        elif args.method in ['mean', 'median']:
            # no weights needed
            pass
        else:
            raise(ValueError(f'method {args.method} not support'))

        self._process_sfs()
        # output array only stored the results and don't has other operation, use 'float32'
        self.init_out(DataType='float32')
        # data init as nan not zeros
        self.pixel_data[:] = np.nan
        self.pixel_weis_sum_chan[:] = np.nan

        # determine the file and indices of the spec in
        #
        nspecs = np.array([len(sf.ra) for sf in self.specfiles])
        lens_cum = np.hstack([0, np.cumsum(nspecs)])
        self.ifile =  np.searchsorted(lens_cum, self.ind_cata + 1, side='left') - 1
        self.ind_ifile = self.ind_cata - lens_cum[self.ifile]
        # self.ifile, self.ind_ifile, self.ind_g, self.inds_cata, self.d2d
        #
        # find the start and stop index for each grid
        ind_g_uni = np.unique(self.ind_g)
        start, stop = _get_start_stop(ind_g_uni, self.ind_g)
        # distribute to each process
        ind_g_uni_list = np.array_split(ind_g_uni, args.nproc)
        start_list = np.array_split(start, args.nproc)
        stop_list = np.array_split(stop, args.nproc)
        ps = []
        qs = []
        verbose = True
        for start_, stop_, ind_g_uni_ in zip(start_list, stop_list, ind_g_uni_list):
                q = Queue()
                p = Process(target=self.calc, args=(q, self, start_, stop_, ind_g_uni_, verbose))
                ps += [p]
                qs += [q]
                verbose = False
        for p in ps:
            p.start()
        #self.res = [q.get() for q in qs]
        for p in ps:
            p.join()  # need after q.get()
        #print('gathering results')

        # saving

        _shape = self.ra_grid.shape + self.pixel_data.shape[-1:]
        self.data_cube = self.pixel_data.reshape(_shape).transpose(2,0,1)
        self.weis_chan_cube = self.pixel_weis_sum_chan.reshape(_shape).transpose(2,0,1)
        self.nums_chan_cube = self.pixel_finite_nums_chan.reshape(_shape).transpose(2,0,1)
        self.save()


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
args.share_mem = True

args.key = None
args.wcs_from = '/home/jyj/jingyj/M33/M33_local_areciob.fits'

args.scale_beams_file = None

args.polar = 'M'

args.frac_finite_min = 0.01
