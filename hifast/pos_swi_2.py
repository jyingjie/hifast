

__all__ = ['IO']


from .utils.io import *


sep_line = '##'+'#'*70+'##'
parser = ArgumentParser(prog=f"python -m hifast.{os.path.basename(sys.argv[0])[:-3]}",
                        formatter_class=formatter_class, allow_abbrev=False,
                        description='Src - Ref', )
add_common_argument(parser)
parser.add_argument('fpath',
                    help='input spectra temperature or flux file path.')
parser.add_argument('--frange', type=float, nargs=2, default=[0, float('inf')],
                    help='Limit frequence range')

parser.add_argument('--t_src', type=float, required=True,
                   help='on-source time[second]')
parser.add_argument('--t_ref', type=float, required=True,
                   help='off-source time[second]')
parser.add_argument('--t_change', type=float, required=True, choices=[30,60],
                    help="switching time[second]. 30s for sepatation of src and ref less than 20'; 60s for sepatation between 20' and 60'.")
parser.add_argument('--n_repeat', type=int, required=True,
                    help='The number of on-source off-source cycles')

parser.add_argument('--only_off', type=bool_fun, choices=[True, False], default='False',
                    help='only noise off')
parser.add_argument('--only_src', type=bool_fun, choices=[True, False], default='False',
                    help='no ref')


class IO(BaseIO):
    ver = 'old'

    def _get_fpart(self,):
        """
        need modify this function
        """
        fpart = '-ps'
        return fpart

    def _import_m(self,):
        super()._import_m()
        global np
        import numpy as np

    def _prepare(self, *args, **kwargs):
        # freq: self.freq_use
        self.mjd_a = self.mjd
        # noise-diode on
        is_on = self.fs['is_on'][:]

        self.inds_on = np.where(is_on)[0]
        self.inds_off = np.where(~is_on)[0]
        self.mjd_on = self.mjd_a[self.inds_on]
        # noise-diode off
        self.mjd_off = self.mjd_a[self.inds_off]
        # load T
        self.p_on = self.s2p[self.inds_on]
        self.p_off = self.s2p[self.inds_off]

#         self.Tcal_s = self.get_Tcal_s()


    @staticmethod
    def in_range(arr, starts, length):
        """
        arr: shape (m,)
        starts: shape (n,)
        length: scalar
        """
        is_in =  (arr[:,None] > starts[None,:]) & (arr[:,None] < (starts+length)[None,:])
        return np.logical_or.reduce(is_in, axis=1)

    def sep(self, t_src, t_ref, n_repeat, t_change=30):
        """
        second
        """
        #######################
        # on: noise-diode on
        # off: noise-diode off
        # src: on-source
        # ref: off-source
        #######################

        # second to day
        t_src = t_src/60/60/24
        t_ref = t_ref/60/60/24
        n_repeat = n_repeat
        t_change = t_change/60/60/24

        mjd_src_start = self.mjd_a[0] + np.arange(n_repeat)*(t_src + t_change + t_ref + t_change)
        mjd_ref_start = mjd_src_start + (t_src + t_change)

        is_ = self.in_range(self.mjd_on, mjd_src_start, t_src)
        self.p_on_src = self.p_on[is_]
        self.inds_on_src = self.inds_on[is_]

        is_ = self.in_range(self.mjd_off, mjd_src_start, t_src)
        self.p_off_src = self.p_off[is_]
        self.inds_off_src = self.inds_off[is_]

        is_ = self.in_range(self.mjd_on, mjd_ref_start, t_ref)
        self.p_on_ref = self.p_on[is_]
        self.inds_on_ref = self.inds_on[is_]

        is_ = self.in_range(self.mjd_off, mjd_ref_start, t_ref)
        self.p_off_ref = self.p_off[is_]
        self.inds_off_ref = self.inds_off[is_]

        self.mjd_src_start = mjd_src_start
        self.mjd_ref_start = mjd_ref_start
        
        
    def gen_integrated_time(self, ):
        args = self.args
        
        tint_src = np.sum(~np.all(np.isnan(self.p_off_src[:,:,-1]), axis = 1))
        tint_ref = np.sum(~np.all(np.isnan(self.p_off_ref[:,:,-1]), axis = 1))
        
        if not args.only_off:
            tint_src += np.sum(~np.all(np.isnan(self.p_on_src[:,:,-1]), axis = 1))
            tint_ref += np.sum(~np.all(np.isnan(self.p_on_ref[:,:,-1]), axis = 1))
        
        from astropy.time import Time
        mjd = Time(self.mjd, format='mjd')
        delta_mjd = (mjd[1] - mjd[0]).sec
        
        self.tint_src = tint_src * delta_mjd
        self.tint_ref = tint_ref * delta_mjd

    def gen_Ta(self, only_off=False):
        """
        using

        self.p_on_src
        self.p_off_src
        self.p_on_ref
        self.p_off_ref

        to gen

        self.Ta as shape (1,chan,polar)

        """
        args = self.args

        Ta_on_src = np.nanmean(self.p_on_src, axis=0, keepdims=True)
        Ta_off_src = np.nanmean(self.p_off_src, axis=0, keepdims=True)
        Ta_on_ref = np.nanmean(self.p_on_ref, axis=0, keepdims=True)
        Ta_off_ref = np.nanmean(self.p_off_ref, axis=0, keepdims=True)

        # use interaged time as weight
        if not args.only_off:
            Ta_src = np.average(np.r_[Ta_on_src, Ta_off_src], axis=0, weights=[len(self.inds_on_src),len(self.inds_off_src)])
            Ta_ref = np.average(np.r_[Ta_on_ref, Ta_off_ref], axis=0, weights=[len(self.inds_on_ref),len(self.inds_off_ref)])
        else:
            Ta_src = Ta_off_src[0]
            Ta_ref = Ta_off_ref[0]

        self.Ta_src = Ta_src[None,]
        self.Ta_ref = Ta_ref[None,]

        if args.only_src:
            self.Ta = self.Ta_src
        else:
            self.Ta = self.Ta_src - self.Ta_ref

        # also change mjd
        inds_tmp = self.inds_off_src[:1]
        self.mjd = self.mjd_a[inds_tmp]

    def gen_radec(self,):
        """
        try to calculation for each record: self.ra_a, self.dec_a
        use the position of first src record as self.ra, self.dec
        also set self.mjd and self.mjd_a

        """

        self.ra_a, self.dec_a = self.ra, self.dec
        inds_tmp = self.inds_off_src[:1]
        self.ra, self.dec = self.ra_a[inds_tmp], self.dec_a[inds_tmp]
        self.mjd = self.mjd_a[inds_tmp]

    def plot_radec(self, figsize=[10,10], outname=None):
        """
        plot the ra dec to check the separation of src and ref
        """

        from matplotlib import pyplot as plt
        plt.figure(figsize=figsize)
        plt.scatter(self.ra_a[self.inds_on_src], self.dec_a[self.inds_on_src], s=1, color='r')
        plt.scatter(self.ra_a[self.inds_off_src], self.dec_a[self.inds_off_src], s=1, color='r', label='src')
        plt.scatter(self.ra_a[self.inds_on_ref], self.dec_a[self.inds_on_ref], s=1, color='b')
        plt.scatter(self.ra_a[self.inds_off_ref], self.dec_a[self.inds_off_ref], s=1, color='b', label='ref')
        plt.xlabel('ra')
        plt.ylabel('dec')
        plt.grid()
        plt.minorticks_on()
        plt.legend()
        if outname is not None:
            print(f'Saving ra dec plot to {outname}')
            plt.savefig(outname)

    def gen_s2p_out(self,):
        args = self.args
        self._prepare()
        self.sep(args.t_src, args.t_ref, args.n_repeat, args.t_change)
        self.gen_integrated_time()
        self.gen_Ta()
        self.gen_radec()
        self.plot_radec(outname=self.fpath_out + '-radec.png')
        self.s2p_out = self.Ta

    def __call__(self, save=True):
        # store output data in self.dict_out
        self.gen_s2p_out()
        self.gen_dict_out()
        for key in ['is_delay', 'is_on', 'next_to_cal']:
            if key in self.dict_out.keys():
                 self.dict_out[key] = np.array([False,])
        self.dict_out['Ta_src'] = MjdChanPolar_to_PolarMjdChan(self.Ta_src)
        self.dict_out['Ta_ref'] = MjdChanPolar_to_PolarMjdChan(self.Ta_ref)
        
        
        for key in ['tint_src', 'tint_ref']:
            print(key, getattr(self, key))
            self.dict_out[key] = np.array([getattr(self, key)])
        
        # save to hdf5 file
        if save:
            self.save()


if __name__ == '__main__':
    args_ = parser.parse_args()
    # print(parser.format_help())
    # print("----------")
    print('#'*35+'Args'+'#'*35)
    args_from = parser.format_values()
    print(args_from)
    print('#'*35+'####'+'#'*35)

    HistoryAdd = {'args_from': args_from} if args_.my_config is not None else None
    io = IO(args_, HistoryAdd=HistoryAdd)
    io()
