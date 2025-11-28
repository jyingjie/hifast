

__all__ = ['sep_line', 'parser', 'IO']



from .utils.io import *

sep_line = '##'+'#'*70+'##'
parser = ArgumentParser(prog=f"python -m hifast.{os.path.basename(sys.argv[0])[:-3]}",
                        formatter_class=formatter_class, allow_abbrev=False,
                        description='T Cal, Src - Ref', )
parser.add_argument('fpath',
                    help='only need input the first chunk file path (e.g. XXX_0001.fits) of RAW spectra data')
parser.add_argument('-f', dest='force', action='store_true',
                    help='overwriting file if out file exists')
parser.add_argument('--outdir', required=True,
                   help='the directory to store output files.')
parser.add_argument('-g', is_write_out_config_file_arg=True,
                        help='save config to file path')
parser.add_argument('-c', '--my-config', is_config_file_arg=True,
                        help='config file path')
parser.add_argument('--frange', type=float, nargs=2,
                    help='freq range')


parser.add_argument('-d', '--n_delay', type=int, required=True,
                   help='time of delay divided by sampling time')
parser.add_argument('-m', '--n_on', type=int, required=True,
                   help='time of Tcal_on divided by sampling time')
parser.add_argument('-n', '--n_off', type=int, required=True,
                    help='time of Tcal_off divided by sampling time')

parser.add_argument('--smooth', choices=['mean','poly','gaussian'], required=True,
                    help="smooth method: 'mean','poly',gaussian',.. default:gaussian")
parser.add_argument('--s_sigma', type=float, default=5,
                    help='sigma for gaussian smooth, default 5MHz')

parser.add_argument('--noise_mode', default='high', choices=['high','low'],
                    help='noise_mode, high or low')
parser.add_argument('--noise_date', default='auto',
                    help='noise obs date, default auto')
parser.add_argument('--med_filter_size_cal', type=int, default=5,
                    help='median filter kernel size for power of cal; odd number; default 5')


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



class IO(Path_IO):
    def __init__(self, args, inplace_args=False):
        """
        args: class
              including attributes: fpath, outdir, frange
        """
        self.args = args if inplace_args else copy.deepcopy(args)
        self._gen_fpath_out()
        self._check_fout()
        self._import_m()
        self.nB = self.get_nB(self.args.fpath)
        self.load_and_add_Header()

    def _import_m(self,):
        global np
        import numpy as np

    def _gen_fpath_out(self,):
        """
        """
        args = self.args

        # replace patten in outdir
        nB = self.get_nB(args.fpath)
        project = get_project(args.fpath)
        ## use the dirname of the fits file as "date"
        date = os.path.basename(os.path.dirname(os.path.abspath(args.fpath)))
        args.outdir = sub_patten(args.outdir, date=date, nB=f'{nB:02d}', project=project)
        args.outdir = os.path.expanduser(args.outdir)
        print(f'outdir: {args.outdir}')
        if not os.path.exists(args.outdir):
            print(f'outdir {args.outdir} not exists. Create it now')
            os.makedirs(args.outdir, exist_ok=True)

        fname_part = re.sub('[0-9]{4}\.fits\Z', '', args.fpath)
        fname_add = os.path.basename(os.path.dirname(os.path.abspath(fname_part)))
        out_name_base = os.path.join(args.outdir, f"{os.path.basename(fname_part)[:-1]}-{fname_add}")
        self.out_name_base = out_name_base
        fpart = '-S_ps'
        self.fpath_out = f'{out_name_base}{fpart}.hdf5'

    def load_and_add_Header(self, ):
        import json
        Header = {}
        Header.update(rec_his(args=json.dumps(self.args.__dict__)))
        self.Header = Header

    def gen_S(self,):

        from .core.cal import PositionSwitch

        args = self.args
        S = PositionSwitch(args.fpath, args.n_delay, args.n_on, args.n_off,
                 frange=args.frange, verbose=True,
                 smooth=args.smooth, s_para={'s_sigma':args.s_sigma},
                 noise_mode=args.noise_mode, noise_date=args.noise_date,
                 med_filter_size_cal=args.med_filter_size_cal)

        S.sep(args.t_src, args.t_ref, args.n_repeat, args.t_change)
        # S.gen_out_name_base(args.outdir)

        self.S = S

    def __call__(self,):
        args = self.args
        self.gen_S()
        S = self.S

        S.plot = True
        S.out_name_base = self.out_name_base
        S.plot_sep()

        S.gen_Ta(only_off=args.only_off)
        S.gen_radec()
        S.plot_radec(outname=self.fpath_out + '-radec.png')

        dict_out = {}
        dict_out['Header'] = self.Header
        dict_out['Ta'] = MjdChanPolar_to_PolarMjdChan(S.Ta)
        dict_out['freq'] = S.freq_use

        dict_out['ra'] = S.ra
        dict_out['dec'] = S.dec
        dict_out['mjd'] = S.mjd

        # for compatible

        dict_out['is_delay'] = np.array([False,])
        dict_out['is_on'] = np.array([False,])
        dict_out['next_to_cal'] = np.array([False,])
        dict_out['Tcal'] = S.Tcal_s


        self.outfield = 'Ta'
        self.dict_out = dict_out

        print("Saving...")
        save_specs_hdf5(self.fpath_out, self.dict_out, wcs_data_name=self.outfield)
        print(f"Saved to {self.fpath_out}")



if __name__ == '__main__':
    args_ = parser.parse_args()

    # print(parser.format_help())
    # print("----------")
    # print(parser.format_values())  # useful for logging where different settings came from
    print('#'*35+'Args'+'#'*35)
    print(parser.format_values())  # useful for logging where different settings came from
    print('#'*35+'####'+'#'*35)
    io = IO(args_)
    io()
