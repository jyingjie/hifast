

__all__ = ['add_ness_sdfits_header', 'IO']


from .utils.io import *


sep_line = '##'+'#'*70+'##'
parser = ArgumentParser(prog=f"python -m hifast.{os.path.basename(sys.argv[0])[:-3]}", formatter_class=formatter_class, allow_abbrev=False,
                        description='Convert the hifast spectra file to sdfits format', )
add_common_argument(parser)
parser.add_argument('fpath',
                    help='input spectra temperature file path.')
parser.add_argument('--frange', type=float, nargs=2, default=[0, float('inf')],
                    help='Limit frequence range')
parser.add_argument('--no_radec', action='store_true', not_in_write_out_config_file=True,
                    help="if set, don't check or add ra dec")
parser.add_argument('--polar', choices=['XX', 'YY', 'All'], default='All',
                   help='Chose polarization')


def add_ness_sdfits_header(header):
    """
    add necessary keys to the header of a SDFITS Table
    https://fits.gsfc.nasa.gov/registry/sdfits.html
    """

    items = {
    'XTENSION': 'BINTABLE', # FITS binary table
    'BITPIX' : 8, # MANDATORY-- Binary data
    'NAXIS': 2, # MANDATORY-- A 2D TABLE
    'NAXIS1': -1, # MANDATORY-- width of table in bytes
    'NAXIS2': -1, # MANDATORY-- Number of rows in table
    'PCOUNT': 0, # MANDATORY-- Recommended value 0, no heap
    'GCOUNT': 1, # MANDATORY--
    'TFIELDS': 1, # MANDATORY-- Number of fields per row
    'EXTNAME': 'SINGLE DISH', # SINGLE DISH
    'EXTVER': 1, # ENTIRELY OPTIONAL
    'EXTLEVEL': 1, # OPTIONAL SHOULD DEFAULT TO 1
    'NMATRIX': 1, # MANDATORY--  MATRIX PER ROW
    'TELESCOP': 'FAST',
    # 'SPECSYS': 'LSRK',
    # 'SSYSOBS': 'TOPOCENT',
    # 'EQUINOX': 2000.0,
    # 'RADESYS': 'FK5',
    }
    ## 'NAXIS1', 'NAXIS2' need set later
    for key in items.keys():
        header[key] =  items[key]


class IO(BaseIO):
    ver = 'old'

    def _get_fpart(self,):
        """
        need modify this function
        """
        fpart = '-to'
        return fpart
    def _gen_fpath_out(self):
        super()._gen_fpath_out( )
        self.fpath_out = os.path.splitext(self.fpath_out)[0] + '.sdfits'

    def init_table(self,):
        """ initialize the binary table"""
        import numpy as np
        from astropy.io import fits

        dtype_descr = [('SCAN', '<i2'),
                        ('CYCLE', '<i4'),
                        ('DATE-OBS', '|S10'),
                        ('TIME', '<f8'),
                        ('EXPOSURE', '<f4'),
                        ('OBJECT', '|S16'),
                        ('OBJ-RA', '<f8'),
                        ('OBJ-DEC', '<f8'),
                        ('RESTFRQ', '<f8'),
                        ('OBSMODE', '|S16'),
                        ('BEAM', '<i2'),
                        ('IF', '<i2'), # need for aoflagger
                        ('FREQRES', '<f8'),
                        ('BANDWID', '<f8'),
                        ('CRPIX1', '<f4'),
                        ('CRVAL1', '<f8'),
                        ('CDELT1', '<f8'),
                        ('CRPIX3', '<f4'),
                        ('CRVAL3', '<f8'),
                        ('CDELT3', '<f8'),
                        ('SCANRATE', '<f4', (2,)),
                        ('TSYS', '<f4', (2,)),
                        ('CALFCTR', '<f4', (2,)),
                        ('DATA', '<f4', (1, 1, *self.s2p.shape[1:][::-1]),),
                        ('FLAGGED', '|u1', (1, 1, *self.s2p.shape[1:][::-1]),),
                        ('TCAL', '<f4', (2,)),
                        ('TCALTIME', '|S16'),
                        ('AZIMUTH', '<f4'),
                        ('ELEVATIO', '<f4'),
                        ('PARANGLE', '<f4'),
                        ('FOCUSAXI', '<f4'),
                        ('FOCUSTAN', '<f4'),
                        ('FOCUSROT', '<f4'),
                        ('TAMBIENT', '<f4'),
                        ('PRESSURE', '<f4'),
                        ('HUMIDITY', '<f4'),
                        ('WINDSPEE', '<f4'),
                        ('WINDDIRE', '<f4')]

        dtype = np.dtype(dtype_descr)
        self.table = fits.BinTableHDU.from_columns(np.zeros(len(self.s2p), dtype=dtype))

    def __call__(self):
        # determine the polarization
        if self.args.polar == 'XX':
            self.s2p = self.s2p[...,0:1]
        elif self.args.polar == 'YY':
            self.s2p = self.s2p[...,1:2]
        elif self.args.polar == 'All':
            pass

        from astropy.io import fits

        hdulist = fits.HDUList()
        # Create a new primary HDU
        primary_hdu = fits.PrimaryHDU()
        # Add the primary and table HDUs to the HDU list
        hdulist.append(primary_hdu)
        self.init_table()
        self.table.data['DATA'][:] = self.s2p.transpose(0,2,1)[:, None, None]

        # set CTYPE1
        self.table.data['CRPIX1'][:]  = 1
        self.table.data['CDELT1'][:]  = (self.freq[1] - self.freq[0]) * 1e6
        self.table.data['CRVAL1'][:]  = self.freq[0] * 1e6
        self.table.data['FREQRES'][:]  = 1420.405751768 * 1e6
        self.table.data['BANDWID'][:]  = 500 * 1e6 # not sure

        # self.table.data['CRPIX3'][:]  = 1
        # self.table.data['CDELT3'][:]  = 0.033
        # self.table.data['CRVAL3'][:]  = 1

        hdulist.append(self.table)
        # Add metadata to the primary HDU
        pass
        # Add metadata to the table HDU
        add_ness_sdfits_header(hdulist[1].header)
        ## change some items in the header
        hdulist[1].header['NAXIS1'] = self.table.data.nbytes // len(self.table.data)
        hdulist[1].header['NAXIS2'] = len(self.table.data)

        hdulist[1].header['CTYPE1']  = 'FREQ'

        hdulist[1].header['CTYPE2']  = 'STOKES'
        ## not sure
        # hdulist[1].header['CRPIX2']  = 1
        # hdulist[1].header['CDELT2']  = -5
        # hdulist[1].header['CRVAL2']  = -1

        hdulist[1].header['CTYPE3']  = 'RA'

        hdulist[1].header['CTYPE4']  = 'DEC'
        ## not sure
        # hdulist[1].header['CRPIX4']  = 1
        # hdulist[1].header['CDELT4']  = 0.02
        # hdulist[1].header['CRVAL4']  = 1

        self.hdulist = hdulist
        # Write the SDFITS file to disk
        hdulist.writeto(self.fpath_out, overwrite=self.args.force)


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
