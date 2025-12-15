

__all__ = []


from .utils.io import *


sep_line = '##'+'#'*70+'##'
parser = ArgumentParser(prog=f"python -m hifast.{os.path.basename(sys.argv[0])[:-3]}",
                        formatter_class=formatter_class, allow_abbrev=False,
                        description='Convert file to new format', )
add_common_argument(parser)
parser.add_argument('fpath',
                    help='input hdf5 file path.')


if __name__ == '__main__':
    args = parser.parse_args()
    dict_in = load_hdf5_to_dict_old(args.fpath)

    if 'T' in dict_in.keys():
        wcs_data_name = 'T'
    elif 'Ta' in dict_in.keys():
        wcs_data_name = 'Ta'
    elif 'flux' in dict_in.keys():
        wcs_data_name = 'flux'
    else:
        raise(ValueError('can not find spec'))

    if 'frame' in dict_in['Header'].keys():
        if dict_in[wcs_data_name].ndim == 2:
            dict_in[wcs_data_name] = dict_in[wcs_data_name][..., None]

    dict_in[wcs_data_name] = MjdChanPolar_to_PolarMjdChan(dict_in[wcs_data_name])

    save_specs_hdf5(os.path.join(args.outdir, os.path.basename(args.fpath)), dict_in)
