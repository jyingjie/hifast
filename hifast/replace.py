from .utils.io import *

sep_line = '##'+'#'*70+'##'
parser = ArgumentParser(prog=f"python -m hifast.{os.path.basename(sys.argv[0])[:-3]}", formatter_class=formatter_class, allow_abbrev=False,
                        description='replace fields in fpath_ori with that from fpath_new', )
parser.add_argument('fpath_ori',
                    help='file path')
parser.add_argument('fpath_new',
                    help='file path')
parser.add_argument('--fields', nargs='+', required=True,
                    help='fields will be replaced, e.g. mjd')
parser.add_argument('-y', '--yes', action='store_true', not_in_write_out_config_file=True,
                    help="Do not ask for confirmation.")


if __name__ == '__main__':
    import h5py
    args_ = parser.parse_args()
    # print(parser.format_help())
    # print("----------")
    print('#'*35+'Args'+'#'*35)
    print(parser.format_values())  # useful for logging where different settings came from
    print('#'*35+'####'+'#'*35)
    args = args_

    f_ori =  h5py.File(args.fpath_ori, 'r+')
    S_ori = f_ori['S']
    f_new =  h5py.File(args.fpath_new, 'r')
    S_new = f_new['S']

    print(f"""replace {','.join(args.fields)} in
{args.fpath_ori}
with that in
{args.fpath_new}
""")
    if not args.yes:
        QA = input('input y[es] to continue\n')
        if not (QA.lower() == 'y' or QA.lower() == 'yes'):
            print('aborting...')
            f_ori.close()
            f_new.close()
            sys.exit(0)
    for key in args.fields:
        print(f"replacing {key}")
        S_ori[key][...] = S_new[key][:]

    his = rec_his()
    for key in his.keys():
        if 'HISTORY' in key:
            pass
            f_ori['Header'].attrs[key] = his[key]

    f_ori.close()
    f_new.close()
