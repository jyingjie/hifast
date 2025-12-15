

__all__ = ['gen_radec_file', 'get_radec']


import os
import re
import sys
from .utils.io import get_nB, replace_nB


def gen_radec_file(file_spec, paras=[]):
    """
    generating ra dec for a specs_T file
    """
    import subprocess
    command = [sys.executable, '-m', 'hifast.radec', file_spec] + paras
    print('run:')
    print(' '.join(command))
    sys.stdout.flush()
    process = subprocess.Popen(
        command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # process.wait()
    print(*process.communicate())
    if process.returncode != 0:
        raise(ValueError(f'fail to generate the radec file of {file_spec}'))


import h5py

def get_radec(file_spec, mjd, nB_radec=1, file_radec = None):
    """
    load the ra dec of spectra in file_spec (*specs_T.hdf5)

    file_spec: str; specs fpath, used to guess radec filename
    mjd: np.array; used to check

    """
    nB = get_nB(file_spec)
    if file_radec is None:
        # ra dec file
        file_radec = file_spec.rsplit('specs_T', 1)[0] + 'specs_T-radec.hdf5'
        #file_radec = re.findall(r'.*-specs_T', file_radec)+'-radec.hdf5'
        file_radec_ori = file_radec
        # default is in beam `nB_radec`
        file_radec = replace_nB(file_radec, nB_radec)

    if os.path.exists(file_radec):
        f = h5py.File(file_radec, 'r')
    elif os.path.exists(file_radec_ori):
        f = h5py.File(file_radec_ori, 'r')
    else:
        raise(OSError(f"can not find the RA DEC file \n{file_radec}\n, please generate it by using 'python -m hifast.radec'"))

    # check if mjd match
    mjd_match = False if len(mjd) != len(
        f['S']['mjd'][:]) else (mjd == f['S']['mjd'][:]).all()
    if not mjd_match:
        print('the mjd of spec file is not same with that in ', file_radec)
        sys.stdout.flush()
        file_radec_2 = replace_nB(file_radec, nB)
        if not os.path.exists(file_radec_2):
            print('try to generate ', file_radec_2)
            sys.stdout.flush()
            header = f['Header'].attrs
            import json
            for key in header.keys():
                argv = json.loads(header[key])['argv']
                paras = argv.split()
                if os.path.basename(paras[0]) not in ['cli_radec.py', 'radec.py']:
                    continue
                else:
                    paras.pop(0)
                    for s in paras:
                        if '.hdf5' in s:
                            paras.remove(s)
            gen_radec_file(file_spec, paras)
            sys.stdout.flush()
        f.close()
        f = h5py.File(file_radec_2, 'r')
        print('using', file_radec_2)
        sys.stdout.flush()
    # check again
    mjd_match = False if len(mjd) != len(
        f['S']['mjd'][:]) else (mjd == f['S']['mjd'][:]).all()
    if not mjd_match:
        raise(ValueError(
            f'the mjd of spec file is not same with that in {file_radec_2} Abort...'))
    # load
    ra = f['S']['ra'+'%d' % nB][:]
    dec = f['S']['dec'+'%d' % nB][:]
    if 'is_extrapo' in f['S'].keys():
        is_extrapo = f['S']['is_extrapo'][:]
    else:
        is_extrapo = None
    f.close()
    return ra, dec, is_extrapo


if __name__ == '__main__':
    import h5py
    import argparse
    parser = argparse.ArgumentParser(allow_abbrev=False)
    parser.add_argument('fname',
                        help='file name')
#     parser.add_argument('-f', '--force', action='store_true',
#                         help='overwriting file if out file exists')
    parser.add_argument('--nB_radec', type=int, default=1,
                        help='Beam number of the radec file name')
    parser.add_argument('--file_radec',
                        help='file radec name if it exists')
#     parser.add_argument('--outdir',
#                         help='default is same with the input file')
    args = parser.parse_args()
    f = h5py.File(args.fname, 'r+')
    g = f['S']

    values = get_radec(args.fname, g['mjd'][:], args.nB_radec, args.file_radec)

    keys = ['ra', 'dec', 'is_extrapo']

    for key, value in zip(keys, values):
        if value is None:
            continue
        if key in g.keys():
            g[key][:] = value
        else:
            g[key] = value
    f.close()
