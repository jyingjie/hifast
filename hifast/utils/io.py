

__all__ = ['formatter_class', 'os', 'sys', 'c', 're', 'copy', 'ArgumentParser', 'argparse', 'sub_patten', 'bool_fun',
           'add_common_argument', 'hide_paras', 'del_paras_in_string', 'rec_his', 'MjdChanPolar_to_PolarMjdChan',
           'PolarMjdChan_to_MjdChanPolar', 'save_dict_hdf5', 'gen_carta_group', 'save_specs_hdf5', 'load_hdf5_to_dict',
           'load_hdf5_to_dict_old', 'get_project', 'get_date_from_path', 'get_nB', 'replace_nB', 'Path_IO', 'BaseIO',
           'HFDataT', 'HFGroup', 'H5HDU', 'H5FitsRead']


import os
import sys
import re
import copy
from .colors import Colors as c
from .conf_arg import ArgumentParser, argparse
# from configargparse import ArgumentParser
# import argparse
formatter_class = argparse.ArgumentDefaultsHelpFormatter


#nbdev_comment _all_ = ['os', 'sys', 'c', 're', 'copy', 'ArgumentParser', 'argparse']


def sub_patten_1(string, **kwargs):
    """
    substitute '%(...)s' in string with value in kwargs
    """
    keys = re.findall(r"%\(([^)]+)\)s", string)
    for key in keys:
        if key not in kwargs.keys():
            raise ValueError(f'can not find value to replace %({key})s')
        string = re.sub(r"%\("+ key + r"\)s", str(kwargs[key]), string)
    return string

def sub_patten_2(string, **kwargs):
    """
    substitute '%[...]s' in string with value in kwargs
    """
    keys = re.findall(r"%\[(.*?)\]s", string)
    for key in keys:
        if key not in kwargs.keys():
            raise ValueError(f'can not find value to replace %[{key}]s')
        string = re.sub(r"%\[" + key + r"\]s", str(kwargs[key]), string)
    return string




# export
def sub_patten(string, **kwargs):
    if '%(' in string:
        return sub_patten_1(string, **kwargs)
    elif '%[' in string:
        return sub_patten_2(string, **kwargs)
    else:
        return string


def bool_fun(s):
    """
    used as the type of parser.add_argument
    """
    s = str(s)
    if s == 'True' or s.lower() == 'yes':
        return True
    elif s == 'False' or s.lower() == 'no':
        return False
    else:
        raise(ValueError("input must be 'True'('yes') or 'False'('no')"))


def add_common_argument(parser):
    # common
    parser.add_argument('--outdir', default='default',
                        help='The directory used to save output file, default is same with the input file')
    parser.add_argument('-f', action='store_true', dest='force',
                        help='if set, overwriting file if output file exists')
    parser.add_argument('-g', is_write_out_config_file_arg=True,
                        help='save config to file path')
    parser.add_argument('-c', '--my-config', is_config_file_arg=True,
                        help='config file path')


def hide_paras(parser, dests, also_not_in_write_out_config_file=True):
    """
    hiding parameters in dests list from the help message for clarity and testing purpose
    """

    for g in parser._action_groups:
        is_ = []
        for a in g._actions:
            if a.dest in dests:
                a.help = argparse.SUPPRESS
                if also_not_in_write_out_config_file:
                    a.not_in_write_out_config_file = True
                is_ += [True]
            else:
                is_ += [False]
        if len(is_) !=0 and len(is_) == sum(is_):
            g.title = ''

def del_paras_in_string(s, dests):
    """
    s: string
    dests: list of string
    """
    import re
    for dest in dests:
        s = re.sub(re.compile(rf'^.*{dest}:.*\n', re.M), '', s)
    return s


def rec_his(**kwargs):
    """
    return OrderedDict {'HISTORY+{time}': '{"key:value" in kwargs}'}
    """
    import json
    from datetime import datetime
    from collections import OrderedDict
    try:
        from ..__init__ import __version__
    except:
        __version__ = 'unknown'

    history = OrderedDict()
    history['version'] = __version__
    history['cwd'] = os.getcwd()
    history['argv'] = ' '.join(sys.argv)
    for key in kwargs.keys():
        history[key] = kwargs[key]
    current_time = datetime.now().strftime("%Y%m%d-%H:%M:%S")
    return OrderedDict({"HISTORY-"+current_time: json.dumps(history, indent=2)})


def add_extra(fin, _dict=None, fields_add=[]):
    """
    add some fields of fin to _dict
    """
    fields = ['is_on', 'next_to_cal', 'is_delay', 'Tcal', 'is_extrapo', 'vel']
    fields += fields_add
    out_add = {}
    for field in fields:
        if field in fin.keys():
            try:
                out_add[field] = fin[field][:]
            except:
                pass
    if _dict is None:
        return out_add
    else:
        _dict.update(out_add)


def MjdChanPolar_to_PolarMjdChan(mcp, check_polar=True):
    if check_polar and mcp.shape[2] != 2 and mcp.shape[2] != 1:
        raise(ValueError('need input (Mjd,Chan,Polar)'))
    return mcp.transpose((2, 0, 1))


def PolarMjdChan_to_MjdChanPolar(pmc, check_polar=True):
    if check_polar and pmc.shape[0] != 2 and pmc.shape[0] != 1:
        raise(ValueError('need input (Polar,Mjd,Chan)'))
    return pmc.transpose((1, 2, 0))


def save_dict_hdf5(fname, dict_in, mode='w', spec2float32=False):
    """
    fname: str
    dict_in: dict; keys of the dict_in are str, value needs support slice, i.e. value[:].
             if 'Header':{...} in dict_in, save in attrs of group 'Header'
    spec2float32: default: True
             if True, convert 'Ta', 'flux', 'T' to float32
    """
    import h5py
    import numpy as np
    try:
        f = h5py.File(fname, mode)
    except OSError:
        if mode == 'w':
            from datetime import datetime
            os.rename(fname, fname+'.del.empty.'+datetime.now().strftime("%Y%m%d-%H%M%S"))
            f = h5py.File(fname, mode)
    for key in dict_in.keys():
        if key == 'Header':
            # sometimes hdf5 raises error if track_order = True
            header = dict_in['Header']
            f.create_group('Header', track_order=False)
            for key2 in header.keys():
                f['Header'].attrs[key2] = header[key2]
        else:
            if isinstance(dict_in[key], h5py.Dataset):
                f[key] = dict_in[key][:]
            elif spec2float32 and isinstance(dict_in[key], np.ndarray) and key in ['Ta', 'flux', 'T', 'Power']:
                f[key] = dict_in[key].astype('float32')
            else:
                f[key] = dict_in[key]
    f.close()


def gen_carta_group(f, data_shape, wcs_data_name, axis1=None, axis2=None,
                    wcs_data_group='S',
                    image_group='Waterfall'):
    import h5py
    import numpy as np
    wcs = {'SIMPLE': 1,
           'WCSAXES': len(data_shape),
           'NAXIS': len(data_shape),
           'NAXIS1': data_shape[2],
           'NAXIS2': data_shape[1],
           'NAXIS3': data_shape[0],
           'CDELT1': 1.,
           'CDELT2': 1.,
           'CDELT3': 1.,
           'CRPIX1': 1,
           'CRPIX2': 1,
           'CRPIX3': 1,
           'CRVAL1': 1420.,
           'CRVAL2': 1.,
           'CRVAL3': 1.,
           'CTYPE1': np.bytes_('OFFSET'),
           'CTYPE2': np.bytes_('MJD-OBS'),
           'CTYPE3': np.bytes_('POLAR'),
           'CUNIT1': np.bytes_('MHz'),
           'CUNIT2': np.bytes_('d'),
           'CUNIT3': np.bytes_('s'),
           'EQUINOX': 2000.0,
           # 'LINE': b'HI 1420 1420.555',
           # 'BITPIX': -32,
           # 'LATPOLE': 0.0,
           # 'LONPOLE': 0.0,
           }
    try:
        # CDELT * NAXIS affect the height to width ratio of the image if AXISES are same type
        wcs['CRVAL1'] = axis1[0]
        wcs['CDELT1'] = axis1[1] - axis1[0]
        wcs['CRVAL2'] = axis2[0]
        wcs['CDELT2'] = axis2[1] - axis2[0]
        # delts = f"{dict_in['freq'][1] - dict_in['freq'][0]}" +\
        #         f"|{dict_in['mjd'][1] - dict_in['mjd'][0]}"
        # wcs['LINE'] = np.bytes_(delts)
        # print(wcs)
    except:
        pass
    f.create_group(image_group)
    f[image_group]['DATA'] = h5py.SoftLink(f'/{wcs_data_group}/{wcs_data_name}')
    # for key in ['MipMaps', 'PermutedData', 'Statistics', 'SwizzledData']:
    #     f['0'].create_group(key)
    for key in wcs.keys():
        f[image_group].attrs[key] = wcs[key]


def save_specs_hdf5(fname, dict_in, mode='w', spec2float32=True, wcs_data_name=None):
    """
    fname: str
    dict_in: dict; keys of the dict_in are str, value needs support slice, i.e. value[:].
             if 'Header':{...} in dict_in, save in attrs of group 'Header'
    spec2float32: default: True
             if True, convert 'Ta', 'flux', 'T' to float32
    """
    import h5py
    import numpy as np
    # handle hdf5 is opened
    try:
        f = h5py.File(fname, mode)
    except OSError:
        if mode == 'w':
            from datetime import datetime
            os.rename(fname, fname+'.del.empty.'+datetime.now().strftime("%Y%m%d-%H%M%S"))
            f = h5py.File(fname, mode)
    # save items in dict_in to group 'S'
    f.create_group('S')
    for key in dict_in.keys():
        if key == 'Header':
            # sometimes hdf5 raises error if track_order = True
            header = dict_in['Header']
            f.create_group('Header', track_order=False)
            for key2 in header.keys():
                f['Header'].attrs[key2] = header[key2]
        else:
            if isinstance(dict_in[key], h5py.Dataset):
                f['S'][key] = dict_in[key][:]
            elif spec2float32 and isinstance(dict_in[key], np.ndarray) and key in ['Ta', 'flux', 'T', 'Power']:
                f['S'][key] = dict_in[key].astype('float32')
            else:
                f['S'][key] = dict_in[key]
    # save group for carta read
    # try:

    if wcs_data_name is None:
        if 'T' in dict_in.keys():
            wcs_data_name = 'T'
        elif 'Ta' in dict_in.keys():
            wcs_data_name = 'Ta'
        elif 'flux' in dict_in.keys():
            wcs_data_name = 'flux'
        else:
            wcs_data_name = 'none'
            # not data to gen_carta_group, return
    if wcs_data_name == 'none':
        f.close()
        return
    try:
        axis1 = dict_in['freq'] if 'freq' in dict_in.keys() else None
        axis2 = dict_in['mjd'] if 'mjd' in dict_in.keys() else None
        if not isinstance(dict_in[wcs_data_name], h5py.ExternalLink):
            data_shape = dict_in[wcs_data_name].shape
        else:
            # beacuse ExternalLinked file has been opend as 'r', need use the same mode ('r') to open fname to access h5py.ExternalLink
            f.close()
            f = h5py.File(fname, 'r')
            data_shape = f['S'][wcs_data_name].shape
            f.close()
            # reopen as 'r+'
            f = h5py.File(fname, 'r+')
        gen_carta_group(f, data_shape, wcs_data_name, axis1, axis2)
    except Exception as Err:
        print(Err)
    f.close()


def load_hdf5_to_dict(fpath):
    """
    fpath: str
    """
    import h5py
    from collections import OrderedDict
    dict_out = {}
    fs = h5py.File(fpath, 'r')
    # load Header
    Header = OrderedDict()
    if 'Header' in fs.keys():
        try:
            Header.update(OrderedDict(fs['Header'].attrs.items()))
        except:
            print('input file have no Header')
    dict_out['Header'] = Header
    # load
    for key in fs['S'].keys():
        if key != 'Header':
            dict_out[key] = fs['S'][key][:]
    fs.close()
    return dict_out


def load_hdf5_to_dict_old(fpath):
    """
    fpath: str
    """
    import h5py
    from collections import OrderedDict
    dict_out = {}
    fs = h5py.File(fpath, 'r')
    # load Header
    Header = OrderedDict()
    if 'Header' in fs.keys():
        try:
            Header.update(OrderedDict(fs['Header'].attrs.items()))
        except:
            print('input file have no Header')
    dict_out['Header'] = Header
    # load
    for key in fs.keys():
        if key != 'Header':
            dict_out[key] = fs[key][:]
    fs.close()
    return dict_out


def get_project(path):
    """
    return project name deduced from "path"
    """
    basename = os.path.basename(path)
    ind = re.search(r'[0-9][0-1]M-', basename[::-1]).end()
    return basename[:len(basename)-ind]


def get_date_from_path(path):
    """
    e.g. input path: 'data/G15_6_arcdrift-M01_W-20210808-specs_T.hdf5'
         return '20210808'
    """
    basename = os.path.basename(path)
    ind = re.search(r'[0-9][0-1]M-', basename[::-1]).start()
    return basename[-ind:].split('-')[1]


def get_nB(path):
    """
    return beam number deduced from "path"
    """
    return int(re.findall(r'-M[0-1][0-9]', os.path.basename(path))[-1][2:])


def replace_nB(path, nB):
    """
    replace the beam number of 'path' with the input nB
    path: str
    nB: int
    """
    return re.sub(r'[0-9][0-1]M-', f"-M{nB:02d}"[::-1], path[::-1], count=1)[::-1]


class Path_IO(object):
    get_nB = staticmethod(get_nB)
    replace_nB = staticmethod(replace_nB)

    def __init__(self, args, inplace_args=False):
        self.args = args if inplace_args else copy.deepcopy(args)

    def _check_fout(self,):
        """
        check whether self.fout exists
        """
        args = self.args
        if os.path.exists(self.fpath_out):
            if args.force:
                print(f"will overwrite the existing output file {self.fpath_out}")
            else:
                print(f"File exists {self.fpath_out}")
                print("exit... Use ' -f ' to overwrite it.")
                sys.exit(0)

    def _gen_fpath_out(self,):
        """
        add self.fpath_out
        """
        args = self.args
        nB = get_nB(args.fpath)
        project = get_project(args.fpath)
        date = get_date_from_path(args.fpath)
        self.nB = nB
        self.project = project
        self.date = date

        if args.outdir is None or args.outdir == 'default':
            args.outdir = os.path.dirname(args.fpath)
        else:
            # replace patten in outdir
            args.outdir = sub_patten(args.outdir, date=date, nB=f'{nB:02d}', project=project)
            # expand '~' as outdir may be a string in bash
            args.outdir = os.path.expanduser(args.outdir)

        print(f'outdir: {args.outdir}')
        if args.outdir!='' and (not os.path.exists(args.outdir)):
            print(f'outdir {args.outdir} not exists. Create it now')
            os.makedirs(args.outdir, exist_ok=True)

        self.fpath_out = os.path.join(args.outdir, self._get_out_basename())

    def _get_out_basename(self):
        """
        modify this fun to change out basename format
        """
        args = self.args
        fpart = self._get_fpart()
        return '.'.join(os.path.basename(args.fpath).split('.')[:-1]) + f'{fpart}.hdf5'

    def _get_fpart(self,):
        """
        need modify this function
        """
        return '-example'

    def _set_show_prog(self,):
        args = self.args
        if hasattr(args, 'show_prog'):
            from tqdm import tqdm
            from functools import partialmethod
            disable = not args.show_prog if args.show_prog is not None else None
            tqdm.__init__ = partialmethod(tqdm.__init__, disable=disable)

    def _set_pre_output(self,):
        from .output import set_output
        pre_str = f'[hifast.{os.path.basename(sys.argv[0])[:-3]}]['
        try:
            pre_str += self.project
        except:
            pass
        try:
            pre_str += f"-M{self.nB:02d}"
        except:
            pass
        try:
            pre_str += f"-{self.date}"
        except:
            pass
        pre_str += '] '
        set_output(pre_str)


class BaseIO(Path_IO):
    """
    methods may need be replaced: _get_fpart, _import_m, gen_s2p_out, __call__, _get_out_basename
    """
    # set ver as 'new' or 'old', if old, PolarMjdChan_to_MjdChanPolar() when reading, MjdChanPolar_to_PolarMjdChan() when saving
    ver = 'new'

    def __init__(self, args, dict_in=None, inplace_args=False, HistoryAdd=None):
        """
        args: class
              including attributes: fpath, outdir, frange
        dict_in: if None, load data from args.fpath, if set, omit data in args.fpath
        History_Add: dict, {'key': str}
        """
        self.args = args if inplace_args else copy.deepcopy(args)
        self.dict_in = dict_in
        self.HistoryAdd = HistoryAdd
        self._gen_fpath_out()
        if self.dict_in is None:
            self._check_fout()
        self._set_show_prog()
        self.nB = self.get_nB(self.args.fpath)
        self._import_m()
        self.open_fpath()
        self.load_specs()
        self.load_radec()
        self.load_and_add_Header()
        self._set_pre_output()

    def _import_m(self,):
        """
        need modify this function
        """
        global h5py, OrderedDict, np
        import h5py
        import numpy as np
        from collections import OrderedDict

    def open_fpath(self,):
        if self.dict_in is None:
            self.fin = h5py.File(self.args.fpath, 'r')
            self.fs = self.fin['S']
        else:
            self.fs = self.dict_in

    def close_fpath(self,):
        if self.dict_in is None:
            try:
                self.fin.close()
            except Exception as Err:
                print(Err)

    def load_radec(self,):
        args = self.args
        fs = self.fs
        if 'ra' in fs.keys() and 'dec' in fs.keys():
            self.ra = fs['ra'][:]
            self.dec = fs['dec'][:]
        elif not getattr(args, 'no_radec', False):
            from ..add_radec import get_radec
            nB_radec = getattr(args, 'nB_radec', 1)
            self.ra, self.dec, self.is_extrapo = get_radec(args.fpath, self.mjd, nB_radec)
        else:
            pass

    def load_specs(self,):
        """
        read mjd, freq, spectra(T or flux ) from self.fs
        load and add header
        """
        args = self.args
        fs = self.fs

        self.mjd = fs['mjd'][:]
        self.freq = fs['freq'][:]
        self.freq_ori = self.freq
        if getattr(args, 'frange', False) and (args.frange[0] > 0. or args.frange[1] < float('inf')):
            self.is_use_freq = (self.freq >= args.frange[0]) & (self.freq <= args.frange[1])
            self.freq = self.freq[self.is_use_freq]
        else:
            self.is_use_freq = None

        if 'T' in fs.keys():
            s2p = fs['T']
            infield = 'T'
            outfield = 'Ta'
        elif 'Ta' in fs.keys():
            s2p = fs['Ta']
            infield = 'Ta'
            outfield = 'Ta'
        elif 'Power' in fs.keys():
            s2p = fs['Power']
            infield = 'Power'
            outfield = 'Power'
        elif 'flux' in fs.keys():
            s2p = fs['flux']
            infield = 'flux'
            outfield = 'flux'
            if getattr(args, 'flux', False):
                raise(ValueError(f'flux already exists, please remove \"--flux\" or set as False'))
        else:
            raise(ValueError('can not find spec'))
        if getattr(args, 'flux', False):
            outfield = 'flux'
        #
        if self.is_use_freq is not None:
            inds = np.where(self.is_use_freq)[0]
            if len(inds) == 0:
                raise(ValueError('please check --frange'))
            s2p = s2p[..., inds[0]:inds[-1]+1]  # freq axis is at end; inds is continuous
        if self.ver == 'old':
            s2p = PolarMjdChan_to_MjdChanPolar(s2p[:])
        self.s2p = s2p
        self.outfield = outfield
        self.infield = infield

    def load_and_add_Header(self,):
        Header = OrderedDict()
        if self.dict_in is None:
            fin = self.fin
            if 'Header' in fin.keys():
                try:
                    Header.update(OrderedDict(fin['Header'].attrs.items()))
                except:
                    print('input file have no Header')
        else:
            if 'Header' in self.dict_in.keys():
                Header.update(self.dict_in['Header'])
        import json
        if hasattr(self, 'HistoryAdd') and self.HistoryAdd is not None:
            his = rec_his(args=json.dumps(self.args.__dict__), **self.HistoryAdd)
        else:
            his = rec_his(args=json.dumps(self.args.__dict__))
        Header.update(his)
        self.Header = Header

    def gen_dict_out(self, *args, **kwargs):
        """
        add self.ra, self.dec, self.mjd, self.freq
        add MjdChanPolar_to_PolarMjdChan(self.s2p_out) to dict_out[outfield]
        add added field:
            self.add_fields = ['is_on', 'next_to_cal', 'is_delay', 'Tcal', 'is_extrapo', 'vel']
        add args in fpath file
        add key:value in kwargs
        add Header
        """
        dict_out = {}
        dict_out['mjd'] = self.mjd
        dict_out['freq'] = self.freq
        if self.ver == 'old':
            self.s2p_out = MjdChanPolar_to_PolarMjdChan(self.s2p_out)
        dict_out[self.outfield] = self.s2p_out
        # add field in add_fields and args from self.fs
        if not hasattr(self, 'add_fields'):
            self.add_fields = ['is_on', 'next_to_cal', 'is_delay', 'Tcal', 'is_extrapo', 'vel', 'is_rfi']
            self.add_fields += ['is_excluded']
        self.add_fields += args
        for field in self.add_fields:
            if field in self.fs.keys():
                dict_out[field] = self.fs[field][()]
        # process if set frange
        if self.is_use_freq is not None:
            if 'vel' in dict_out.keys():
                dict_out['vel'] = dict_out['vel'][self.is_use_freq]
            if 'Tcal' in dict_out.keys():
                # use ``try`` for backwards compatible
                try:
                    dict_out['Tcal'] = dict_out['Tcal'][:, self.is_use_freq]
                except:
                    pass
            for key in ['is_excluded', 'is_rfi']:
                if key in dict_out.keys():
                    dict_out[key] = dict_out[key][:, self.is_use_freq]
        # add ra dec
        for key in ['ra', 'dec', 'is_extrapo']:
            if hasattr(self, key):
                val_tmp = getattr(self, key)
                if val_tmp is not None:
                    dict_out[key] = val_tmp
        # add or replace from kwargs
        for key in kwargs.keys():
            dict_out[key] = kwargs[key]
        dict_out['Header'] = self.Header
        self.dict_out = dict_out

    def save(self,):
        print("Saving...")
        save_specs_hdf5(self.fpath_out, self.dict_out, wcs_data_name=self.outfield)
        print(f"Saved to {self.fpath_out}")

    def gen_s2p_out(self):
        """
        need modify this function
        """
        # process spectra
        self.s2p_out = self.s2p

    def __call__(self, save=True):
        # store output data in self.dict_out
        self.gen_s2p_out()
        self.gen_dict_out()
        # save to hdf5 file
        if save:
            self.save()


class HFDataT:
    """
    h5py dataset with Ta(flux, DATA...) transpose
    """
    def __init__(self, dataset):
        assert dataset.ndim == 3
        self.dataset = dataset
        self.ndim = dataset.ndim
        self.size = dataset.size
        self.dtype = dataset.dtype

    @property
    def shape(self):
        shape = self.dataset.shape
        return (shape[1], shape[2], shape[0])

    def __len__(self,):
        return len(self.dataset)

    def __repr__(self,):
        return f"shape {self.shape}, type {self.dtype}"

    def __getitem__(self, arg_in):
        """"""
        arg = [slice(None, None, None),]*3
        if isinstance(arg_in, int) or isinstance(arg_in, slice):
            arg[0] = arg_in
        elif isinstance(arg_in, type(Ellipsis)):
            pass
        elif isinstance(arg_in, tuple):
            ell_ind = None
            for i, a in enumerate(arg_in):
                if isinstance(a, type(Ellipsis)):
                    ell_ind = i
                    break
            if ell_ind is not None:
                ind = ell_ind
                arg[:ind] = arg_in[:ind]
                arg[ind-len(arg_in):] = arg_in[ind-len(arg_in):]
            else:
                arg[:len(arg_in)] = arg_in
        else:
            raise(ValueError(f'slic err {arg_in}'))


        arg_T = [slice(None, None, None),]*3
        for i in range(len(arg)):
            if isinstance(arg[i], int):
                arg[i] = slice(arg[i], arg[i]+1)
                arg_T[i] = 0
        arg = (arg[2], arg[0], arg[1])
        arg_T = tuple(arg_T)

        res = self.dataset.__getitem__(arg)
        res = res.transpose((1, 2, 0))[arg_T]

        return res

class HFGroup:

    def __init__(self, group, kT, DATA_key=None):
        """
        group: input h5py group
        kT: key need transpose
        """
        self.group = group
        self._key_trans_ = kT
        self.keys = group.keys
        self.DATA_key = DATA_key

    def __getitem__(self, arg):
        if arg == 'DATA' and 'DATA' not in self.group.keys():
            if self.DATA_key is None:
                raise(ValueError('no DATA in HFGroup'))
            else:
                arg = self.DATA_key
        r = self.group.__getitem__(arg)
        if arg == self._key_trans_:
            r = HFDataT(r)
        return r

    def __repr__(self,):
        return f"keys: {self.group.keys()}"


from dataclasses import dataclass
@dataclass
class H5HDU:
    name: str
    data: HFGroup
    header: dict

class H5FitsRead:

    def __init__(self, fpath, kT='DATA', hdu_names=['0', '1'], DATA_key=None):
        import h5py
        f = h5py.File(fpath, 'r')
        self.f = f
        self.close = f.close

        self.hdus = {}
        for name in hdu_names:
            try:
                if name.isdigit():
                    self.hdus[int(name)] = self._read_hdu(name, kT, DATA_key)
                else:
                    self.hdus[name] = self._read_hdu(name, kT, DATA_key)
            except exception as e:
                print(e)
                pass

    def _read_hdu(self, hdu_name, kT, DATA_key=None):
        return H5HDU(hdu_name,
                     HFGroup(self.f[hdu_name], kT, DATA_key),
                     dict(self.f[hdu_name].attrs.items())
                     )

    def __getitem__(self, arg):
        return self.hdus[arg]
