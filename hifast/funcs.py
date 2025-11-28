

__all__ = []


import h5py
import warnings


from .core.cal import CalOnOff
from .core.cal import FastRawSpec
from .core.radec import get_radec

from .core.flux import FluxCali

from .core.corr_vel import frame_correct_freq
from .core.corr_vel import freq2vel, vel2freq

from .utils.io import MjdChanPolar_to_PolarMjdChan
from .utils.io import PolarMjdChan_to_MjdChanPolar
from .utils.io import get_nB
from .utils.io import replace_nB
from .utils.io import load_hdf5_to_dict

try:
    from .interaction import bld_i as interact_bld
    from .interaction import sw_i as interact_sw
except ImportError as err:
    print(err)
    warnings.warn(f"{err}. Failed to import interactive modules")

from .utils.io import HFDataT

from .utils import obs_log


class HFSpec:

    def __init__(self, s, *,
                 kT=('Ta', 'flux', 'Power', 'DATA'),
                 DATA_cand=('Ta', 'flux', 'Power')):
        """
        s: str
        kT: field needed to be transposed by PolarMjdChan_to_MjdChanPolar
        DATA_cand: key candidates as 'DATA'
        """
        if isinstance(s, str):
            self.fin = h5py.File(s, 'r')
        else:
            raise(ValueError('input'))
        self._group = self.fin['S']
        self._keys_trans_ = kT
        self._DATA_cand = DATA_cand

        self._gen_DATA()
        self._set_attrs()


    def _gen_DATA(self,):
        self._DATA = None
        keys_ = self._group.keys()
        if 'DATA' in keys_:
            self._DATA = 'DATA'
        else:
            for key in self._DATA_cand:
                if key in keys_:
                    self._DATA = key
                    break
        # if self._DATA is None:
        #     raise(KeyError('can not find spec data'))

    def keys(self,):
        keys = list(self._group.keys())
        if self._DATA is not None and 'DATA' not in keys:
            keys.append('DATA')
        return tuple(keys)

    @property
    def Header(self,):
        Header = dict(self.fin['Header'].attrs.items())
        return Header

    def __getitem__(self, arg):
        if arg == 'DATA':
            arg = self._DATA
        r = self._group.__getitem__(arg)
        if arg in self._keys_trans_:
            r = HFDataT(r)
        if r.ndim == 1:
            r = r[:]
        return r

    def __repr__(self,):
        return f"keys:\n{repr(self.keys())}"

    def _set_attrs(self,):
        keys_forbid = ['keys', 'fin']
        for key in self.keys():
            if key in keys_forbid:
                warnings.warn(f'can not set {key} as attrs')
            else:
                setattr(self, key, self.__getitem__(key))
