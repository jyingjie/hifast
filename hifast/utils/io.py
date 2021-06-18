import os

## io
def rec_his(**kwargs):
    import sys
    import json
    from datetime import datetime
    from collections import OrderedDict
    from ._version import get_versions
    
    history= OrderedDict()
    history['version']= get_versions()['version']
    history['cwd']= os.getcwd()
    history['argv']= ' '.join(sys.argv)
    for key in kwargs.keys():
        history[key]= repr(kwargs[key])
    current_time= datetime.now().strftime("%Y%m%d-%H:%M:%S")
    return OrderedDict({"HISTORY-"+current_time : json.dumps(history,indent=2)})

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

def save_dict_hdf5(fname, dict_in, header=None, mode='w'):
    """
    fname: str
    dict_in: dict; keys of the dict_in are str, value are numpy array like.
    header: dict
    """
    import h5py
    try:
        f= h5py.File(fname,mode)
    except OSError:
        if mode=='w':
            from datetime import datetime
            os.rename(fname,fname+'.bak.empty.'+datetime.now().strftime("%Y%m%d-%H%M%S"))
            f= h5py.File(fname,mode)
    if header is not None:
        # sometimes hdf5 raises error if track_order = True
        f.create_group('Header',track_order=False)
        for key in header.keys():
            f['Header'].attrs[key]= header[key]
    for key in dict_in.keys():
        f[key]= dict_in[key]
    f.close()
