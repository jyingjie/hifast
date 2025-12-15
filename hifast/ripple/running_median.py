# author: Xu Chen, Li Fujia, 2021.06
# code: Xu Chen

import numpy as np
from astropy import log
from copy import deepcopy
from tqdm import tqdm
# from matplotlib import pyplot as plt

from multiprocessing import Process, Queue
from threadpoolctl import ThreadpoolController
controller = ThreadpoolController()

def check_bottleneck():
    try:
        import bottleneck as bn
    except ModuleNotFoundError:
        bn = np
        log.warning("Install the package 'bottleneck' will help to speed up. Now use 'numpy' instead.")
    return bn

bn = check_bottleneck()


def _running_block(q, data, n, func, method,):
    
    if func == 'iter':
#         raise ValueError(np.sum(np.isnan(data)))
        
        tlen = data.shape[0]

        t1 = np.arange(tlen)-n
        t2 = np.arange(tlen)+n
        t1[t1 < 0] = 0
        t2[t2 < 2*n] = 2*n
        t1[t1 > tlen - 2*n - 1] = tlen - 2*n - 1
        t2[t2 > tlen - 1] = tlen - 1

        sw = np.zeros_like(data)
        for i in tqdm(range(tlen)):
            if method == 'median':
                sw[i] = bn.nanmedian(data[t1[i]:t2[i]], axis = 0)
            elif method == 'mean':
                sw[i] = bn.nanmean(data[t1[i]:t2[i]], axis = 0)
            
    elif func == 'smooth':
        from ..utils.misc import smooth1d
        if method == 'median':
            s_method_t = 'median'
        elif method == 'mean':
            s_method_t = 'boxcar' 
        print('Smooth ing ...')
        sw = smooth1d(data, axis=0, sigma=n, method=s_method_t)
    
    q.put(sw)

@controller.wrap(limits=1, user_api='blas')
def running_median_mp(data, n, func, method, nproc):

    n_proc = min(nproc, data.shape[1]) # freq

    ps = []
    qs = []
    i_list = np.arange(data.shape[1])
    d_list = np.array_split(i_list, n_proc, axis = 0)
    for d in d_list:
        q = Queue()
        p = Process(target = _running_block, 
                    args=(q, data[:, d, :], n, func, method,))

        ps += [p]
        qs += [q]
    for p in ps:
        p.start()
    res = np.hstack([q.get()for q in qs])

    for p in ps:
        p.join()  # need after q.get()

    return res


def running_median(data, nspec, func='iter', method = 'median', nproc = 1):
    """
    data: (mjd, freq, polar)
    nspec: moving window length
    func: 'iter' using numpy or bottleneck
          'smooth' boxcar filter
    """
    if method not in ['median', 'mean']:
        raise ValueError(f"{method} should be 'median' or 'mean'!")
    
    n = nspec // 2
    
    if func == 'iter':
        tlen = data.shape[0]
        if tlen < 2*n:
            print(f"use all specs to {method}")
        else:
            print(f"use {2*n} specs to {method}")
    elif func == 'smooth':
        print(f"use {2*n +1} specs to {method}")
        
    sw = running_median_mp(data, n, func, method, nproc)

    return sw


### minmed

def get_trange(N, npart = None, dpart = None,):
    if dpart is None:
        dpart = N // npart
    p1 = np.arange(0, N, dpart)
    p2 = np.arange(dpart, N + dpart, dpart)
    p2[p2 > N] = N
    
    if (N - p1[-1]) < 0.75 * dpart and (N - p1[-1]) > 0:
        p1 = np.delete(p1, -1)
        p2 = np.delete(p2, -2)
    return p1, p2


def minmed(data, nsection = None, nspec = None, npart = 1, method = 'MedMed'):
    """
    data: 2D
    npart: divide the data into n parts
    nsection: divide each part into n sections
    then calculate the median of each section.
    use the 'median' or 'min' as the baseline of this part.
    Ref: Putman et al. 2002 
    https://ui.adsabs.harvard.edu/link_gateway/2002AJ....123..873P/doi:10.1086/338088
    """
    bn = check_bottleneck()
    
    N = data.shape[0]
    print(f"Divided the data into {npart} part(s).")
    if nsection is None: 
        print(f"Each part has {N//npart//nspec} sections. Each section has {nspec} specs.")
    elif nspec is None:
        print(f"Each part has {nsection} sections. Each section has {N//npart//nsection} specs.")
    
    p1s, p2s = get_trange(N, npart)
    
    bsl = np.zeros_like(data)
    for p1, p2 in zip(p1s, p2s):
        print("part:",[p1, p2])
        s1, s2 = get_trange(p2 - p1, nsection, nspec)
        ind1, ind2 = s1 + p1, s2 + p1

        meds = np.zeros(np.hstack((len(ind1),data.shape[1:])))
        for i in tqdm(range(len(ind1))):
            meds[i] = bn.nanmedian(data[ind1[i]:ind2[i]], axis = 0)
        
        if method == 'MinMed':
            med = bn.nanmin(meds, axis = 0)
        elif method == 'MedMed':
            med = bn.nanmedian(meds, axis = 0)
            
        bsl[p1:p2] = med
        
    return bsl
