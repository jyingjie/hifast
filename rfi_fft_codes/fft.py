#!/usr/bin/env python
# coding: utf-8

# author: Xu Chen, Li Fujia, 2021.06
# code: Xu Chen

import numpy as np
import numpy.fft as fft


class FFT(object):
    def __init__(self, data,freq):
        # 采样点数
        N = freq.shape[0]
        # 采样周期T
        fdelta = freq[1] - freq[0]
        # 采样频率
        #fs = 1/fdelta
        
        # fft
        if len(data.shape) == 1:
            A_data = fft.rfft(data)
        else:
            A_data = fft.rfft(data,axis = 1)
            
        self.complex_num = A_data
        self.x = fft.rfftfreq(N,fdelta)
        self.amp = np.abs(A_data)
        self.phi = np.angle(A_data)

        