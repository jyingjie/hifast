

__all__ = ['ConvFun']


import numpy as np
from scipy import special
from functools import partial


class ConvFun(object):
    def __init__(self, kernel_type, beam_fwhw, *,
                gaussian_fwhw=None,
                bsize=None, gsize=None,
                wei_min=None):
        """
        kernel_type: str
            'gaussian':
            'bessel_gaussian':
            'sinc_gaussian':

        """
        if kernel_type == 'gaussian':
            if gaussian_fwhw is None:
                gaussian_fwhw = beam_fwhw/2
            self.sigma = gaussian_fwhw / (2*(2*np.log(2))**0.5)
            self.dis_cut_suggest = 3*self.sigma
            self.gaussian_fwhw = gaussian_fwhw

        elif kernel_type == 'bessel_gaussian' or kernel_type == 'sinc_gaussian':
            self.bsize = 1.55 * beam_fwhw / 3. if bsize is None else bsize
            self.gsize = 2.52 * beam_fwhw / 3. if gsize is None else gsize
            if kernel_type == 'bessel_gaussian':
                # first zero special.jn_zeros(1, 1) = 3.8317059702075
                self.dis_v0 = 3.8317059702075*self.bsize/np.pi
            elif kernel_type == 'sinc_gaussian':
                self.dis_v0 = self.bsize
            self.dis_cut_suggest = self.dis_v0
        else:
            raise(ValueError(f'kernel \"{kernel_type}\" is not supported'))

        self.kernel_type = kernel_type
        self.wei_min = wei_min
        self.beam_fwhw = beam_fwhw

        self.beam_change()

    def beam_change(self):
        # calculate new beam fwhw
        if self.kernel_type == 'gaussian':
            self.new_beam_fwhw = (self.beam_fwhw**2 + self.gaussian_fwhw**2)**0.5
            # volume under the Gaussian function: 2*pi*A*sigma_X*sigma_Y
            self.beam_factor = self.new_beam_fwhw**2 / self.beam_fwhw**2

        else:
            # To do: use Numerical method to calculate
            pass

    def __call__(self, dis):

        wei = getattr(self, f"{self.kernel_type}_kernel")(dis)

        if self.wei_min is not None:
            wei[wei < self.wei_min] = 0

        return wei

    def gaussian_kernel(self, dis):
        sigma = self.sigma
        w = np.exp(-dis**2/(2*sigma**2))
        return w

    def bessel_gaussian_kernel(self, dis):
        """
        """
        a, b = self.bsize, self.gsize
        x = np.pi*dis/a
        wei = np.empty_like(dis)
        not_zero = dis != 0
        wei[~not_zero] = 1 # peak value depend on the not_zero calculating method
        wei[not_zero] = 2*special.j1(x[not_zero])/(x[not_zero])*np.exp(-(dis[not_zero]/b)**2)

        wei[abs(dis) > self.dis_v0] = 0
        return wei

    def sinc_gaussian_kernel(self, dis):

        a, b = self.bsize, self.gsize
        x = np.pi*dis/a
        wei = np.empty_like(dis)
        not_zero = dis != 0
        wei[~not_zero] = 1
        wei[not_zero] = np.sin(x[not_zero])/(x[not_zero])*np.exp(-(dis[not_zero]/b)**2)

        wei[abs(dis) > self.dis_v0] = 0
        return wei
