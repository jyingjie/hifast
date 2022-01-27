import mpl_interactions.ipyplot as iplt
import ipywidgets as widgets
from ipywidgets import interact

__all__ = ['iplt',
           'widgets',
           'interact',
           'w_conf',
           '_IntSlider',
           '_BoundedIntText',
           '_BoundedFloatText',
           '_IntText',
           '_FloatSlider',
           '_FloatLog10Slider',
           '_RadioButtons',
           '_Dropdown',
           '_FloatRangeSlider',
          ]


w_conf = {}
w_conf['continuous_update'] = False
w_conf['orientation'] = 'horizontal'
w_conf['readout'] = True
w_conf['layout'] = {}


def _IntSlider(readout_format='d', **kwargs):
    """
    kwargs:
    slider = (value=7, min=0, max=10, step=1)
    """
    res = {}
    for key in kwargs.keys():
        res[key] = widgets.IntSlider(*kwargs[key], description=f'{key}',
                                     **w_conf, readout_format=readout_format)
    return res




def _BoundedIntText(readout_format='d', **kwargs):
    """
    kwargs:
    slider = (value=7, min=0, max=10, step=1)
    """
    res = {}
    for key in kwargs.keys():
        res[key] = widgets.BoundedIntText(*kwargs[key], description=f'{key}',
                                          **w_conf, readout_format=readout_format)
    return res

def _BoundedFloatText(readout_format='.4f', **kwargs):
    """
    kwargs:
    slider = (value=7, min=0, max=10, step=1)
    """
    res = {}
    for key in kwargs.keys():
        v = kwargs[key]
        res[key] = widgets.BoundedFloatText(value=v[0], min=v[1], max=v[2], step=v[3], description=f'{key}',
                                          **w_conf, readout_format=readout_format)
    return res

def _FloatRangeSlider(readout_format='.4f', **kwargs):
    """
    kwargs:
    slider = (value=[5, 7.5], min=0, max=10.0, step=0.1)
    """
    res = {}
    for key in kwargs.keys():
        v = kwargs[key]
        res[key] = widgets.FloatRangeSlider(value=v[0], min=v[1], max=v[2], step=v[3], description=f'{key}',
                                          **w_conf, readout_format=readout_format)
    return res

def _IntText(readout_format='d', **kwargs):
    """
    kwargs:
    slider = (value=7, min=0, max=10, step=1)
    """
    res = {}
    for key in kwargs.keys():
        res[key] = widgets.IntText(*kwargs[key], description=f'{key}',
                                   **w_conf, readout_format=readout_format)
    return res


def _FloatSlider(readout_format='.2f', **kwargs):
    """
    kwargs:
    slider = (value=7, min=0, max=10, step=1)
    """
    res = {}
    for key in kwargs.keys():
        v = kwargs[key]
        res[key] = widgets.FloatSlider(value=v[0], min=v[1], max=v[2], step=v[3], description=f'{key}',
                                       **w_conf, readout_format=readout_format)
    return res


def _FloatLog10Slider(readout_format='.2f', **kwargs):
    """
    log10
    kwargs:
    slider = (value, min exponent of base, max exponent of base, exponent step)
    """
    res = {}
    for key in kwargs.keys():
        v = kwargs[key]
        res[key] = widgets.FloatLogSlider(value=v[0], base=10, min=v[1], max=v[2], step=v[3], description=f'{key}',
                                          **w_conf, readout_format=readout_format)
    return res


def _RadioButtons(**kwargs):
    """
    kwargs:
    slider = options
    """
    res = {}
    for key in kwargs.keys():
        v = kwargs[key]
        res[key] = widgets.RadioButtons(options=v, value=v[0], description=f'{key}',
                                        **w_conf,)
    return res
# widgets.RadioButtons(
#     options=['pepperoni', 'pineapple', 'anchovies'],
# #    value='pineapple', # Defaults to 'pineapple'
# #    layout={'width': 'max-content'}, # If the items' names are long


def _Dropdown(**kwargs):
    """
    kwargs:
    slider = options
    """
    res = {}
    for key in kwargs.keys():
        v = kwargs[key]
        res[key] = widgets.Dropdown(options=v, value=v[0], description=f'{key}',
                                    **w_conf,)
    return res
