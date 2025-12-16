import __main__
if hasattr(__main__, '__file__'):
    import matplotlib
    matplotlib.use('Agg')
    from ..utils import set_err
