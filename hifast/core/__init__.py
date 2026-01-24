from ..utils import is_notebook
if not is_notebook():
    import matplotlib
    matplotlib.use('Agg')
    from ..utils import set_err
