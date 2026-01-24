
def is_notebook():
    """
    Check if the code is running in a Jupyter Notebook or IPython kernel.
    Returns True if interactive (Notebook/Jupyter Lab/QtConsole/IPython Terminal).
    Returns False if running as a standard script.
    """
    try:
        shell = get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell':
            return True   # Jupyter notebook or qtconsole
        elif shell == 'TerminalInteractiveShell':
            return True   # Terminal running IPython
        else:
            return False  # Other type (?)
    except NameError:
        return False      # Probably standard Python interpreter
