
def check_backend():
    import matplotlib as mpl
    # check if it is already ipympl
    if 'ipympl' in mpl.get_backend():
        return

    # try to switch
    try:
        from IPython import get_ipython
        ipython = get_ipython()
        if ipython:
            # equivalent to %matplotlib ipympl
            ipython.run_line_magic('matplotlib', 'ipympl')
            print("Switched to 'ipympl' backend automatically.")
    except Exception as e:
        # print(f"Failed to switch backend: {e}")
        pass

    # check again
    if 'ipympl' not in mpl.get_backend():
         raise RuntimeError("Please run `%matplotlib ipympl` manually.")
