

__all__ = ['flag_data']


import aoflagger

def flag_data(data, strategy_file):
    """
    Parameters
    ----------
    data: np.narray
        shape: (n_times, n_channel, n_polar)
    strategy_file: str
        file path of strategy file

    Returns
    -------
    flagvalues: bools
        shape: (n_times, n_channel)
    """
    # https://aoflagger.readthedocs.io/en/latest/python_interface.html
    # https://www.andreoffringa.org/aoflagger/doxygen/classaoflagger_1_1AOFlagger.html
    # overview
    flagger = aoflagger.AOFlagger()
    strategy = flagger.load_strategy_file(strategy_file)

    # make date set
    data_ao = flagger.make_image_set(*data.shape)
    for imgindex in range(data.shape[-1]):
        data_ao.set_image_buffer(imgindex, data[:, :, imgindex].T)
    # run
    flags = strategy.run(data_ao)
    # or run with existingFlags
    # strategy.run(data_ao, flag)
    flagvalues = flags.get_buffer()
    flagvalues = flagvalues.T
    return flagvalues
