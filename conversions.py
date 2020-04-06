import numpy as np
import constants as const


def nm_to_wavenumber(data):
    """ Convert wavelength from nm to wavenumbers (cm^-1). """
    if type(data) is list:
        return [nm_to_wavenumber(s) for s in data]
    else:
        return 1e7 / data


def linewidth_to_nsamples(wavelengths, linewidth):
    """Given an array of wavelengths like emission or excitation, return the
    number of samples that corresponds to the appropriate linewidth.
    """
    return int(
        np.ceil(linewidth * wavelengths.size / (wavelengths[-1] - wavelengths[0]))
    )
