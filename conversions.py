import numpy as np
import constants as const


def nm_to_wavenumber(wavelength):
    """ Convert wavelength from nm to wavenumbers (cm^-1). """
    return 1e7 / wavelength


def linewidth_to_nsamples(wavelengths, linewidth):
    """Given an array of wavelengths like emission or excitation, return the
    number of samples that corresponds to the appropriate linewidth.
    """
    return int(
        np.ceil(linewidth * wavelengths.size / (wavelengths[-1] - wavelengths[0]))
    )
