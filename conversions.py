import numpy as np
import constants as const


def nm_to_wavenumber(data):
    """ Convert wavelength from nm to wavenumbers (cm^-1). """
    if type(data) is list:
        return [nm_to_wavenumber(s) for s in data]
    else:
        return 1e7 / data


def kelvin_to_wavenumber(temp):
    """Convert kelvin to wavenumbers."""
    return temp * 0.695028
