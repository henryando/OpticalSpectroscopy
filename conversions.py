import numpy as np
import constants as const


def nm_to_wavenumber(wavelength):
    """ Convert wavelength from nm to wavenumbers (cm^-1). """
    return 1e7 / wavelength
