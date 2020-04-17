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


def wx15_towx13(w, x):
    f4 = 60
    f6 = 13860
    b4 = w * x / f4
    b6 = w * (1 - abs(x)) / f6

    f4 = 60
    f6 = 7560
    c = f4 * b4 / (f6 * b6)
    x = np.sign(x) / (c + 1)
    w = b4 * f4 / x

    return w, x


def wx13_towx15(w, x):
    f4 = 60
    f6 = 7560
    b4 = w * x / f4
    b6 = w * (1 - abs(x)) / f6

    f4 = 60
    f6 = 13860
    c = f4 * b4 / (f6 * b6)
    x = np.sign(x) / (c + 1)
    w = b4 * f4 / x

    return w, x
