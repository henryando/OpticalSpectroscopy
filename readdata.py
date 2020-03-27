import numpy as np
from scipy.io import loadmat
import yaml


def read_metadata(relative_filepath):
    """ Scrapes the metadata we care about: temperatures, acquisition times,
    lockin sensitivity.
    """
    ymlpath = relative_filepath[:-3] + "yml"
    stream = open(ymlpath, "r")
    metadict = yaml.safe_load(stream)
    return (metadict["temp"], metadict["acquisition_time"], metadict["sens"])


def read_excitation(relative_filepath):
    """ Returns a dictionary with:

    'ex' (excitation wavelengths),
    'sig' (the lockin signal, rotated for maximum signal strength),
    'temp' (the temperature),
    'sens' (the lockin sensitivity).
    """
    rawdata = loadmat(relative_filepath)
    temp, _, sens = read_metadata(relative_filepath)
    sig = np.sqrt(rawdata["lockinI"] ** 2 + rawdata["lockinQ"] ** 2)

    return {"ex": rawdata["ActualWavelength"], "sig": sig, "temp": temp, "sens": sens}


# fp = "Data/Fa/Fa_PbWO4,EDFA=867mA,gain=110dB,ExcitationScan,spectrometer=0.0nm,1490.0-1560.0nm,scan_rate=0.1nm_s.mat"
