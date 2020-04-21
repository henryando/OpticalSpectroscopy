import numpy as np
from scipy.io import loadmat
import yaml
import constants as const
import glob
from spectrumtools import Spectrum
from spectrumtools import Excitation
import conversions as conv


class ScanTypeError(Exception):
    """An error type to describe when you try to process a scan of one
    type as a scan of a different type by accident.
    """

    def __init__(self, relative_filepath):
        self.relative_filepath = relative_filepath


def get_filepaths(relative_filepath):
    """Returns all the .mat filepaths for a given folder."""
    return glob.glob("%s/*.mat" % relative_filepath)


def read_metadata(relative_filepath):
    """Scrapes the metadata we care about: temperatures, acquisition times,
    lockin sensitivity.
    """
    with open(relative_filepath[:-3] + "yml", "r") as stream:
        try:
            metadict = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    return (metadict["temp"], metadict["acquisition_time"], metadict["sens"])


def read_excitation(relative_filepath):
    """Returns a dictionary with:

    'ex' (excitation wavelengths),
    'sig' (the lockin signal, rotated for maximum signal strength),
    'temp' (the temperature),
    'sens' (the lockin sensitivity).
    """
    rawdata = loadmat(relative_filepath)
    temp, _, sens = read_metadata(relative_filepath)

    try:
        # Rotate lockin phase to get maximum signal strength relative to noise
        theta = np.arctan(np.mean(rawdata["lockinI"]) / np.mean(rawdata["lockinQ"]))
        sig = np.cos(theta) * rawdata["lockinI"] + np.sin(theta) * rawdata["lockinQ"]
        # Is this the right formula?

        return Excitation(
            conv.nm_to_wavenumber(rawdata["NominalWavelength"][0]),
            sig[0],
            float(temp),
            sens,
        )
    except KeyError:
        raise ScanTypeError(relative_filepath)


def read_all_excitation(folder):
    """Returns a list of data dictionaries for all the 2d spectra in a folder."""
    fnames = get_filepaths(folder)
    spectra = []
    for f in fnames:
        try:
            spectra.append(read_excitation(f))
        except ScanTypeError:
            pass
    return spectra


def read_2dspectrum(relative_filepath):
    """Returns a dictionary with:

    'ex' (excitation wavelengths),
    'em' (emission wavelengths),
    'spec' (the spectrum, in normalized counts per second),
    'temp' (the temperature),
    'time' (the acquisition time in seconds).
    """
    rawdata = loadmat(relative_filepath)

    temp, time, _ = read_metadata(relative_filepath)
    try:
        ex = rawdata["ActualWavelength"][0]
        # Convert signal to counts per second
        spec = (rawdata["Spectrum2D"] - const.READOUT_NOISE_COUNTS) / time
        spec = spec[:, ex > 0]
        ex = conv.nm_to_wavenumber(ex[ex > 0])
        # print(temp)
        # print(relative_filepath)
        if len(ex) == 0:
            raise ScanTypeError(relative_filepath)
        return Spectrum(
            ex,
            conv.nm_to_wavenumber(const.CALIBRATED_EMISSION),
            spec,
            float(temp),
            time,
        )
    except KeyError:
        raise ScanTypeError(relative_filepath)


def read_all_2dspectra(folder):
    """Returns a list of data dictionaries for all the 2d spectra in a folder."""
    fnames = get_filepaths(folder)
    spectra = []
    for f in fnames:
        try:
            spectra.append(read_2dspectrum(f))
        except ScanTypeError:
            pass
    return spectra


# fp = "Data/Fa/Fa_PbWO4,EDFA=867mA,gain=110dB,ExcitationScan,spectrometer=0.0nm,
# 1490.0-1560.0nm,scan_rate=0.1nm_s.mat"
