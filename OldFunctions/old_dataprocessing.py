import numpy as np
import scipy.io as sio
import yaml
import matplotlib.pyplot as plt
import os


def read_2d_data(fname):

    """ Reads a 2D data set from a given file. Returns a data
    dictionary with keys:

    1. excitation
    2. emission
    3. spectrum
    4. temperature
    5. acquistion_time

    Are there other parameters we should care about? """

    rawdata = sio.loadmat(fname)
    # the emission axis is this in most cases
    # if it changes I'll have to redo this file
    emission = (np.linspace(1, 1024, 1024) -
                -1.46581267e+04)/9.91909400e+00 + 20.
    spectrum = rawdata['Spectrum2D']
    excitation = rawdata['ActualWavelength'][0]

    with open(fname[:-3] + 'yml', 'r') as stream:
        try:
            metadata = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    datadict = {
        "emission": emission,
        "excitation": excitation,
        "spectrum": spectrum,
        "temp": metadata['temp'],
        "acquistion_time": metadata['acquisition_time']
    }

    return datadict


def correct_for_acqusition(data):

    """ Normalizes for acqusition time. Returns spectrum in units
    of counts per second. """

    data["spectrum"] = (data["spectrum"] - 700)/data["acqusition_time"]

    return data










fstring = "PTO_DSO/PTO_DSO,EDFA=867mA,gain=110dB2Dsweep_%d.mat"
nfiles = 4

fnames = [fstring % i for i in range(nfiles)]
data = [read_2d_data[fnames[i]] for i in range(nfiles)]
