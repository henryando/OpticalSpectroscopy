import spectrumtools as st
import readdata as rd
import peakfindinggui as pfg
import peakgroupinggui as pgg


def read_all_2dspectra(folder):
    """Return a list of dictionaries for all the 2d spectra in a given folder.
    Also print the number of spectra and the mean temperature of the spectra.
    Each dictionary has entries: 'em': emission wavelengths, 'ex': excitation
    wavelengths, 'spec': spectrum2d data, 'temp': temperature, 'time': acquisition
    time for spectrometer.
    """
    spectra = rd.read_all_2dspectra(folder)
    print("There are %d spectra in this folder." % len(spectra))
    print("The temperature is %.1f." % spectra[0]["temp"])
    return spectra


def iterative_smooth(spectra, N=2, Wx=2, Wy=2):
    """Takes a list of spectrum dictionaries (or a single dictionary). Smoothes
    the data with a simple moving window average. Smoothing alternates x and y,
    a total of N times. The windows in x and y can be chosen independently.
    """
    return st.iterative_smooth(spectra, N=N, Wx=Wx, Wy=Wy)


def find_peaks(spectra, linewidth=0.07, noisefraction=1 / 10):
    """Locates potential peaks based on the minimum permissible peak height
    (background noise time noisefraction) and the minimum permissible permissible
    linewidth over which the peak must be maximal. Then opens a GUI to select the
    potential peaks which we actually want. Returns a dict with two entries,
    'ex': excitation coordinates of peaks and 'em': emission coordinates of peaks.
    """
    goodpeaks = pfg.find_peaks(spectra, linewidth=linewidth, nf=noisefraction)
    goodpeaks["ex"] = goodpeaks["ex"][0]  # Not really sure why I have to do this.
    goodpeaks["em"] = goodpeaks["em"][0]
    return goodpeaks


def find_lines(spectra, peaks, Wx=0.3, Wy=1.5):
    """Takes a list of spectra dictionaries and peak dictionaries. Group the peaks
    into lines by hand using a GUI. The rows and columns of the GUI are set by the
    linewidths Wx and Wy. Returns the lines in excitation and emission.
    """
    return pgg.find_lines(spectra, peaks, grouping_Wx=Wx, grouping_Wy=Wy)


def plot_spectra(spectra, clevels=None, figsize=(7, 7), figure=None):
    """Plots a 2d spectrum (or list of spectra). Returns fig and ax handles."""
    return pfg.plot2d(spectra, clevels=clevels, figsize=figsize, inputfigure=figure)


def plot_lines(lines):
    """Takes a tuple of (excitation lines, emission lines). Plots them to the axis.
    """
    pfg.plot_lines(lines)


def plot_peaks(peaks):
    """Takes a peak dictionary. Plots the peaks with small red dots."""
    pfg.plot_peaks(peaks["ex"], peaks["em"])


def plot_good_peaks(peaks):
    """Takes a peak dictionary. Plots the peaks with big green dots."""
    pfg.plot_good_peaks(peaks["ex"], peaks["em"])
