import numpy as np
import spectrumtools as st
import readdata as rd
import peakfindinggui as pfg
import peakgroupinggui as pgg
import matplotlib.pyplot as plt


def read_all_2dspectra(folder):
    """Return a list of Spectrum objects for all the 2d spectra in a given folder.
    Also print the number of spectra and the mean temperature of the spectra.
    Each object has fields: 'em': emission wavelengths, 'ex': excitation
    wavelengths, 'spec': spectrum2d data, 'temp': temperature, 'time': acquisition
    time for spectrometer.
    """
    spectra = rd.read_all_2dspectra(folder)
    print("There are %d spectra in this folder." % len(spectra))
    print("The temperature is %.1f." % spectra[0].temp)
    return spectra


def iterative_smooth(spectra, N=2, Wx=2, Wy=2):
    """Takes a list of spectrum objects (or a single object). Smoothes
    the data with a simple moving window average. Smoothing alternates x and y,
    a total of N times. The windows in x and y can be chosen independently.
    """
    return st.iterative_smooth(spectra, N=N, Wx=Wx, Wy=Wy)


def find_peaks(spectra, linewidth=0.07, noisefraction=1 / 10):
    """Locates potential peaks based on the minimum permissible peak height
    (background noise time noisefraction) and the minimum permissible permissible
    linewidth over which the peak must be maximal. Then opens a GUI to select the
    potential peaks which we actually want. Returns a tuple with two entries,
    excitation coordinates of peaks and emission coordinates of peaks.
    """
    peaks = pfg.find_peaks(spectra, linewidth=linewidth, nf=noisefraction)
    return (peaks.ex, peaks.em)


def find_lines(spectra, peaks, Wx=0.3, Wy=1.5):
    """Takes a list of spectra objects and peak objects. Group the peaks
    into lines by hand using a GUI. The rows and columns of the GUI are set by the
    linewidths Wx and Wy. Returns the lines in excitation and emission.
    """
    return pgg.find_lines(spectra, peaks, grouping_Wx=Wx, grouping_Wy=Wy)


def plot_spectra(spectra, clevels=None, figsize=(7, 7), figure=None):
    """Plots a 2d spectrum (or list of spectra). Returns fig and ax handles."""
    return pfg.plot2d(spectra, clevels=clevels, figsize=figsize, inputfigure=figure)


def plot_lines(lines, fmt="g", ax=None):
    """Plots the emission and excitation lines."""
    if ax is None:
        ax = plt.gca()
    exlines = lines[0]
    emlines = lines[1]
    xlim = plt.xlim()
    ylim = plt.ylim()
    ax.plot([exlines, exlines], ylim, fmt, linewidth=0.5)
    ax.plot(xlim, [emlines, emlines], fmt, linewidth=0.5)


def plot_peaks(peaks, color="r", size=2, ax=None):
    """Takes a peak dictionary. Plots the peaks with small red dots."""
    if type(peaks) is st.Peaks:
        (x, y) = (peaks.ex, peaks.em)
    else:
        (x, y) = (peaks[0], peaks[1])
    plt.plot(
        x,
        y,
        "o",
        markersize=size,
        markeredgecolor="w",
        markeredgewidth=0.5,
        markerfacecolor="k",
    )


# def _em_at_given_ex(data, ex):
#     if type(data) is list:
#         for d in data:
#             result = _em_at_given_ex(d, ex)
#             if result is not None:
#                 return result
#     else:
#         if (min(data.ex) < ex) and (ex < max(data.ex)):
#             return data.spec[:, sum(data.ex < ex) + 1]


# def peak_heights_at_lines(data, exlines, emlines):
#     em = data[0].em
#     peaks = np.zeros(len(exlines), len(emlines))


def ex_line_nearest_peak(data, ex):
    """Takes a spectrum and an excitation line, and returns the
    emission line of the local max closest in excitation to the given
    line.
    """
    if type(data) is list:
        for s in data:
            emspec = ex_line_nearest_peak(s, ex)
            if emspec is not None:
                # print("returning em")
                return emspec
    else:
        if (ex < min(data.ex)) or (ex > max(data.ex)):
            # print("out of bounds")
            return None
        else:
            # print("in bounds")
            ind = sum(data.ex >= ex)
            return data.spec[:, ind + 1]


def peak_heights_at_lines(data, exlines, emlines):
    """Takes a spectrum and a set of lines, and returns an array with the spectral
    strength at the intersection of all the lines.
    """
    peaks = np.zeros((len(exlines), len(emlines)))
    em = data[0].em
    for i in range(len(exlines)):
        curr = ex_line_nearest_peak(data, exlines[i])
        if curr is None:
            peaks[i] = np.nan
        else:
            for j in range(len(emlines)):
                peaks[i, j] = curr[sum(em > emlines[j]) + 1]

    # peaks = peaks - np.amin(peaks)
    # peaks = peaks / peaks[0, 0]
    return peaks
