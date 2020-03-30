import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import cm
import conversions as conv


def logify_spectrum(spectrum):
    """Converts a spectrum to log scale and fixes error values."""
    spectrum[spectrum <= 0] = np.nan
    return np.log(spectrum)


#######################################################################################
# A key front-end function
#######################################################################################


def iterative_smooth(data, N=2, Wx=2, Wy=2):
    """Smooths a spectrum in both directions N times with a moving
    window of 2W+1 points. Window can be specified for both x and y.
    """
    if type(data) is dict:
        spectrum = logify_spectrum(data["spec"])
    else:
        spectrum = logify_spectrum(data)

    avg = np.mean(np.mean(spectrum))
    sz = spectrum.shape
    target = np.ones((sz[0] + 2 * Wx, sz[1] + 2 * Wy)) * avg
    target[Wx:-Wx, Wy:-Wy] = spectrum

    for i in range(N):
        for j in range(sz[0]):
            target[j + Wx, :] = np.mean(target[j : (j + 1 + 2 * Wx), :], axis=0)
        for j in range(sz[1]):
            target[:, j + Wy] = np.mean(target[:, j : (j + 1 + 2 * Wy)], axis=1)

    if type(data) is dict:
        data["spec"] = np.exp(target[Wx:-Wx, Wy:-Wy])
        return data
    else:
        return np.exp(target[Wx:-Wx, Wy:-Wy])


#######################################################################################


def find_local_maxes_2d(data):
    """Finds the peaks (local max in both directions) of a 2 dimensional data set.
    Use np.nonzero to get indices.
    """
    xpeaks = np.bitwise_and(data[:-2, :] < data[1:-1, :], data[1:-1, :] > data[2:, :])
    ypeaks = np.bitwise_and(data[:, :-2] < data[:, 1:-1], data[:, 1:-1] > data[:, 2:])
    peaks = np.zeros_like(data, dtype=bool)
    peaks[1:-1, 1:-1] = np.bitwise_and(xpeaks[:, 1:-1], ypeaks[1:-1, :])
    return peaks


def filter_peaks_square_perimeter(data, peaks, xrad, yrad, height=0):
    """Given a data array and a set of peaks in that array, filter the peaks
    essentially by whether they are large enough to be plausible (are they greater
    than the mean of a certain square surrounding them?)
    """
    (xpeaks, ypeaks) = np.nonzero(peaks)
    (xmax, ymax) = peaks.shape

    for i in range(xpeaks.size):
        r_left = xpeaks[i] > xrad
        r_right = xpeaks[i] < xmax - xrad
        r_bottom = ypeaks[i] > yrad
        r_top = ypeaks[i] < ymax - yrad
        if r_left and r_right and r_bottom and r_top:
            b = np.mean(data[xpeaks[i] - xrad : xpeaks[i] + xrad, ypeaks[i] - yrad])
            t = np.mean(data[xpeaks[i] - xrad : xpeaks[i] + xrad, ypeaks[i] + yrad])
            le = np.mean(data[xpeaks[i] - xrad, ypeaks[i] - yrad : ypeaks[i] + yrad])
            r = np.mean(data[xpeaks[i] + xrad, ypeaks[i] - yrad : ypeaks[i] + yrad])
            square_mean = (b + t + le + r) / 4
            if square_mean > data[xpeaks[i], ypeaks[i]] - height:
                peaks[xpeaks[i], ypeaks[i]] = 0
        else:
            peaks[xpeaks[i], ypeaks[i]] = 0

    return peaks


def filter_peaks_square_area(data, peaks, xrad, yrad):
    """Given a data array and a set of peaks in that array, filter the peaks based
    on whether they are larger than every other point in the square surrounding
    them.
    """
    (xpeaks, ypeaks) = np.nonzero(peaks)
    (xmax, ymax) = peaks.shape

    for i in range(xpeaks.size):
        r_left = xpeaks[i] > xrad
        r_right = xpeaks[i] < xmax - xrad
        r_bottom = ypeaks[i] > yrad
        r_top = ypeaks[i] < ymax - yrad
        if r_left and r_right and r_bottom and r_top:
            square_max = np.amax(
                data[
                    xpeaks[i] - xrad : xpeaks[i] + xrad,
                    ypeaks[i] - yrad : ypeaks[i] + yrad,
                ]
            )
            if square_max > data[xpeaks[i], ypeaks[i]]:
                peaks[xpeaks[i], ypeaks[i]] = 0
        else:
            peaks[xpeaks[i], ypeaks[i]] = 0

    return peaks


#######################################################################################
# A key front-end function
#######################################################################################


def find_peaks(datadict, linewidth=0.2, noisefraction=(1 / 8)):
    """Takes the emission wavelength vector, excitation wavelength vector,
    spectrum (in the original scale, not log), and returns two vectors, expeaks
    and empeaks, which are the coordinates of the peaks that it found. Optional
    arguments control the linewidth to force the points to be maximal over
    (2*linewidth), and the fraction of the noise level to force the peaks to be
    greater by (noisefraction).
    """
    spectrum = datadict["spec"]
    emission = datadict["em"]
    excitation = datadict["ex"]
    spectrum = np.log(spectrum)
    peaks = find_local_maxes_2d(spectrum)
    emrad = conv.linewidth_to_nsamples(emission, linewidth)
    exrad = conv.linewidth_to_nsamples(excitation, linewidth)
    noise = np.std(spectrum)
    peaks = filter_peaks_square_perimeter(
        spectrum, peaks, exrad, emrad, height=(noise * noisefraction)
    )
    peaks = filter_peaks_square_area(spectrum, peaks, exrad, emrad)
    (xpeaks, ypeaks) = np.nonzero(peaks)
    peakdict = {"ex": excitation[ypeaks], "em": emission[xpeaks]}

    return peakdict


#######################################################################################


def plot_peaks(ax, peakdict):
    """Plots peaks on the given axis. """
    ax.scatter(peakdict["ex"], peakdict["em"], s=20, facecolors="none", edgecolors="r")


def manual_peak_sort(ax, peakdict):
    """Lets the user filter the peaks manually. """
    n = peakdict["ex"].size
    sitenums = np.zeros(n)
    for i in range(n):
        p = ax.scatter(
            peakdict["ex"][i],
            peakdict["em"][i],
            s=20,
            facecolors="none",
            edgecolors="r",
        )
        plt.draw()
        answer = input("Which site?")
        try:
            site = int(answer)
            sitenums[i] = site
            print("Peak added to site %d.\n" % site)
            p.remove()
        except ValueError:
            print("Peak removed.")
            p.remove()

    indexlist = [np.equal(sitenums, i) for i in range(max(sitenums) + 1)]
    sitewisepeakdicts = [
        {"ex": peakdict["ex"][ind], "em": peakdict["em"][ind]} for ind in indexlist[1:]
    ]

    return sitewisepeakdicts


def plot2d(datadict, clevels=None, figsize=(6, 6), inputfigure=None):
    """Plots a 2d spectrum without sideplots."""
    excitation = datadict["ex"]
    emission = datadict["em"]
    spectrum = datadict["spec"]

    while len(excitation) > 2000:
        excitation = excitation[0::2]
        spectrum = spectrum[0::2, :]

    spectrum = np.log(spectrum)

    if inputfigure is None:
        fig = plt.figure(figsize=figsize)
        ax1 = plt.gca()
    else:
        fig = inputfigure
        ax1 = fig.axes[0]

    if clevels is None:
        clevels = np.linspace(np.amin(spectrum), np.amax(spectrum), 100)

    ax1.contourf(excitation, emission, spectrum, levels=clevels, cmap="jet")
    ax1.tick_params(
        axis="both",
        which="both",
        bottom=True,
        top=False,
        left=False,
        right=True,
        labelbottom=True,
        labelleft=True,
        labelright=False,
    )
    ax1.set_xlabel("Excitation Wavelength (nm)")
    ax1.set_ylabel("Emission Wavelength (nm)")

    ax1.margins(x=0, y=0)

    ax1.set_facecolor(cm.jet(0))

    return fig, ax1
