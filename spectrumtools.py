import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


def iterative_smooth(spectrum, N, W):
    """ Smooths a spectrum in both directions N times with a moving
    window of 2W+1 points.
    """
    avg = np.mean(np.mean(spectrum))
    sz = spectrum.shape
    target = np.ones((sz[0] + 2 * W, sz[1] + 2 * W)) * avg
    target[W:-W, W:-W] = spectrum

    for i in range(N):
        for j in range(sz[0]):
            target[j + W, :] = np.mean(target[j : (j + 1 + 2 * W), :], axis=0)
        for j in range(sz[1]):
            target[:, j + W] = np.mean(target[:, j : (j + 1 + 2 * W)], axis=1)

    return target[W:-W, W:-W]


def find_local_maxes_2d(data):
    """ Finds the peaks (local max in both directions) of a 2 dimensional data set.
    Use np.nonzero to get indices.
    """
    xpeaks = np.bitwise_and(data[:-2, :] < data[1:-1, :], data[1:-1, :] > data[2:, :])
    ypeaks = np.bitwise_and(data[:, :-2] < data[:, 1:-1], data[:, 1:-1] > data[:, 2:])
    peaks = np.zeros_like(data, dtype=bool)
    peaks[1:-1, 1:-1] = np.bitwise_and(xpeaks[:, 1:-1], ypeaks[1:-1, :])
    return peaks


def linewidth_to_nsamples(wavelengths, linewidth):
    """ Given an array of wavelengths like emission or excitation, return the
    number of samples that corresponds to the appropriate linewidth.
    """
    return int(
        np.ceil(linewidth * wavelengths.size / (wavelengths[-1] - wavelengths[0]))
    )


def filter_peaks_square_mean(data, peaks, xrad, yrad, height=0):
    """ Given a data array and a set of peaks in that array, filter the peaks
    essentially by whether they are large enough to be plausible (are they greater
    than the mean of a certain square surrounding them?)
    """
    (xpeaks, ypeaks) = np.nonzero(peaks)
    (xmax, ymax) = peaks.shape

    for i in range(xpeaks.size):
        if ((xpeaks[i] > xrad) and (xpeaks[i] < xmax - xrad)) and (
            (ypeaks[i] > yrad) and (ypeaks[i] < ymax - yrad)
        ):
            # is it larger than the mean of the square perimeter?
            perimeter = (data[xpeaks[i], ypeaks[i]] - height) > (
                np.mean(data[xpeaks[i] - xrad : xpeaks[i] + xrad, ypeaks[i] - yrad])
                + np.mean(data[xpeaks[i] - xrad : xpeaks[i] + xrad, ypeaks[i] + yrad])
                + np.mean(data[xpeaks[i] - xrad, ypeaks[i] - yrad : ypeaks[i] + yrad])
                + np.mean(data[xpeaks[i] + xrad, ypeaks[i] - yrad : ypeaks[i] + yrad])
            ) / 4
            # is it a max relative to the corners of the box?
            linelocalmax = data[xpeaks[i], ypeaks[i]] > max(
                (
                    data[xpeaks[i] - xrad, ypeaks[i] - yrad],
                    data[xpeaks[i] - xrad, ypeaks[i] + yrad],
                    data[xpeaks[i] + xrad, ypeaks[i] - yrad],
                    data[xpeaks[i] + xrad, ypeaks[i] + yrad],
                )
            )
            # filter the peak
            peaks[xpeaks[i], ypeaks[i]] = perimeter and linelocalmax

        else:
            peaks[xpeaks[i], ypeaks[i]] = 0

    return peaks


def integrate_emission(spectrum):
    """ Integrates along the emission axis for plot. Assumes spectrum is log
    already.
    """
    return np.sum(np.exp(spectrum), axis=1)


def integrate_excitation(spectrum):
    """ Integrates along the excitation axis for plot. Assumes spectrum is log
    already.
    """
    return np.sum(np.exp(spectrum), axis=0)


def plot2d_sideplots(
    excitation,
    emission,
    spectrum,
    clevels=None,
    figsize=(6, 6),
    nrows=5,
    ncols=5,
    inputfigure=None,
):
    """ Plots a 2d spectrum with integrated side plots showing the total excitation
    and emission spectra.
    """
    while len(excitation) > 2000:
        excitation = excitation[0::2]
        spectrum = spectrum[0::2, :]

    if clevels is None:
        clevels = np.linspace(np.amin(spectrum), np.amax(spectrum), 100)

    if inputfigure is None:
        fig = plt.figure(figsize=figsize)
        spec = gridspec.GridSpec(ncols=ncols, nrows=nrows, figure=fig)
        ax1 = fig.add_subplot(spec[1:, 1:])
        ax2 = fig.add_subplot(spec[0, 1:], sharex=ax1)
        ax3 = fig.add_subplot(spec[1:, 0], sharey=ax1)
        fig.subplots_adjust(hspace=0, wspace=0)
    else:
        fig = inputfigure
        axes = fig.axes
        ax1 = axes[0]
        ax2 = axes[1]
        ax3 = axes[2]

    ax1.contourf(excitation, emission, spectrum, levels=clevels, cmap="jet")
    ax1.tick_params(
        axis="both",
        which="both",
        bottom=True,
        top=False,
        left=False,
        right=True,
        labelbottom=True,
        labelleft=False,
        labelright=True,
    )
    ax1.set_xlabel("Excitation Wavelength (nm)")
    ax3.set_ylabel("Emission Wavelength (nm)")

    ax2.plot(excitation, integrate_excitation(spectrum))
    ax2.tick_params(
        axis="both",
        which="both",
        bottom=False,
        top=True,
        right=False,
        left=False,
        labelbottom=False,
        labeltop=True,
        labelleft=False,
        labelright=False,
    )

    ax3.plot(-integrate_emission(spectrum), emission)
    ax3.tick_params(
        axis="both",
        which="both",
        right=False,
        left=True,
        bottom=False,
        top=False,
        labelbottom=False,
        labeltop=False,
        labelleft=True,
        labelright=False,
    )

    ax1.margins(x=0, y=0)

    return fig, ax1


def plot2d(
    excitation, emission, spectrum, clevels=None, figsize=(6, 6), inputfigure=None
):
    """ Plots a 2d spectrum without sideplots. """
    if inputfigure is None:
        fig = plt.figure(figsize=figsize)
        ax1 = plt.gca()
    else:
        fig = inputfigure
        ax1 = fig.axes[0]

    while len(excitation) > 2000:
        excitation = excitation[0::2]
        spectrum = spectrum[0::2, :]

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

    return fig, ax1


def find_peaks(
    excitation,
    emission,
    spectrum,
    window=3,
    iterations=2,
    linewidth=0.2,
    noisefraction=(1 / 16),
):
    """ Takes the emission wavelength vector, excitation wavelength vector,
    spectrum (in the original scale, not log), and returns two vectors, expeaks
    and empeaks, which are the coordinates of the peaks that it found. Optional
    arguments control the moving window size (2*window+1), the number of times to
    smooth (iterations), the linewidth to force the points to be maximal over
    (2*linewidth), and the fraction of the noise level to force the peaks to be
    greater by (noisefraction).
    """
    spectrum = np.log(spectrum)
    spectrum = iterative_smooth(spectrum, iterations, window)
    peaks = find_local_maxes_2d(spectrum)
    emrad = linewidth_to_nsamples(emission, linewidth)
    exrad = linewidth_to_nsamples(excitation, linewidth)
    noise = np.std(spectrum)
    peaks = filter_peaks_square_mean(
        spectrum, peaks, exrad, emrad, height=(noise * noisefraction)
    )
    (xpeaks, ypeaks) = np.nonzero(peaks)
    peakdict = {"ex": excitation[ypeaks], "em": emission[xpeaks]}

    return peakdict, spectrum


def plot_peaks(ax, peakdict):
    """ Plots peaks on the given axis. """
    ax.plot(peakdict["ex"], peakdict["em"], "ro", markersize=5)


def manual_peak_sort(ax, peakdict):
    """ Lets the user filter the peaks manually. """
    n = peakdict["ex"].size
    sitenums = np.zeros(n)
    for i in range(n):
        ax.plot(peakdict["ex"][i], peakdict["em"][i], "ro", markersize=5)
        plt.draw()
        answer = input("Which site?\n")
        try:
            site = int(answer)
            sitenums[i] = site
            print("Peak added to site %d.\n" % site)
        except ValueError:
            print("Peak removed.\n")

    indexlist = [np.equal(sitenums, i) for i in range(max(sitenums) + 1)]
    sitewisepeakdicts = [
        {"ex": peakdict["ex"][ind], "em": peakdict["em"][ind]} for ind in indexlist[1:]
    ]

    return sitewisepeakdicts
