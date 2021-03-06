import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
import spectrumtools as st


class PeakClicker:
    def __init__(self, data, peaks):
        self.datadict = data
        self.xpeaks = peaks.ex
        self.ypeaks = peaks.em
        self.N = self.xpeaks.size
        self.goodpeaks = np.zeros(self.N, dtype=bool)
        self.fig, self.ax = plot2d(data)
        (xmin, xmax) = plt.xlim()
        (ymin, ymax) = plt.ylim()
        self.ar = (xmax - xmin) ** 2 / (ymax - ymin) ** 2
        plot_peaks(self.xpeaks, self.ypeaks, ax=self.ax)
        self.cid = self.fig.canvas.mpl_connect("button_press_event", self)

    def __call__(self, event):
        x = event.xdata
        y = event.ydata
        print("x = %.2f, y = %.2f" % (x, y))
        dist = np.zeros(self.N)
        for i in range(self.N):
            dist[i] = (x - self.xpeaks[i]) ** 2 + self.ar * (y - self.ypeaks[i]) ** 2
        dist = dist - np.amin(dist)
        self.goodpeaks[dist == 0] = 1
        _plot_good_peaks(self.xpeaks[self.goodpeaks], self.ypeaks[self.goodpeaks])
        self.fig.canvas.draw()

    def good_peaks(self):
        return st.Peaks(self.xpeaks[self.goodpeaks], self.ypeaks[self.goodpeaks])


def plot2d(data, clevels=None, figsize=(7, 7), inputfigure=None, xpixels=2000):
    """Plots a 2d spectrum without sideplots."""
    if type(data) is list:
        fig = plt.figure(figsize=figsize)
        ax = plt.subplot(111)
        if clevels is None:
            a = min([np.min(d.spec) for d in data])
            b = max([np.max(d.spec) for d in data])
            a = max((a, 100))
            b = 0.8 * (b - a) + a
            clevels = np.linspace(np.log(a), np.log(b), 100)
        for d in data:
            plot2d(
                d, clevels=clevels, inputfigure=fig, xpixels=int(xpixels / len(data)),
            )
        return fig, ax

    excitation = data.ex
    emission = data.em
    spectrum = data.spec

    while len(excitation) > xpixels:
        excitation = excitation[0::2]
        spectrum = spectrum[:, 0::2]

    spectrum = np.log(spectrum)

    if inputfigure is None:
        fig = plt.figure(figsize=figsize)
        ax = plt.gca()
    else:
        fig = inputfigure
        ax = fig.axes[0]

    if clevels is None:
        clevels = np.linspace(np.amin(spectrum), np.amax(spectrum), 100)

    ax.contourf(excitation, emission, spectrum, levels=clevels, cmap="jet")
    ax.tick_params(
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
    ax.set_xlabel(r"Excitation Wavelength / cm$^{-1}$")
    ax.set_ylabel(r"Emission Wavelength / cm$^{-1}$")
    ax.set_facecolor(cm.jet(0))
    return fig, ax


def plot_peaks(xdata, ydata, ax=None):
    """Formatting in only one spot."""
    if ax is None:
        ax = plt.gca()
    ax.scatter(xdata, ydata, s=1, facecolors="none", edgecolors="r")


def _plot_good_peaks(xdata, ydata, ax=None):
    """Formatting in only one spot."""
    if ax is None:
        ax = plt.gca()
    ax.scatter(xdata, ydata, s=20, facecolors="g", edgecolors="g")


####################################################################################
# The main function for this module.
####################################################################################
def find_peaks(spectra, linewidth=0.07, nf=1 / 12):
    goodpeaks = st.Peaks()
    for s in spectra:
        peaks = st.find_potential_peaks(s, linewidth=linewidth, noisefraction=nf)
        print("There are %d peaks." % peaks.em.size)
        pc = PeakClicker(s, peaks)
        print("Double click the good peaks, then quit the plot.")
        plt.show()
        goodpeaks.ex = np.concatenate((goodpeaks.ex, pc.good_peaks().ex), axis=0)
        goodpeaks.em = np.concatenate((goodpeaks.em, pc.good_peaks().em), axis=0)
    plot2d(spectra)
    _plot_good_peaks(goodpeaks.ex, goodpeaks.em)
    plt.show()
    return goodpeaks
