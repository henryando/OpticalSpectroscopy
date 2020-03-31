import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm


class PeakClicker:
    def __init__(self, datadict, peakdict):
        self.datadict = datadict
        self.xpeaks = peakdict["ex"]
        self.ypeaks = peakdict["em"]
        self.N = self.xpeaks.size
        self.goodpeaks = np.zeros(self.N, dtype=bool)
        self.fig, self.ax = plot2d(datadict)
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
        plot_good_peaks(self.xpeaks[self.goodpeaks], self.ypeaks[self.goodpeaks])
        self.fig.canvas.draw()

    def good_peaks(self):
        return {"ex": self.xpeaks[self.goodpeaks], "em": self.ypeaks[self.goodpeaks]}


def plot2d(datadict, clevels=None, figsize=(7, 7), inputfigure=None):
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
    ax.set_xlabel("Excitation Wavelength / nm")
    ax.set_ylabel("Emission Wavelength / nm")
    ax.set_facecolor(cm.jet(0))
    return fig, ax


def plot2d_list(spectra, clevels=None, figsize=(7, 7), inputfigure=None):
    """Does plot2d for many spectra."""
    fig = plt.figure(figsize=figsize)
    ax = plt.subplot(111)
    for s in spectra:
        plot2d(s, clevels=clevels, inputfigure=fig)
    return fig, ax


def plot_peaks(xdata, ydata, ax=None):
    """Formatting in only one spot."""
    if ax is None:
        ax = plt.gca()
    ax.scatter(xdata, ydata, s=1, facecolors="none", edgecolors="r")


def plot_good_peaks(xdata, ydata, ax=None):
    """Formatting in only one spot."""
    if ax is None:
        ax = plt.gca()
    ax.scatter(xdata, ydata, s=20, facecolors="g", edgecolors="g")
