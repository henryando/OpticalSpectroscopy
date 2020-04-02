import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
import spectrumtools as st


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


def plot2d(datadict, clevels=None, figsize=(7, 7), inputfigure=None, xpixels=2000):
    """Plots a 2d spectrum without sideplots."""
    if type(datadict) is list:
        fig = plt.figure(figsize=figsize)
        ax = plt.subplot(111)
        for d in datadict:
            plot2d(
                d,
                clevels=clevels,
                inputfigure=fig,
                xpixels=int(xpixels / len(datadict)),
            )
        return fig, ax

    excitation = datadict["ex"]
    emission = datadict["em"]
    spectrum = datadict["spec"]

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
    ax.set_xlabel("Excitation Wavelength / nm")
    ax.set_ylabel("Emission Wavelength / nm")
    ax.set_facecolor(cm.jet(0))
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


def plot_lines(lines, ax=None):
    """Plots the emission and excitation lines."""
    if ax is None:
        ax = plt.gca()
    exlines = lines[0]
    emlines = lines[1]
    xlim = plt.xlim()
    ylim = plt.ylim()
    ax.plot((exlines, exlines), ylim, "r", linewidth=0.5)
    ax.plot(xlim, (emlines, emlines), "r", linewidth=0.5)


####################################################################################
# The main function for this module.
####################################################################################
def find_peaks(spectra, linewidth=0.07, nf=1 / 12):
    goodpeaks = {"ex": [], "em": []}
    for s in spectra:
        peakdict = st.find_peaks(s, linewidth=linewidth, noisefraction=nf)
        print("There are %d peaks." % peakdict["em"].size)
        pc = PeakClicker(s, peakdict)
        print("Double click the good peaks, then quit the plot.")
        plt.show()
        goodpeaks["ex"] = np.concatenate(
            (goodpeaks["ex"], pc.good_peaks()["ex"]), axis=0
        )
        goodpeaks["em"] = np.concatenate(
            (goodpeaks["em"], pc.good_peaks()["em"]), axis=0
        )
    plot2d(spectra)
    plot_good_peaks(goodpeaks["ex"], goodpeaks["em"])
    plt.show()
    goodpeaks["ex"] = goodpeaks["ex"][0]
    goodpeaks["em"] = goodpeaks["em"][0]
    return goodpeaks
