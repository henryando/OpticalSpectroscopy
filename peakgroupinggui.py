import numpy as np
import matplotlib.pyplot as plt
import peakfindinggui as pfg
import spectrumtools as st


def _get_bins(wavelengths, linewidth):
    """ Returns a list of bin edges covering the given wavelength
    ranges, with bins of the given linewidth. """
    return np.linspace(
        min(wavelengths),
        max(wavelengths),
        num=int((max(wavelengths) - min(wavelengths)) / linewidth),
    )


def _occupancy_matrix(exbins, embins, expeaks, empeaks):
    """ Returns a matrix with binary occupancy values to say whether
    there is a peak in each bin.
    """
    omatrix = np.zeros((exbins.size - 1, embins.size - 1), dtype=bool)
    for i in range(expeaks.size):
        x = int(0)
        y = int(0)
        while expeaks[i] > exbins[x]:
            x += 1
        while empeaks[i] > embins[y]:
            y += 1
        omatrix[x - 1, y - 1] = 1

    return omatrix


def _get_omatrix(data, peaks, Wx, Wy):
    """Gets the occupancy matrix, given either a single spectrum or a list
    of many spectra.
    """
    if type(data) is list:
        exmins = np.zeros(len(data))
        exmaxs = np.zeros(len(data))
        emmins = np.zeros(len(data))
        emmaxs = np.zeros(len(data))
        for i in range(len(data)):
            exmins[i] = min(data[i].ex)
            exmaxs[i] = max(data[i].ex)
            emmins[i] = min(data[i].em)
            emmaxs[i] = max(data[i].em)
        exmin = min(exmins)
        exmax = max(exmaxs)
        emmin = min(emmins)
        emmax = max(emmaxs)
    else:
        exmin = min(data.ex)
        exmax = max(data.ex)
        emmin = min(data.em)
        emmax = max(data.em)

    exbins = _get_bins((exmin, exmax), Wx)
    embins = _get_bins((emmin, emmax), Wy)
    omatrix = _occupancy_matrix(exbins, embins, peaks.ex, peaks.em)
    return omatrix, exbins, embins


def _get_occupied_cols(omatrix):
    """Returns a list of the occupied columns in the omatrix."""
    return (np.sum(omatrix, axis=1) > 0).astype(int)


def _get_occupied_rows(omatrix):
    """Returns a list of the occupied rows in the omatrix."""
    return (np.sum(omatrix, axis=0) > 0).astype(int)


class PeakGrouper:
    def __init__(self, data, peaks, Wx, Wy, ax=None):
        self.data = data
        self.peaks = peaks
        self.omatrix, exbins, embins = _get_omatrix(data, peaks, Wx, Wy)
        self.exbins = exbins
        self.embins = embins
        self.cols = _get_occupied_cols(self.omatrix)
        self.rows = _get_occupied_rows(self.omatrix)
        if ax is None:
            self.fig = plt.gcf()
            self.ax = plt.gca()
        else:
            self.ax = ax
            self.fig = ax.figure
        self.state = None  # True if grouping excitation, False for emission
        self.imhandle = None  # Handle to the image currently displayed.

        self.cid = self.fig.canvas.mpl_connect("button_press_event", self)

    def __call__(self, event):
        if self.state is None:
            return
        elif self.state is True:
            col = self._mousex_to_col(event.xdata)
            self.cols[col] = 2
            self.show_cols()
        elif self.state is False:
            row = self._mousey_to_row(event.ydata)
            self.rows[row] = 2
            self.show_rows()
        self.fig.canvas.draw()

    def pick_cols(self):
        """Lets you pick the columns you want."""
        self.state = True
        self.show_cols()
        plt.show()
        self.state = None
        lines = self.remove_extra_peaks_col()
        self.rows = _get_occupied_rows(self.omatrix)
        self.cols = _get_occupied_cols(self.omatrix)
        return lines

    def pick_rows(self):
        """Lets you pick the rows you want."""
        self.state = False
        self.show_rows()
        plt.show()
        self.state = None
        lines = self.remove_extra_peaks_row()
        self.rows = _get_occupied_rows(self.omatrix)
        self.cols = _get_occupied_cols(self.omatrix)
        return lines

    def remove_extra_peaks_col(self):
        """Remove the peaks except those in the columns you selected.
        Also returns the central wavelengths of the peaks in each bin.
        """
        goodcols = np.nonzero(self.cols * (self.cols - 1))[0]
        N = len(goodcols)
        wavelengths = np.zeros(N)
        goodpeaks = st.Peaks()
        for i in range(N):
            peakinds = np.bitwise_and(
                self.peaks.ex > self.exbins[goodcols[i]],
                self.peaks.ex < self.exbins[goodcols[i] + 1],
            )
            if sum(peakinds > 0):
                wavelengths[i] = np.mean(self.peaks.ex[peakinds])
            else:
                wavelengths[i] = np.nan
            goodpeaks.ex = np.concatenate(
                (goodpeaks.ex, self.peaks.ex[peakinds]), axis=0
            )
            goodpeaks.em = np.concatenate(
                (goodpeaks.em, self.peaks.em[peakinds]), axis=0
            )
        self.peaks = goodpeaks
        self.omatrix[self.cols < 2] = 0
        return wavelengths

    def remove_extra_peaks_row(self):
        """Remove the peaks except those in the rows you selected.
        Also returns the central wavelengths of the peaks in each bin.
        """
        goodrows = np.nonzero(self.rows * (self.rows - 1))[0]
        N = len(goodrows)
        wavelengths = np.zeros(N)
        goodpeaks = st.Peaks()
        for i in range(N):
            peakinds = np.bitwise_and(
                self.peaks.em > self.embins[goodrows[i]],
                self.peaks.em < self.embins[goodrows[i] + 1],
            )
            if sum(peakinds > 0):
                wavelengths[i] = np.mean(self.peaks.em[peakinds])
            else:
                wavelengths[i] = np.nan
            goodpeaks.ex = np.concatenate(
                (goodpeaks.ex, self.peaks.ex[peakinds]), axis=0
            )
            goodpeaks.em = np.concatenate(
                (goodpeaks.em, self.peaks.em[peakinds]), axis=0
            )
        self.peaks = goodpeaks
        self.omatrix[:, self.rows < 2] = 0
        return wavelengths

    def show_omatrix(self, ax=None):
        """Plot the omatrix on top of a preexisting axis."""
        if ax is None:
            ax = plt.gca()
        ax.imshow(
            self.omatrix.transpose(),
            cmap="winter",
            alpha=0.3,
            zorder=1,
            origin="lower",
            aspect="auto",
            extent=(
                min(self.exbins),
                max(self.exbins),
                min(self.embins),
                max(self.embins),
            ),
        )

    def show_cols(self):
        """Plot the unclicked cols in red, and the clicked ones in green."""
        if self.imhandle is not None:
            self.imhandle.remove()
        self.imhandle = self.ax.imshow(
            [self.cols, self.cols],
            cmap="brg",
            vmin=0,
            vmax=2,
            alpha=0.2,
            zorder=2,
            origin="lower",
            aspect="auto",
            extent=(
                min(self.exbins),
                max(self.exbins),
                min(self.embins),
                max(self.embins),
            ),
        )

    def show_rows(self):
        """Plot the unclicked rows in red, and the clicked ones in green."""
        if self.imhandle is not None:
            self.imhandle.remove()
        self.imhandle = self.ax.imshow(
            np.transpose([self.rows, self.rows]),
            cmap="brg",
            vmin=0,
            vmax=2,
            alpha=0.2,
            zorder=2,
            origin="lower",
            aspect="auto",
            extent=(
                min(self.exbins),
                max(self.exbins),
                min(self.embins),
                max(self.embins),
            ),
        )

    def _mousex_to_col(self, x):
        """Convert mouse coordinate to column number."""
        col = 0
        while x > self.exbins[col]:
            col += 1
        return col - 1

    def _mousey_to_row(self, y):
        """Convert mouse coordinate to row number."""
        row = 0
        while y > self.embins[row]:
            row += 1
        return row - 1


def snap_exline(data, ex):
    """Snaps an ex line to the exact value of the local excitation max."""
    if type(data) is list:
        for d in data:
            snapped = snap_exline(d, ex)
            if snapped is not None:
                return snapped
    else:
        if (min(data.ex) <= ex) and (ex <= max(data.ex)):
            return data.ex[(sum(data.ex > ex) + 1)]


def snap_emline(data, em):
    """Snaps an em line to the exact value of the nearest emission point."""
    if type(data) is list:
        for d in data:
            snapped = snap_emline(d, em)
            if snapped is not None:
                return snapped
    else:
        return data.em[(sum(data.em > em) + 1)]


####################################################################################
# The main function for this module.
####################################################################################


def find_lines(spectra, peaks, grouping_Wx=0.5, grouping_Wy=1):
    if type(peaks) is not st.Peaks:
        peaks = st.Peaks(peaks[0], peaks[1])

    fig, ax = pfg.plot2d(spectra)
    pfg.plot_peaks(peaks.ex, peaks.em)
    pg = PeakGrouper(spectra, peaks, grouping_Wx, grouping_Wy)
    pg.show_omatrix()
    try:
        exlines = pg.pick_cols()
    except RuntimeWarning:
        pass

    peaks = pg.peaks
    fig, ax = pfg.plot2d(spectra)
    pfg.plot_peaks(peaks.ex, peaks.em)
    pg = PeakGrouper(spectra, peaks, grouping_Wx, grouping_Wy)
    pg.show_omatrix()
    try:
        emlines = pg.pick_rows()
    except RuntimeWarning:
        pass

    exlines = exlines[np.bitwise_not(np.isnan(exlines))]
    emlines = emlines[np.bitwise_not(np.isnan(emlines))]

    # exlines = [snap_exline(spectra, l) for l in exlines]
    # emlines = [snap_emline(spectra, l) for l in emlines]

    return exlines, emlines
