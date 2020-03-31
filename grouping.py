import numpy as np
import matplotlib.pyplot as plt


def _dist(a, b):
    """Returns the distance between two boolean vectors a and b."""
    return sum(np.bitwise_xor(a, b))


def get_bins(wavelengths, linewidth):
    """ Returns a list of bin edges covering the given wavelength
    ranges, with bins of the given linewidth. """
    return np.linspace(
        min(wavelengths),
        max(wavelengths),
        num=int((max(wavelengths) - min(wavelengths)) / linewidth),
    )


def occupancy_matrix(exbins, embins, expeaks, empeaks):
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


def group_by_columns(o, d):
    """ Groups the columns of the omatrix into clusters, where col1
    and col2 belong to the same cluster if they are separated by less
    than dist (or connected via a path like this.)
    """
    N = o.shape[0]
    col_groups = -np.ones(N, dtype=int)
    groups = []
    for i in range(N):
        if col_groups[i] == -1:
            # the column doesn't have a group yet
            g = []  # start a new group
            ind = len(groups)  # the index for the new group
            for j in range(N):
                if _dist(o[i], o[j]) < d:
                    g = g + [j]
                    col_groups[j] = ind
            groups = groups + [g]
        else:
            # the column does have a group
            ind = col_groups[i]
            g = groups[ind]
            for j in range(N):
                if col_groups[j] < 0:
                    if _dist(o[i], o[j]) < d:
                        g = g + [j]
                        col_groups[j] = ind

    return groups, col_groups


# def group_by_columns(o, d):
#     """ Groups the columns of the omatrix into clusters, where col1
#     and col2 belong to the same cluster if they are separated by less
#     than dist (or connected via a path like this.)
#     """
#     N = o.shape[0]
#     col_groups = -np.ones(N, dtype=int)
#     groups = []
#     for i in range(N):
#         if col_groups[i] == -1:
#             # the column doesn't have a group yet
#             g = [] # start a new group
#             ind = len(groups) # the length


def get_omatrix(data, peakdict, Wx, Wy):
    """Gets the occupancy matrix, given either a single spectrum or a list
    of many spectra.
    """
    if type(data) is list:
        exmins = np.zeros(len(data))
        exmaxs = np.zeros(len(data))
        emmins = np.zeros(len(data))
        emmaxs = np.zeros(len(data))
        for i in range(len(data)):
            exmins[i] = min(data[i]["ex"])
            exmaxs[i] = max(data[i]["ex"])
            emmins[i] = min(data[i]["em"])
            emmaxs[i] = max(data[i]["em"])
        exmin = min(exmins)
        exmax = max(exmaxs)
        emmin = min(emmins)
        emmax = max(emmaxs)
    else:
        exmin = min(data["ex"])
        exmax = max(data["ex"])
        emmin = min(data["em"])
        emmax = max(data["em"])

    exbins = get_bins((exmin, exmax), Wx)
    embins = get_bins((emmin, emmax), Wy)
    omatrix = occupancy_matrix(exbins, embins, peakdict["ex"], peakdict["em"])
    return omatrix, {"exbins": exbins, "embins": embins}


class PeakGrouper:
    """Just a convenient way to store all the relevant data associated with a
    grouping procedure.
    """

    def __init__(self, data, peakdict, Wx, Wy, bool_dist):
        self.omatrix, self.bins = get_omatrix(data, peakdict, Wx, Wy)
        self.groups, self.col_groups = group_by_columns(self.omatrix, bool_dist)

    def trim_omatrix(self, nmissing):
        """ Reduces the True entries of the omatrix. All omatrix squares
        in Group 0 columns are set to False. Then, for each of the other
        Groups, the peaks that don't appear in enough of the columns are
        removed as well. "Enough" is defined by the number of columns minus
        the number you can miss (nmissing). """
        omatrix = self.omatrix
        groups = self.groups
        h = omatrix.shape[1]  # get the height of the omatrix
        # omatrix[groups[0]] = 0  # set all the group 0 columns to 0
        for g in groups[1:]:
            for i in range(h):
                omatrix[g, i] = sum(omatrix[g, i]) > (len(g) - nmissing)
        self.omatrix = omatrix
        return

    def show_omatrix(self, ax=None):
        """ Plot the omatrix on top of a preexisting axis. """
        if ax is None:
            ax = plt.gca()
        exbins = self.bins["exbins"]
        embins = self.bins["embins"]
        ax.imshow(
            self.omatrix.transpose(),
            cmap="winter",
            alpha=0.3,
            zorder=1,
            origin="lower",
            aspect="auto",
            extent=(min(exbins), max(exbins), min(embins), max(embins)),
        )

    def show_groups(self, ax=None):
        """ Plot the groups on top of everything else. """
        if ax is None:
            ax = plt.gca()
        exbins = self.bins["exbins"]
        embins = self.bins["embins"]
        ax.imshow(
            [self.col_groups, self.col_groups],
            cmap="tab10",
            alpha=0.3,
            zorder=2,
            origin="lower",
            aspect="auto",
            extent=(min(exbins), max(exbins), min(embins), max(embins)),
        )
