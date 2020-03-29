import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
import os
import peakfinding


def hist_top(fig, excitation, expeaks, linewidth):

    """ Plots a histogram of peaks on the top axis of a figure. """

    axes = fig.axes
    topax = axes[1]
    topax.cla()
    exbins = get_bins(excitation, linewidth)
    topax.hist(expeaks, bins=exbins)
    topax.margins(0)

    return exbins


def hist_side(fig, emission, empeaks, linewidth):

    """ Plots a histogram of peaks on the left axis of a figure. """

    axes = fig.axes
    sideax = axes[2]
    sideax.cla()
    embins = get_bins(emission, linewidth)
    sideax.hist(empeaks, bins=embins, orientation='horizontal')
    a, b = sideax.get_xlim()
    sideax.set_xlim(b, a)
    sideax.margins(0)

    return embins


def get_bins(wavelengths, linewidth):

    """ Returns a list of bin edges covering the given wavelength
    ranges, with bins of the given linewidth. """

    return np.linspace(min(wavelengths),
                       max(wavelengths),
                       num=int((max(wavelengths)
                                - min(wavelengths))/linewidth))


def occupancy_matrix(exbins, embins, expeaks, empeaks):

    """ Returns a matrix with binary occupancy values to say whether
    there is a peak in each bin. """

    omatrix = np.zeros((exbins.size-1, embins.size-1), dtype=bool)
    for i in range(expeaks.size):
        x = int(0)
        y = int(0)
        while (expeaks[i] > exbins[x]):
            x += 1
        while (empeaks[i] > embins[y]):
            y += 1
        omatrix[x-1, y-1] = 1

    return omatrix


def show_omatrix(ax, omatrix, exbins, embins):

    """ Plot the omatrix on top of a preexisting axis. """

    ax.imshow(omatrix.transpose(),
              cmap='winter',
              alpha=0.3,
              zorder=1,
              origin='lower',
              aspect='auto',
              extent=(min(exbins), max(exbins),
                      min(embins), max(embins)))


def show_groups(ax, col_groups, exbins, embins):

    """ Plot the groups on top of everything else. """

    ax.imshow([col_groups, col_groups],
              cmap="tab10",
              alpha=0.3,
              zorder=2,
              origin='lower',
              aspect='auto',
              extent=(min(exbins), max(exbins),
                      min(embins), max(embins)))


def trim_omatrix(omatrix, groups, nmissing):

    """ Reduces the True entries of the omatrix. All omatrix squares
    in Group 0 columns are set to False. Then, for each of the other
    Groups, the peaks that don't appear in enough of the columns are
    removed as well. "Enough" is defined by the number of columns minus
    the number you can miss (nmissing). """

    h = omatrix.shape[1]  # get the height of the omatrix
    omatrix[groups[0]] = 0  # set all the group 0 columns to 0

    for g in groups[1:]:
        for i in range(h):
            # omatrix[g, i] = 

    return omatrix


def dist(a, b):

    """ Returns the distance between two boolean vectors a and b. """

    return sum(np.bitwise_xor(a, b))


def group_by_columns(o, d):

    """ Groups the columns of the omatrix into clusters, where col1
    and col2 belong to the same cluster if they are separated by less
    than dist (or connected via a path like this.) """

    N = o.shape[0]
    col_groups = -np.ones(N, dtype=int)
    groups = []
    for i in range(N):
        if col_groups[i] == -1:
            # the column doesn't have a group yet
            g = []  # start a new group
            ind = len(groups)  # the index for the new group
            for j in range(N):
                if dist(o[i], o[j]) < d:
                    g = g + [j]
                    col_groups[j] = ind
            groups = groups + [g]
        else:
            # the column does have a group
            ind = col_groups[i]
            g = groups[ind]
            for j in range(N):
                if col_groups[j] < 0:
                    if dist(o[i], o[j]) < d:
                        g = g + [j]
                        col_groups[j] = ind

    return groups, col_groups


#  __main__ #

# hard coded variables
# folder = r'Jin'
# fname = r'Jin_MgO,EDFA=867mA,gain=100dB2Dsweep_5.mat'
folder = r'Anthemius'
fname = r'Anthemius_ZnS,EDFA=867mA,gain=100dB2Dsweep_1.mat'
lwidth = 0.2
nfraction = 0.05
iterations = 2
window = 2

# read data in
rawdata = sio.loadmat(os.path.join(folder, fname))

# trim and simplify data
spec_mask = (350, 100)  # (0, 1) for no trimming
calibratedemissionaxis = (np.linspace(1, 1024, 1024) -
                          -1.46581267e+04)/9.91909400e+00 + 20.
emission = calibratedemissionaxis[spec_mask[0]:-spec_mask[1]]
spectrum = rawdata['Spectrum2D'][spec_mask[0]:-spec_mask[1], :]
excitation = rawdata['ActualWavelength'][0]

# plot
expeaks, empeaks, spectrum = peakfinding.find_peaks(
    excitation,
    emission,
    spectrum,
    linewidth=lwidth,
    iterations=iterations,
    window=window,
    noisefraction=nfraction)

fig, ax = peakfinding.plot2d_sideplots(excitation, emission, spectrum)
peakfinding.plot_peaks(ax, expeaks, empeaks)

exbins = hist_top(fig, excitation, expeaks, lwidth)
embins = hist_side(fig, emission, empeaks, 1)
omatrix = occupancy_matrix(exbins, embins, expeaks, empeaks)

show_omatrix(ax, omatrix, exbins, embins)

plt.close()
# plt.show()
plt.close()

groups, col_groups = group_by_columns(omatrix, 5)

fig, ax = peakfinding.plot2d(excitation, emission, spectrum)
peakfinding.plot_peaks(ax, expeaks, empeaks)
show_omatrix(ax, omatrix, exbins, embins)
show_groups(ax, col_groups, exbins, embins)

# plt.show()
plt.close()

omatrix = trim_omatrix(omatrix, groups, 0)
