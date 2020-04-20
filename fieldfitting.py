import numpy as np
import matplotlib.pyplot as plt
import spectroscopymain as sm
import spectrumtools as st
import lineanalysis as la
import cubicfields as cf


def _processed_fieldmodel_levels(x, W, J, spacing, minimum, ignorelevels=0):
    eigvals = np.unique(np.round(cf.energies(J, W, x), decimals=8))
    eigvals = eigvals * spacing / (eigvals[1 + ignorelevels] - eigvals[0])
    eigvals = eigvals - eigvals[0] + minimum
    return eigvals


def energy_dist(obs, calc):
    d = 0
    for o in obs:
        d += min((o - calc) ** 2)
    return d


def field_fit(levels, J, wsign=1, ignorelevels=0):
    """Plots field fit, and returns best fit parameters W and x. Takes as inputs
    the energy levels as a list and the total angular momentum J.
    """
    xmin = -1
    xmax = 1
    xvals = np.linspace(xmin, xmax, num=2000)
    levels = np.sort(levels)
    spacing = levels[1] - levels[0]
    minimum = levels[0]
    energies = [
        _processed_fieldmodel_levels(x, wsign, J, spacing, minimum, ignorelevels)
        for x in xvals
    ]
    length = round(np.mean([len(e) for e in energies]))
    good = np.zeros_like(xvals, dtype=bool)
    for i in range(len(good)):
        good[i] = len(energies[i]) == length
    xvals = xvals[good]
    energies = [e for e in energies if len(e) == length]
    energies = np.asarray(energies)
    dists = [energy_dist(levels, e) for e in energies]

    g = np.argmin(np.asarray(dists))
    plt.figure()
    plt.plot(xvals, energies, "k", linewidth=1)
    for l in levels:
        plt.plot((xmin, xmax), (l, l), "b--", linewidth=1)
    plt.ylim(
        (
            min(levels) - (levels[1] - levels[0]) / 5,
            max(levels + 50) + (levels[1] - levels[0]) / 4,
        )
    )
    plt.plot((xvals[g], xvals[g]), plt.ylim(), "r:")
    plt.show()

    x = xvals[g]
    dmin = dists[g]
    eigvals = np.unique(np.round(cf.energies(J, wsign, x), decimals=8))
    w = wsign * spacing / (eigvals[1 + ignorelevels] - eigvals[0])
    return w, x, dmin


def get_levels(W, x, J):
    eigvals = np.unique(np.round(cf.energies(J, W, x), decimals=8))
    return eigvals - min(eigvals)
