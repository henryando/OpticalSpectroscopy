import numpy as np
import spectrumtools as st
import conversions as conv
import matplotlib.pyplot as plt
import pickle
import spectroscopymain as sm
from scipy.optimize import curve_fit


def gauss(x, a, b, w, x0, s):
    return a * np.exp(-((x - x0) ** 2) / (2 * w ** 2)) + b + x * s


def lorentz(x, a, b, w, x0, s):
    return a / (1 + ((x - x0) ** 2 / w)) + b + x * s


def rsquared(f, popt, x, y):
    return 1 - np.sum((f(x, *popt) - y) ** 2) / np.sum((y - np.mean(y)) ** 2)


def fit_expeaks(excitation, spectrum, ela, width=1.5):
    i = 1
    for x in ela.ylevels + ela.z1y1:
        goodinds = np.bitwise_and(excitation > x - width, excitation < x + width)
        ex = excitation[goodinds] - x
        spec = spectrum[goodinds]
        p0 = (max(spec), min(spec), 0.05, 0, 0)
        gopt, _ = curve_fit(gauss, ex, spec, p0=p0)
        lopt, _ = curve_fit(lorentz, ex, spec, p0=p0)
        gfwhm = 2.355 * gopt[2]
        lfwhm = 2 * np.sqrt(lopt[2])
        gr2 = rsquared(gauss, gopt, ex, spec)
        lr2 = rsquared(lorentz, lopt, ex, spec)
        print("%d & %.2f & %.3f & %.2f & %.3f \\\\" % (i, lfwhm, lr2, gfwhm, gr2))
        if lr2 > gr2:
            print("$%.2f$ &" % lfwhm)
        else:
            print("$(%.2f)$ &" % gfwhm)
        i = i + 1
        plt.figure()
        plt.plot(ex, spec)
        plt.plot(ex, gauss(ex, *gopt))
        plt.plot(ex, lorentz(ex, *lopt))
        plt.legend(("Data", "Gaussian", "Lorentzian"))
        plt.show()

    return
