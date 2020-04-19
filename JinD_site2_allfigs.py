import numpy as np
import matplotlib.pyplot as plt
import spectroscopymain as sm
import lineanalysis as la
import fieldfitting as ff
from scipy.optimize import curve_fit


ela = la.load_object("Data/JinD/site2.pkl")
fileformat = "Figures/JinD_site2_%dK.png"

folders = ("Data/JinD 0", "Data/JinD 1", "Data/JinD 2", "Data/JinD 3")
temps = (13, 32, 60, 106)
N = 2
Wx = 1
Wy = 2


def prep(i):
    folder = folders[i]
    fname = fileformat % temps[i]
    peaks = np.load(folder + "/peaks.npy")
    spectra = sm.read_all_2dspectra(folder)
    spectra = sm.iterative_smooth(spectra, N=N, Wx=Wx, Wy=Wy)
    return fname, peaks, spectra


# print level structure
if True:
    print("Z1Y1:")
    print(ela.z1y1)
    print("Z levels:")
    [print("%.1f" % l) for l in ela.zlevels]
    print("Y levels:")
    [print("%.1f" % l) for l in ela.ylevels]

if False:
    ff.field_fit(ela.zlevels, 15 / 2)
    ff.field_fit(ela.ylevels, 13 / 2)

# 13K
if False:
    fname, peaks, spectra = prep(0)
    sm.plot_spectra(spectra, figsize=(5, 5))

    ela.plot_exline(1, 1)
    ela.plot_exline(1, 2)
    ela.plot_exline(1, 3)
    ela.plot_exline(1, 4)
    ela.plot_emline(1, 1)
    ela.plot_emline(2, 1)
    ela.plot_emline(3, 1)
    ela.plot_emline(4, 1)
    ela.plot_emline(5, 1)

    sm.plot_peaks(peaks)
    plt.savefig(fname, dpi=800, bbox_inches="tight")
    plt.show()


# 30K
if False:
    fname, peaks, spectra = prep(1)
    sm.plot_spectra(spectra, figsize=(5, 5))

    ela.plot_exline(1, 1)
    ela.plot_exline(1, 2)
    ela.plot_exline(1, 3)
    ela.plot_emline(1, 1)
    ela.plot_emline(2, 1)
    ela.plot_emline(3, 1)
    ela.plot_emline(4, 1)
    ela.plot_emline(5, 1)

    ela.plot_emline(4, 2, color=1)
    ela.plot_emline(4, 3, color=1)
    ela.plot_emline(2, 2, color=1)
    ela.plot_emline(2, 3, color=1)
    ela.plot_emline(3, 2, color=1)
    ela.plot_emline(1, 2, color=1)

    sm.plot_peaks(peaks)
    plt.savefig(fname, dpi=800, bbox_inches="tight")
    plt.show()


# 60K
if False:
    fname, peaks, spectra = prep(2)
    sm.plot_spectra(spectra, figsize=(5, 5))

    ela.plot_exline(1, 1)
    ela.plot_exline(1, 2)
    ela.plot_exline(1, 3)
    ela.plot_emline(1, 1)
    ela.plot_emline(2, 1)
    ela.plot_emline(3, 1)
    ela.plot_emline(4, 1)
    ela.plot_emline(5, 1)

    ela.plot_emline(4, 2, color=1)
    ela.plot_emline(4, 3, color=1)
    ela.plot_emline(2, 2, color=1)
    ela.plot_emline(2, 3, color=1)
    ela.plot_emline(3, 2, color=1)
    ela.plot_emline(1, 2, color=1)
    ela.plot_emline(1, 3, color=1)
    ela.plot_emline(2, 4, color=1)

    sm.plot_peaks(peaks)
    plt.savefig(fname, dpi=800, bbox_inches="tight")
    plt.show()


# 100K
if False:
    fname, peaks, spectra = prep(3)
    sm.plot_spectra(spectra, figsize=(5, 5))

    ela.plot_exline(1, 1)
    ela.plot_exline(1, 2)
    ela.plot_exline(1, 3)
    ela.plot_emline(1, 1)
    ela.plot_emline(2, 1)
    ela.plot_emline(3, 1)
    ela.plot_emline(4, 1)
    ela.plot_emline(5, 1)

    ela.plot_emline(4, 2, color=1)
    ela.plot_emline(4, 3, color=1)
    ela.plot_emline(2, 2, color=1)
    ela.plot_emline(2, 3, color=1)
    ela.plot_emline(3, 2, color=1)
    ela.plot_emline(1, 2, color=1)
    ela.plot_emline(1, 3, color=1)
    ela.plot_emline(2, 4, color=1)
    ela.plot_emline(3, 3, color=1)
    ela.plot_emline(5, 4, color=1)

    sm.plot_peaks(peaks)
    plt.savefig(fname, dpi=800, bbox_inches="tight")
    plt.show()


# temperature dependence
if False:
    folders = ("Data/JinD 3", "Data/JinD 2", "Data/JinD 1", "Data/JinD 0")
    temps = (104, 60, 32, 13)
    spectra = [sm.read_all_2dspectra(f) for f in folders]

    ela.plot_temperature_dependence(spectra, temps)
    plt.xlim((6300, 6600))

    ela.add_empeak(1, 1)
    ela.add_empeak(2, 1)
    ela.add_empeak(3, 1)
    ela.add_empeak(4, 1)
    ela.add_empeak(5, 1)

    ela.add_empeak(4, 2)
    ela.add_empeak(4, 3)
    ela.add_empeak(2, 2)
    ela.add_empeak(2, 3)
    ela.add_empeak(3, 2)
    ela.add_empeak(1, 2)
    ela.add_empeak(1, 3)
    ela.add_empeak(2, 4)
    ela.add_empeak(3, 3)
    ela.add_empeak(5, 4)

    plt.savefig("Figures/JinD_site2_temp.png", dpi=400)
    plt.show()


def gauss(x, a, b, w, x0):
    return a * np.exp(-((x - x0) ** 2) / (2 * w ** 2)) + b


def lorentz(x, a, b, w, x0):
    return a / (1 + ((x - x0) ** 2 / w)) + b


if False:
    folders = ("Data/JinD 0", "Data/JinD 1", "Data/JinD 2", "Data/JinD 3")
    temps = (13, 32, 60, 104)
    spectra = [sm.read_all_2dspectra(f) for f in folders]
    w = 2
    d = 25
    for s in spectra:
        for r in s:
            if (min(r.ex) < ela.z1y1) and (max(r.ex) > ela.z1y1):
                goodinds = np.bitwise_and(r.ex > ela.z1y1 - w, r.ex < ela.z1y1 + w)
                spec = r.spec[:, goodinds]
                ex = r.ex[goodinds]
                emrange = 6470
                emind = int(sum(r.em < emrange))
                spec = np.sum(spec[emind - d : emind + d], axis=0)
                ex = ex - np.mean(ex)
                spec = spec - min(spec)
                p0 = (max(spec) - 100, 100, 0.3, 0)
                popt1, pcov1 = curve_fit(gauss, ex, spec, p0=p0)
                popt2, pcov2 = curve_fit(lorentz, ex, spec, p0=p0)

                print("At temperature %.1f:" % r.temp)
                print("Gaussian SSE: %.1f" % sum((gauss(ex, *popt1) - spec) ** 2))
                print("Gaussian FWHM: %.2f" % (popt1[2] * 2.355))
                print("Lorentzian SSE: %.1f" % sum((lorentz(ex, *popt2) - spec) ** 2))
                print("Lorentzian FWHM: %.2f" % (2 * np.sqrt(popt2[2])))

                plt.figure()
                plt.plot(ex, spec)
                plt.plot(ex, gauss(ex, *popt1))
                plt.plot(ex, lorentz(ex, *popt2))
                plt.legend(("Data", "Gaussian", "Lorentzian"))
                plt.show()
