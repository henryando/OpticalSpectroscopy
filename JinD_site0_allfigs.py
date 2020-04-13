import numpy as np
import matplotlib.pyplot as plt
import spectroscopymain as sm
import lineanalysis as la
import fieldfitting as ff

ela = la.load_object("Data/JinD/site0.pkl")
fileformat = "Figures/JinD_site0_%dK.png"
tempfilename = "Figures/JinD_site0_temp.png"

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
if False:
    print("Z levels:")
    [print(l) for l in ela.zlevels]
    print("Y levels:")
    [print(l + ela.z1y1) for l in ela.ylevels]


if True:
    ff.field_fit(ela.zlevels[:-1], 15 / 2)
    ff.field_fit(ela.ylevels, 13 / 2)


# 13K
if False:
    fname, peaks, spectra = prep(0)
    sm.plot_spectra(spectra)

    ela.plot_exline(1, 1)
    ela.plot_exline(1, 2, offset=1)
    ela.plot_exline(1, 3)
    ela.plot_exline(1, 4)
    ela.plot_exline(1, 5)
    ela.plot_emline(1, 1)
    ela.plot_emline(2, 1)
    ela.plot_emline(3, 1)
    ela.plot_emline(4, 1)
    ela.plot_emline(5, 1, fmt=":")

    sm.plot_peaks(peaks)
    plt.savefig(fname, dpi=400)
    plt.show()


# 30K
if False:
    fname, peaks, spectra = prep(1)
    sm.plot_spectra(spectra)

    ela.plot_exline(1, 1)
    ela.plot_exline(1, 2, offset=1)
    ela.plot_exline(1, 3)
    ela.plot_exline(1, 4)
    ela.plot_emline(1, 1)
    ela.plot_emline(2, 1)
    ela.plot_emline(3, 1)
    ela.plot_emline(4, 1)
    ela.plot_emline(5, 1, fmt=":")
    ela.plot_emline(1, 2, color=1, offset=1)
    ela.plot_emline(1, 3, color=1)

    sm.plot_peaks(peaks)
    plt.savefig(fname, dpi=400)
    plt.show()


# 60K
if False:
    fname, peaks, spectra = prep(2)
    sm.plot_spectra(spectra)

    ela.plot_exline(1, 1)
    ela.plot_exline(1, 2, offset=1)
    ela.plot_exline(1, 3)
    ela.plot_exline(1, 4)
    ela.plot_emline(1, 1)
    ela.plot_emline(2, 1)
    ela.plot_emline(3, 1)
    ela.plot_emline(4, 1)
    ela.plot_emline(5, 1, fmt=":")
    ela.plot_emline(1, 2, color=1, offset=1)
    ela.plot_emline(1, 3, color=1)
    ela.plot_emline(2, 4, color=1, offset=1)
    ela.plot_emline(3, 3, color=1)
    ela.plot_emline(5, 3, color=1, fmt=":")

    sm.plot_peaks(peaks)
    plt.savefig(fname, dpi=400)
    plt.show()


# 100K
if False:
    fname, peaks, spectra = prep(3)
    sm.plot_spectra(spectra)

    ela.plot_exline(1, 1)
    ela.plot_exline(1, 2, offset=1)
    ela.plot_exline(1, 3)
    ela.plot_exline(1, 4)
    ela.plot_emline(1, 1)
    ela.plot_emline(2, 1)
    ela.plot_emline(3, 1)
    ela.plot_emline(4, 1)
    ela.plot_emline(5, 1, fmt=":")
    ela.plot_emline(1, 3, color=1)
    ela.plot_emline(3, 3, color=1)

    sm.plot_peaks(peaks)
    plt.savefig(fname, dpi=400)
    plt.show()


# temperature dependence
if False:
    folders = ("Data/JinD 2", "Data/JinD 1", "Data/JinD 0")
    temps = (60, 32, 13)
    spectra = [sm.read_all_2dspectra(f) for f in folders]

    ela.plot_temperature_dependence(spectra, temps)

    plt.xlim((6275, 6600))

    ela.add_empeak(1, 1)
    ela.add_empeak(2, 1)
    ela.add_empeak(3, 1)
    ela.add_empeak(4, 1)
    ela.add_empeak(5, 1)

    ela.add_empeak(2, 4)
    ela.add_empeak(1, 3)
    ela.add_empeak(3, 3)
    ela.add_empeak(1, 2)
    ela.add_empeak(5, 3)

    plt.savefig(tempfilename, dpi=400)
    plt.show()
