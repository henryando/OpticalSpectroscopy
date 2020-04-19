import numpy as np
import matplotlib.pyplot as plt
import spectroscopymain as sm
import lineanalysis as la
import fieldfitting as ff
import conversions as conv
from scipy.optimize import curve_fit
import sys

figtype = sys.argv[1]

folders = ("Data/JinD 0", "Data/JinD 1", "Data/JinD 2", "Data/JinD 3")
temps = (13, 32, 60, 106)
N = 2
Wx = 2
Wy = 2


def prep(folder):
    peaks = np.load(folder + "/peaks.npy")
    spectra = sm.read_all_2dspectra(folder)
    spectra = sm.iterative_smooth(spectra, N=N, Wx=Wx, Wy=Wy)
    return peaks, spectra


def plotspec(data, fmt):
    if fmt == "small":
        sm.plot_spectra(data, figsize=(3, 4.25))
        plt.xlabel("")
        plt.ylabel("")
    else:
        sm.plot_spectra(data, figsize=(6, 8))


def save(fname, resolution):
    if resolution == "highres":
        plt.savefig(fname, dpi=800, bbox_inches="tight")
    else:
        plt.savefig(fname, dpi=300, bbox_inches="tight")


##################################################
if figtype == "spectrum":
    temp = temps[int(sys.argv[2])]
    folder = folders[int(sys.argv[2])]
    size = sys.argv[3]
    resolution = sys.argv[4]

    peaks, spectra = prep(folder)
    plotspec(spectra, size)
    save("Figures/JinD_spectrum_%dK.png" % temp, resolution)
    plt.show()


else:
    site = int(sys.argv[2])
    ela = la.load_object("Data/JinD/site%d.pkl" % site)

    ##################################################
    if figtype == "printlevels":
        print("Z1Y1:")
        print("%.1f" % ela.z1y1)
        print("Z levels:")
        [print("%.1f" % l) for l in ela.zlevels]
        print("Y levels:")
        [print("%.1f" % l) for l in ela.ylevels]

    ##################################################
    elif figtype == "tempdependence":

        fname = "Figures/JinD_site%d_temp.png" % site
        resolution = sys.argv[3]

        if site == 1:
            folders = ("Data/JinD 2", "Data/JinD 1", "Data/JinD 0")
            temps = (60, 32, 13)
        else:
            folders = ("Data/JinD 3", "Data/JinD 2", "Data/JinD 1", "Data/JinD 0")
            temps = (106, 60, 32, 13)

        spectra = [sm.read_all_2dspectra(f) for f in folders]
        ela.plot_temperature_dependence(spectra, temps)

        if site == 1:
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

        elif site == 2:
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

        if resolution == "highres":
            plt.savefig(fname, dpi=800, bbox_inches="tight")
        else:
            plt.savefig(fname, dpi=300, bbox_inches="tight")

        plt.show()

    ##################################################
    elif figtype == "plotlevels":

        temp = int(sys.argv[3])
        folder = folders[temp]
        temp = temps[temp]
        size = sys.argv[4]
        resolution = sys.argv[5]

        peaks, spectra = prep(folder)
        elp = la.EnergyLevelsPlot(ela, spectra)
        elp.generate_figure(size)

        if site == 1:
            if temp == 13:
                elp.plot_exline(1, 1)
                elp.plot_exline(1, 2, offset=1)
                elp.plot_exline(1, 3)
                elp.plot_exline(1, 4)
                elp.plot_exline(1, 5)
                elp.plot_emline(1, 1)
                elp.plot_emline(2, 1)
                elp.plot_emline(3, 1)
                elp.plot_emline(4, 1)
                elp.plot_emline(5, 1, fmt=":")
            if temp == 32:
                elp.plot_exline(1, 1)
                elp.plot_exline(1, 2, offset=1)
                elp.plot_exline(1, 3)
                elp.plot_exline(1, 4)
                elp.plot_emline(1, 1)
                elp.plot_emline(2, 1)
                elp.plot_emline(3, 1)
                elp.plot_emline(4, 1)
                elp.plot_emline(5, 1, fmt=":")
                elp.plot_emline(1, 2, color=1, offset=1)
                elp.plot_emline(1, 3, color=1)
            if temp == 60:
                elp.plot_exline(1, 1)
                elp.plot_exline(1, 2, offset=1)
                elp.plot_exline(1, 3)
                elp.plot_exline(1, 4)
                elp.plot_emline(1, 1)
                elp.plot_emline(2, 1)
                elp.plot_emline(3, 1)
                elp.plot_emline(4, 1)
                elp.plot_emline(5, 1, fmt=":")
                elp.plot_emline(1, 2, color=1, offset=1)
                elp.plot_emline(1, 3, color=1)
                elp.plot_emline(2, 4, color=1, offset=1)
                elp.plot_emline(3, 3, color=1)
                elp.plot_emline(5, 3, color=1, fmt=":")
            if temp == 106:
                elp.plot_exline(1, 1)
                elp.plot_exline(1, 2, offset=1)
                elp.plot_exline(1, 3)
                elp.plot_exline(1, 4)
                elp.plot_emline(1, 1)
                elp.plot_emline(2, 1)
                elp.plot_emline(3, 1)
                elp.plot_emline(4, 1)
                elp.plot_emline(5, 1, fmt=":")
                elp.plot_emline(1, 3, color=1)
                elp.plot_emline(3, 3, color=1)
        if site == 2:
            if temp == 13:
                elp.plot_exline(1, 1)
                elp.plot_exline(1, 2)
                elp.plot_exline(1, 3)
                elp.plot_exline(1, 4)
                elp.plot_emline(1, 1)
                elp.plot_emline(2, 1)
                elp.plot_emline(3, 1)
                elp.plot_emline(4, 1)
                elp.plot_emline(5, 1)
            if temp == 32:
                elp.plot_exline(1, 1)
                elp.plot_exline(1, 2)
                elp.plot_exline(1, 3)
                elp.plot_emline(1, 1)
                elp.plot_emline(2, 1)
                elp.plot_emline(3, 1)
                elp.plot_emline(4, 1)
                elp.plot_emline(5, 1)
                elp.plot_emline(4, 2, color=1)
                elp.plot_emline(4, 3, color=1)
                elp.plot_emline(2, 2, color=1)
                elp.plot_emline(2, 3, color=1)
                elp.plot_emline(3, 2, color=1)
                elp.plot_emline(1, 2, color=1)
            if temp == 60:
                elp.plot_exline(1, 1)
                elp.plot_exline(1, 2)
                elp.plot_exline(1, 3)
                elp.plot_emline(1, 1)
                elp.plot_emline(2, 1)
                elp.plot_emline(3, 1)
                elp.plot_emline(4, 1)
                elp.plot_emline(5, 1)
                elp.plot_emline(4, 2, color=1)
                elp.plot_emline(4, 3, color=1)
                elp.plot_emline(2, 2, color=1)
                elp.plot_emline(2, 3, color=1)
                elp.plot_emline(3, 2, color=1)
                elp.plot_emline(1, 2, color=1)
                elp.plot_emline(1, 3, color=1)
                elp.plot_emline(2, 4, color=1)
            if temp == 106:
                elp.plot_exline(1, 1)
                elp.plot_exline(1, 2)
                elp.plot_exline(1, 3)
                elp.plot_emline(1, 1)
                elp.plot_emline(2, 1)
                elp.plot_emline(3, 1)
                elp.plot_emline(4, 1)
                elp.plot_emline(5, 1)
                elp.plot_emline(4, 2, color=1)
                elp.plot_emline(4, 3, color=1)
                elp.plot_emline(2, 2, color=1)
                elp.plot_emline(2, 3, color=1)
                elp.plot_emline(3, 2, color=1)
                elp.plot_emline(1, 2, color=1)
                elp.plot_emline(1, 3, color=1)
                elp.plot_emline(2, 4, color=1)
                elp.plot_emline(3, 3, color=1)
                elp.plot_emline(5, 4, color=1)

        sm.plot_peaks(peaks)
        save("Figures/JinD_site%d_%dK.png" % (site, temp), resolution)
        plt.show()


linewidthfilename = "Figures/JinD_site0_lineshape.png"
