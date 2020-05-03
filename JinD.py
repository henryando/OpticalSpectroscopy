import numpy as np
import matplotlib.pyplot as plt
import spectroscopymain as sm
import lineanalysis as la
import fieldfitting as ff
import conversions as conv
import linewidthtools as lwt
import sys
from scipy.optimize import curve_fit


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
    else:
        sm.plot_spectra(data, figsize=(6, 8))


def save(fname, resolution):
    if resolution == "highres":
        plt.savefig(fname, dpi=800, bbox_inches="tight")
    else:
        plt.savefig(fname, dpi=300, bbox_inches="tight")


def double_exp(x, a1, w1, a2, w2, c, b):
    return b + np.min([a1 * np.exp((x - c) / w1), a2 * np.exp((c - x) / w2)], axis=0)


def secondstark(E, a1, c1, w1, a2, c2, w2, fwhm, b):
    f1 = (E > c1) * a1 * np.exp(-E / w1)
    f2 = (E < c2) * a2 * np.exp(E / w2)
    f1 = np.convolve(f1, lwt.lorentz(E, 1, 0, fwhm, 0, 0))[
        int(len(E) / 2) : int(3 * len(E) / 2)
    ]
    f2 = np.convolve(f2, lwt.lorentz(E, 1, 0, fwhm, 0, 0))[
        int(len(E) / 2) : int(3 * len(E) / 2)
    ]
    return 10 * (f1 + f2) / len(E) + b


def site3(x, fname, resolution, dofit):
    goodinds = np.bitwise_and(excitation > x - width, excitation < x + width)
    ex = excitation[goodinds] - x
    spec = spectrum[goodinds]
    p0 = (max(spec), min(spec), 0.05, 0, 0)
    gopt, _ = curve_fit(lwt.gauss, ex, spec, p0=p0)
    lopt, _ = curve_fit(lwt.lorentz, ex, spec, p0=p0)
    gr2 = lwt.rsquared(lwt.gauss, gopt, ex, spec)
    lr2 = lwt.rsquared(lwt.lorentz, lopt, ex, spec)

    print(gr2)
    print(lr2)

    plt.figure(figsize=(3.2, 3.2))
    plt.plot(ex, spec)
    plt.xlabel("Detuning / cm$^{-1}$")
    plt.ylabel("PLE Intensity / a.u.")
    save(fname, resolution)
    plt.show()

    if dofit:
        p0 = (
            max(spec) - min(spec),
            0.05,
            0.3,
            max(spec) - min(spec),
            -0.05,
            0.3,
            0.005,
            min(spec),
        )
        eopt, _ = curve_fit(secondstark, ex, spec, p0=p0)
        plt.figure(figsize=(3.2, 3.2))
        plt.plot(ex, spec)
        plt.plot(ex, secondstark(ex, *eopt))
        plt.xlabel("Detuning / cm$^{-1}$")
        plt.ylabel("PLE Intensity / a.u.")
        plt.legend(("Observed", "Fitted"))
        save(fname, resolution)
        print(eopt)
        plt.show()


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

##################################################
if figtype == "sidebands":
    temp = temps[int(sys.argv[2])]
    folder = folders[int(sys.argv[2])]
    size = sys.argv[3]
    resolution = sys.argv[4]

    peaks, spectra = prep(folder)
    plotspec(spectra, size)
    plt.xlim((6564, 6573))
    plt.ylim((6460, 6520))
    save("Figures/JinD_sidebands_%dK.png" % temp, resolution)
    plt.show()

##################################################
elif figtype == "site3starkshift":
    resolution = sys.argv[2]
    folder = folders[2]
    data = sm.read_all_excitation(folder)[0]
    excitation = data.ex + 0.3
    spectrum = data.spec
    width = 0.7
    nofitcenter = 0.05
    ela = la.load_object("Data/JinD/site3.pkl")

    x = ela.z1y1
    site3(x, "Figures/Site3_y1.png", resolution, False)
    x = ela.z1y1 + (ela.ylevels[1] + ela.ylevels[2]) / 2
    site3(x, "Figures/Site3_y2.png", resolution, True)
    x = ela.z1y1 + ela.ylevels[3]
    site3(x, "Figures/Site3_y3.png", resolution, False)


##################################################
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
    if figtype == "fitlevels":
        if not site == 3:
            w, x, d = ff.field_fit(ela.zlevels, 15 / 2, wsign=1)
            print("B4: %.2e, B6: %.2e, dist: %.3f" % (*conv.wx15_tobs(w, x), d))
            if d < 1:
                print("looks good...trying other side...")
                w, x = conv.wx15_towx13(w, x)
                levels = ff.get_levels(w, x, 13 / 2)
                print("dist: %.3f" % ff.energy_dist(ela.ylevels, levels))
            w, x, d = ff.field_fit(ela.zlevels, 15 / 2, wsign=-1)
            print("B4: %.2e, B6: %.2e, dist: %.3f" % (*conv.wx15_tobs(w, x), d))
            if d < 1:
                print("looks good...trying other side...")
                w, x = conv.wx15_towx13(w, x)
                levels = ff.get_levels(w, x, 13 / 2)
                print("dist: %.3f" % ff.energy_dist(ela.ylevels, levels))
        w, x, d = ff.field_fit(ela.ylevels, 13 / 2, wsign=1)
        print("B4: %.2e, B6: %.2e, dist: %.3f" % (*conv.wx13_tobs(w, x), d))
        if d < 1:
            print("looks good...trying other side...")
            w, x = conv.wx13_towx15(w, x)
            levels = ff.get_levels(w, x, 15 / 2)
            print("dist: %.3f" % ff.energy_dist(ela.zlevels, levels))
        w, x, d = ff.field_fit(ela.ylevels, 13 / 2, wsign=-1)
        print("B4: %.2e, B6: %.2e, dist: %.3f" % (*conv.wx13_tobs(w, x), d))
        if d < 1:
            print("looks good...trying other side...")
            w, x = conv.wx13_towx15(w, x)
            levels = ff.get_levels(w, x, 15 / 2)
            print("dist: %.3f" % ff.energy_dist(ela.zlevels, levels))

    ##################################################
    elif figtype == "linewidths":
        temp = int(sys.argv[3])
        folder = folders[temp]
        data = sm.read_all_excitation(folder)[0]
        if (temp < 2) or (temp == 3):
            data.ex = data.ex + 1.2
        else:
            data.ex = data.ex + 0.3
        if len(sys.argv) == 5:
            width = float(sys.argv[4])
        else:
            width = 1.5
        lwt.fit_expeaks(data.ex, data.spec, ela, width=width)

    ##################################################
    elif figtype == "tempdependence":

        fname = "Figures/JinD_site%d_temp.png" % site
        resolution = sys.argv[3]

        if site == 1:
            folders = ("Data/JinD 2", "Data/JinD 1", "Data/JinD 0")
            temps = (60, 32, 13)
        elif site == 4:
            folders = ("Data/JinD 3", "Data/JinD 2", "Data/JinD 1")
            temps = (106, 60, 32)
        else:
            folders = ("Data/JinD 3", "Data/JinD 2", "Data/JinD 1", "Data/JinD 0")
            temps = (106, 60, 32, 13)

        spectra = [sm.read_all_2dspectra(f) for f in folders]
        if not site == 3:
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

        elif site == 5:
            plt.xlim()
            ela.add_empeak(1, 1)
            ela.add_empeak(2, 1)
            ela.add_empeak(3, 1)
            ela.add_empeak(1, 2)
            ela.add_empeak(2, 2)
            ela.add_empeak(3, 2)

        elif site == 4:
            plt.xlim()
            ela.add_empeak(1, 1)
            ela.add_empeak(2, 1)
            ela.add_empeak(3, 1)
            ela.add_empeak(1, 2)
            ela.add_empeak(2, 2)
            ela.add_empeak(3, 2)

        elif site == 3:
            linefilename = "Data/JinD 0/lines1.npy"
            lines = np.load(linefilename, allow_pickle=True)
            ela = la.EnergyLevelsAssignments(lines[0], lines[1])
            spectra[3].reverse()

            ela.plot_temperature_dependence(spectra, temps)
            plt.xlim((6350, 6600))

            ela.add_empeak(1, 1)
            ela.add_empeak(1, 2)
            ela.add_empeak(1, 3)

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
        if site == 3:
            if temp == 13:
                lines = np.load("Data/JinD 0/lines3.npy", allow_pickle=True)
                elp.plot_exline(1, 1)
                elp.plot_exline(1, 2)
                elp.plot_exline(1, 3, offset=1)
                elp.plot_exline(1, 4)
                elp.plot_emline(1, 1)
                elp.plot_emline(1, 1, color="y", name="?", coord=lines[1][1])
                elp.plot_emline(1, 1, color="y", name="?", coord=lines[1][0])
            if temp == 32:
                lines = np.load("Data/JinD 1/lines3.npy", allow_pickle=True)
                elp.plot_exline(1, 1)
                elp.plot_exline(1, 2)
                elp.plot_exline(1, 3, offset=1)
                elp.plot_exline(1, 4)
                elp.plot_emline(1, 1)
                elp.plot_emline(1, 2, color=1)
                elp.plot_emline(1, 1, color="y", name="?", coord=lines[1][1])
                elp.plot_emline(1, 1, color="y", name="?", coord=lines[1][0])
            if temp == 60:
                lines = np.load("Data/JinD 2/lines3.npy", allow_pickle=True)
                elp.plot_exline(1, 1)
                elp.plot_exline(1, 2)
                elp.plot_exline(1, 3, offset=1)
                elp.plot_exline(1, 4)
                elp.plot_emline(1, 1)
                elp.plot_emline(1, 2, color=1)
                elp.plot_emline(1, 4, color=1)
                elp.plot_emline(1, 1, color="y", name="?", coord=lines[1][3])
                elp.plot_emline(1, 1, color="y", name="?", coord=lines[1][2])
                elp.plot_emline(1, 1, color="y", name="?", coord=lines[1][1])
                elp.plot_emline(1, 1, color="y", name="?", coord=lines[1][0])
            if temp == 106:
                lines = np.load("Data/JinD 3/lines3.npy", allow_pickle=True)
                elp.plot_exline(1, 1)
                elp.plot_exline(1, 2)
                elp.plot_exline(1, 3, offset=1)
                elp.plot_exline(1, 4)
                elp.plot_emline(1, 1)
                elp.plot_emline(1, 2, color=1)
                elp.plot_emline(1, 4, color=1)
                elp.plot_emline(1, 1, color="y", name="?", coord=lines[1][3])
                elp.plot_emline(1, 1, color="y", name="?", coord=lines[1][2])
                elp.plot_emline(1, 1, color="y", name="?", coord=lines[1][1])
                elp.plot_emline(1, 1, color="y", name="?", coord=lines[1][0])
        if site == 4:
            if temp == 13:
                elp.plot_exline(1, 2)
                elp.plot_exline(1, 1, fmt="--")
                elp.plot_emline(1, 1)
                elp.plot_emline(2, 1)
                elp.plot_emline(3, 1)
            if temp == 32:
                elp.plot_exline(1, 1)
                elp.plot_exline(1, 2)
                elp.plot_emline(1, 1)
                elp.plot_emline(2, 1)
                elp.plot_emline(3, 1)
                elp.plot_emline(1, 2, color=1)
                elp.plot_emline(2, 2, color=1)
            if temp == 60:
                elp.plot_exline(1, 1)
                elp.plot_exline(1, 2)
                elp.plot_emline(1, 1)
                elp.plot_emline(2, 1)
                elp.plot_emline(3, 1)
                elp.plot_emline(1, 2, color=1)
                elp.plot_emline(2, 2, color=1)
                elp.plot_emline(3, 2, color=1)
            if temp == 106:
                elp.plot_exline(1, 1)
                elp.plot_exline(1, 2)
                elp.plot_emline(1, 1)
                elp.plot_emline(2, 1)
                elp.plot_emline(3, 1)
                elp.plot_emline(1, 2, color=1)
                elp.plot_emline(2, 2, color=1)

        if site == 5:
            if temp == 13:
                elp.plot_exline(1, 1)
                elp.plot_exline(1, 2)
                elp.plot_emline(1, 1)
                elp.plot_emline(2, 1)
                elp.plot_emline(3, 1)
            if temp == 32:
                elp.plot_exline(1, 1)
                elp.plot_exline(1, 2)
                elp.plot_emline(1, 1)
                elp.plot_emline(2, 1)
                elp.plot_emline(3, 1)
                elp.plot_emline(1, 2, color=1)
                elp.plot_emline(3, 2, color=1)
            if temp == 60:
                elp.plot_exline(1, 1)
                elp.plot_exline(1, 2)
                elp.plot_emline(1, 1)
                elp.plot_emline(2, 1)
                elp.plot_emline(3, 1)
                elp.plot_emline(1, 2, color=1)
                elp.plot_emline(2, 2, color=1)
                elp.plot_emline(3, 2, color=1)
            if temp == 106:
                elp.plot_exline(1, 1)
                elp.plot_exline(1, 2)
                elp.plot_emline(1, 1)
                elp.plot_emline(2, 1)
                elp.plot_emline(3, 1)
                elp.plot_emline(1, 2, color=1)
                elp.plot_emline(2, 2, color=1)

        if temp == 13:
            plt.text(6470, 6700, "a)", color="k", fontsize=14)
        if temp == 32:
            plt.text(6470, 6700, "b)", color="k", fontsize=14)
        if temp == 60:
            plt.text(6470, 6700, "c)", color="k", fontsize=14)
        if temp == 106:
            plt.text(6470, 6700, "d)", color="k", fontsize=14)

        sm.plot_peaks(peaks)
        save("Figures/JinD_site%d_%dK.png" % (site, temp), resolution)
        plt.show()


linewidthfilename = "Figures/JinD_site0_lineshape.png"
