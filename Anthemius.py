import numpy as np
import matplotlib.pyplot as plt
import spectroscopymain as sm
import lineanalysis as la


folder = "Data/Anthemius"
N = 1
Wx = 1
Wy = 1

spectra = sm.read_all_2dspectra(folder)
spectra = sm.iterative_smooth(spectra, N=N, Wx=Wx, Wy=Wy)

# site 1
lines = np.load(folder + "/lines1.npy", allow_pickle=True)
peaks = np.load(folder + "/peaks.npy")
ela = la.load_object(folder + "/levels1.pkl")
fname = "Figures/Anthemius1.png"


# correct emission axis
peaks[1] = peaks[1] + ela.z1y1 - ela.y1z1
for s in spectra:
    s.em = s.em + ela.z1y1 - ela.y1z1


# print level structure
if False:
    print("Z1Y1:")
    print(ela.z1y1)
    print("Z levels:")
    [print(l) for l in ela.zlevels]
    print("Y levels:")
    [print(l) for l in ela.ylevels]

# plot lines
if False:
    sm.plot_spectra(spectra)

    plt.ylim((6425, 6525))

    ela.plot_exline(1, 1)
    ela.plot_exline(1, 2)
    ela.plot_exline(1, 3)
    ela.plot_exline(1, 4)
    ela.plot_exline(1, 5)

    ela.plot_emline(1, 1)
    ela.plot_emline(2, 1)
    ela.plot_emline(3, 1)
    ela.plot_emline(4, 1)

    ela.plot_exline(3, 1, color="m")
    ela.plot_emline(1, 2, color="m")
    ela.plot_emline(1, 3, color="m")

    sm.plot_peaks(peaks)
    plt.savefig(fname, dpi=400)
    plt.show()


spectra = sm.read_all_2dspectra(folder)
spectra = sm.iterative_smooth(spectra, N=N, Wx=Wx, Wy=Wy)
# site 2
lines = np.load(folder + "/lines2.npy", allow_pickle=True)
peaks = np.load(folder + "/peaks.npy")
fname = "Figures/Anthemius2.png"
lines[0] = np.sort(lines[0])
lines[1] = np.sort(lines[1])
ela = la.EnergyLevelsAssignments(lines[0][[1, 3, 4]], lines[1][:-1])
# correct emission axis
peaks[1] = peaks[1] + ela.z1y1 - ela.y1z1
for s in spectra:
    s.em = s.em + ela.z1y1 - ela.y1z1

# print level structure
if False:
    print("Z1Y1:")
    print(ela.z1y1)
    print("Z levels:")
    [print(l) for l in ela.zlevels]
    print("Y levels:")
    [print(l) for l in ela.ylevels]

# plot lines
if True:
    sm.plot_spectra(spectra)

    plt.ylim((6425, 6525))

    ela.plot_exline(1, 1)
    ela.plot_exline(1, 2)
    ela.plot_exline(1, 3)

    ela.plot_emline(1, 1)
    ela.plot_emline(2, 1)
    ela.plot_emline(3, 1)
    ela.plot_emline(4, 1)
    ela.plot_emline(5, 1)

    ela.plot_emline(1, 2, color="m")

    sm.plot_peaks(peaks)
    plt.savefig(fname, dpi=400)
    plt.show()
