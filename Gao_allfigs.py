import numpy as np
import matplotlib.pyplot as plt
import spectroscopymain as sm
import lineanalysis as la
import fieldfitting as ff


fname = "Figures/Gao.png"
folder = "Data/Gao"
N = 1
Wx = 1
Wy = 1

lines = np.load(folder + "/lines.npy", allow_pickle=True)
peaks = np.load(folder + "/peaks.npy")

# this is where we pick the lines we really care about:
ex = lines[0][2:]
em = lines[1][:-2]

ela = la.EnergyLevelsAssignments(ex, em)
spectra = sm.read_all_2dspectra(folder)
spectra = sm.iterative_smooth(spectra, N=N, Wx=Wx, Wy=Wy)


# print level structure
if True:
    print("Z levels:")
    [print(l) for l in ela.zlevels]
    print("Y levels:")
    [print(l + ela.z1y1) for l in ela.ylevels]

# plot lines
if True:
    sm.plot_spectra(spectra)

    ela.plot_exline(1, 1)
    ela.plot_exline(1, 2)
    ela.plot_exline(1, 3)
    ela.plot_emline(1, 1)
    ela.plot_emline(2, 1)
    ela.plot_emline(3, 1)
    ela.plot_emline(4, 1)

    sm.plot_peaks(peaks)
    plt.savefig(fname, dpi=400)
    plt.show()
