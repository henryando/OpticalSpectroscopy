import numpy as np
import matplotlib.pyplot as plt
import spectroscopymain as sm
import lineanalysis as la


fname = "Figures/Gao.png"
folder = "Data/Gao"
N = 1
Wx = 1
Wy = 1

lines = np.load(folder + "/lines.npy", allow_pickle=True)
peaks = np.load(folder + "/peaks.npy")
ela = la.load_object(folder + "/levels.pkl")
ela.ylevels[2] = ela.ylevels[2] + 0.2

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
    ela.plot_exline(1, 2, offset=1)
    ela.plot_exline(1, 3)

    ela.plot_emline(1, 1)
    ela.plot_emline(2, 1)
    ela.plot_emline(3, 1)
    ela.plot_emline(4, 1)
    ela.plot_emline(5, 1)

    ela.plot_exline(2, 1, color="m")
    ela.plot_exline(3, 3, color="m")

    ela.plot_emline(1, 2, color="m", offset=1)
    ela.plot_emline(1, 3, color="m")
    ela.plot_emline(3, 3, color="m", offset=1)
    ela.plot_emline(5, 3, color="m")

    sm.plot_peaks(peaks)
    plt.savefig(fname, dpi=400)
    plt.show()
