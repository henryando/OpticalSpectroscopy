import numpy as np
import matplotlib.pyplot as plt
import spectroscopymain as sm
import lineanalysis as la

folders = ("Data/JinD 0", "Data/JinD 1", "Data/JinD 2", "Data/JinD 3")
temps = (13, 32, 61, 106)
tempindex = 2
site = 0

folder = folders[tempindex]
linefilename = folders[0] + "/lines%d.npy"
peakfilename = folder + "/peaks.npy"
total = 3
N = 2
Wx = 1
Wy = 2


spectra = sm.read_all_2dspectra(folder)
spectra = sm.iterative_smooth(spectra, N=N, Wx=Wx, Wy=Wy)
peaks = np.load(peakfilename)
lines = np.load(linefilename % site, allow_pickle=True)
reallines = np.load(folders[tempindex] + "/lines%d.npy" % site, allow_pickle=True)
newline = min(reallines[1])
ela = la.EnergyLevelsAssignments(lines[0], lines[1])
ela.ylevels = ela.ylevels[[0, 1, 2, 4, 5]]
newline = ela.y1z1 - newline + ela.ylevels[2]
ela.zlevels = np.concatenate((ela.zlevels, [newline]), axis=0)
ela.save_object("site0.npy")


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

plt.show()
