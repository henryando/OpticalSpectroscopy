import numpy as np
import matplotlib.pyplot as plt
import spectroscopymain as sm
import lineanalysis as la

folders = ("Data/JinD 0", "Data/JinD 1", "Data/JinD 2", "Data/JinD 3")
temps = (13, 30, 61, 106)
tempindex = 0
site = 1

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
print(lines)
# reallines = np.load(folders[tempindex] + "/lines%d.npy" % site, allow_pickle=True)
ela = la.EnergyLevelsAssignments(lines[0], lines[1])


sm.plot_spectra(spectra)
sm.plot_peaks(peaks, size=1)
ela.plot_exline(1, 1)
ela.plot_exline(1, 2)
ela.plot_exline(1, 3)
ela.plot_emline(1, 1)
ela.plot_emline(2, 1)
plt.show()
