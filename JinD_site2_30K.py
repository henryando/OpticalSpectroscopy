import numpy as np
import matplotlib.pyplot as plt
import spectroscopymain as sm
import lineanalysis as la

folders = ("Data/JinD 0", "Data/JinD 1", "Data/JinD 2", "Data/JinD 3")
temps = (13, 30, 61, 106)
tempindex = 1
site = 2

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
lines[0][0] = lines[0][0] - 0.1
lines[0][2] = lines[0][2] + 0.1
lines[0][3] = lines[0][3] - 0.3

lines[1][0] = lines[1][0] + 0.2
lines[1][1] = lines[1][1] + 0.4
lines[1][2] = lines[1][2] + 0.2
lines[1][3] = lines[1][3] + 0.4

lines[0].append(lines[0][0] + 0.4)
lines[0].append(lines[0][0] + 0.7)
lines[0].append(lines[0][3] + 0.5)
# lines[0].append(lines[0][3] + 0.5)

print(lines)
# reallines = np.load(folders[tempindex] + "/lines%d.npy" % site, allow_pickle=True)
ela = la.EnergyLevelsAssignments(lines[0], lines[1])


sm.plot_spectra(spectra)
sm.plot_peaks(peaks, size=1)
ela.plot_exline(1, 1)
ela.plot_exline(1, 2, offset=1)
ela.plot_exline(1, 3, offset=2)
ela.plot_exline(1, 4)
ela.plot_exline(1, 5, offset=1)
ela.plot_exline(1, 6, offset=2)
ela.plot_exline(1, 7, offset=3)
ela.plot_exline(1, 8)
ela.plot_exline(1, 9, offset=1)
ela.plot_exline(1, 10, offset=2)
# ela.plot_exline(1, 11)
ela.plot_emline(1, 1)
ela.plot_emline(2, 1)
ela.plot_emline(3, 1)
ela.plot_emline(4, 1)

ela.plot_emline(1, 4, color="m")
ela.plot_emline(2, 4, color="m")
ela.plot_emline(3, 4, color="m")
ela.plot_emline(4, 4, color="m")
ela.plot_emline(2, 10, color="m")
ela.plot_emline(4, 10, color="m")
plt.show()
