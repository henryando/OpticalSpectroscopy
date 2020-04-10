import numpy as np
import matplotlib.pyplot as plt
import spectroscopymain as sm
import lineanalysis as la

folders = ("Data/JinD 0", "Data/JinD 1", "Data/JinD 2", "Data/JinD 3")
temps = (13, 30, 61, 106)
tempindex = 1

folder = folders[tempindex]
linefilename = folders[0] + "/lines2a.npy"
peakfilename = folder + "/peaks.npy"
assignmentfilename = "lines2a.npy"
total = 3
N = 2
Wx = 1
Wy = 2


spectra = sm.read_all_2dspectra(folder)
spectra = sm.iterative_smooth(spectra, N=N, Wx=Wx, Wy=Wy)
peaks = np.load(peakfilename)
lines = np.load(linefilename, allow_pickle=True)
# lines[0][0] = lines[0][0] - 0.1
# lines[0][2] = lines[0][2] + 0.1
# lines[0][3] = lines[0][3] - 0.3

# lines[1][0] = lines[1][0] + 1.3
# lines[1][1] = lines[1][1] + 1.3
# lines[1][2] = lines[1][2] + 1
# lines[1][3] = lines[1][3] + 1.3
# lines[1][4] = lines[1][4] + 1.3
# lines[1][5] = lines[1][5] + 1.3

# lines[1] = np.asarray(lines[1])

# lines[0].append(lines[0][0] + 0.4)
# lines[0].append(lines[0][0] + 0.7)
# lines[0].append(lines[0][3] + 0.5)

# lines[1].append(lines[1][1] - 2)
# lines[0].append(lines[0][3] + 0.5)

# np.save("site2a.npy", lines)
# reallines = np.load(folders[tempindex] + "/lines%d.npy" % site, allow_pickle=True)
# lines = np.load("site2a.npy", allow_pickle=True)
# ela = la.EnergyLevelsAssignments(lines[0], lines[1])
ela = la.load_object(assignmentfilename)


sm.plot_spectra(spectra)
sm.plot_peaks(peaks, size=1)
ela.plot_exline(1, 1)
ela.plot_exline(1, 2)
ela.plot_exline(1, 3)
# ela.plot_exline(1, 4)
ela.plot_emline(1, 1)
ela.plot_emline(2, 1)
ela.plot_emline(3, 1)
ela.plot_emline(4, 1)
ela.plot_emline(5, 1)
# ela.plot_emline(6, 1)

ela.plot_emline(4, 2, color="m")
ela.plot_emline(4, 3, color="m")
ela.plot_emline(2, 2, color="m")
ela.plot_emline(2, 3, color="m")
ela.plot_emline(3, 2, color="m")
ela.plot_emline(1, 2, color="m")
plt.show()
