import matplotlib.pyplot as plt
import spectroscopymain as sm
import lineanalysis as la
import numpy as np

folders = ("Data/JinD 0", "Data/JinD 1", "Data/JinD 2", "Data/JinD 3")
temps = (13, 30, 61, 106)
tempindex = 3

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
ela = la.load_object(assignmentfilename)


sm.plot_spectra(spectra)
sm.plot_peaks(peaks, size=1)
ela.plot_exline(1, 1)
ela.plot_exline(1, 2)
ela.plot_exline(1, 3)
ela.plot_emline(1, 1)
ela.plot_emline(2, 1)
ela.plot_emline(3, 1)
ela.plot_emline(4, 1)
ela.plot_emline(5, 1)

ela.plot_emline(4, 2, color="m")
ela.plot_emline(4, 3, color="m")
ela.plot_emline(2, 2, color="m")
ela.plot_emline(2, 3, color="m")
ela.plot_emline(3, 2, color="m")
ela.plot_emline(1, 2, color="m")
ela.plot_emline(1, 3, color="m")
ela.plot_emline(2, 4, color="m")
ela.plot_emline(3, 3, color="m")
ela.plot_emline(5, 4, color="m")
plt.show()
