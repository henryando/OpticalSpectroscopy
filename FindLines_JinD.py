import numpy as np
import matplotlib.pyplot as plt
import spectroscopymain as sm


folder = "Data/JinD 0"
i = 1
total = i

peakfilename = folder + "/peaks_site2.npy"
linefilename = folder + "/lines2a.npy"
figfilename = folder + "/alllines.png"
colors = ("g", "m", "y", "b")
N = 2
Wx = 1
Wy = 2
grouping_Wx = 0.3
grouping_Wy = 0.5

spectra = sm.read_all_2dspectra(folder)
spectra = sm.iterative_smooth(spectra, N=N, Wx=Wx, Wy=Wy)
peaks = np.load(peakfilename)


lines = sm.find_lines(spectra, peaks, Wx=grouping_Wx, Wy=grouping_Wy)
np.save(linefilename, lines)


sm.plot_spectra(spectra)
sm.plot_peaks(peaks, size=1)
lines = np.load(linefilename, allow_pickle=True)
sm.plot_lines(lines, fmt=colors[1])
plt.savefig(figfilename)
plt.show()
