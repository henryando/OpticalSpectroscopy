import numpy as np
import matplotlib.pyplot as plt
import spectroscopymain as sm


# folder = "Data/JinD 3"
folder = "Data/JinD 3"
peakfilename = folder + "/peaks.npy"
linefilename = folder + "/lines%d.npy"
figfilename = folder + "/alllines.png"
colors = ("g", "m", "y", "b")
N = 2
Wx = 2
Wy = 2
filter_lw = 0.07
filter_nf = 1 / 6
grouping_Wx = 0.1
grouping_Wy = 1

spectra = sm.read_all_2dspectra(folder)
spectra = sm.iterative_smooth(spectra, N=N, Wx=Wx, Wy=Wy)

peaks = sm.find_peaks(spectra, linewidth=filter_lw, noisefraction=filter_nf)
# np.save(peakfilename, peaks)
peaks = np.load(peakfilename)

i = 0
lines = sm.find_lines(spectra, peaks, Wx=grouping_Wx, Wy=grouping_Wy)
# np.save(linefilename % i, lines)


sm.plot_spectra(spectra)
sm.plot_peaks(peaks, size=1)
for j in range(i + 1):
    lines = np.load(linefilename % j, allow_pickle=True)
    sm.plot_lines(lines, fmt=colors[j])
plt.savefig(figfilename)
plt.show()
