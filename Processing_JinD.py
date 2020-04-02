import numpy as np
import matplotlib.pyplot as plt
import spectroscopymain as sm


folder = "Data/JinD 3"
peakfilename = "Data/JinD 3/peaks.npy"
linefilename = "Data/JinD 3/lines.npy"
N = 4
Wx = 3
Wy = 3
filter_lw = 0.07
filter_nf = 1 / 10
grouping_Wx = 0.4
grouping_Wy = 1.3

spectra = sm.read_all_2dspectra(folder)
spectra = sm.iterative_smooth(spectra, N=N, Wx=Wx, Wy=Wy)

# peaks = sm.find_peaks(spectra, linewidth=filter_lw, noisefraction=filter_nf)
# np.save(peakfilename, (peaks["ex"], peaks["em"]))
peaks = np.load(peakfilename)

# lines = sm.find_lines(spectra, peaks)
# np.save(linefilename, lines)
lines = np.load(linefilename)

sm.plot_spectra(spectra)
sm.plot_peaks(peaks)
sm.plot_lines(lines)
plt.show()
