import numpy as np
import matplotlib.pyplot as plt
import spectroscopymain as sm


folder = "Data/JinD 0"
peakfilename = folder + "/peaks_site2.npy"
linefilename = folder + "/lines%d.npy"
figfilename = folder + "/alllines.png"
colors = ("g", "m", "y", "b")
N = 2
Wx = 1
Wy = 1
filter_lw = 0.25
filter_nf = 1 / 15

spectra = sm.read_all_2dspectra(folder)
spectra = sm.iterative_smooth(spectra, N=N, Wx=Wx, Wy=Wy)

peaks = sm.find_peaks(spectra, linewidth=filter_lw, noisefraction=filter_nf)
np.save(peakfilename, peaks)
