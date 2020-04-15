import numpy as np
import matplotlib.pyplot as plt
import spectroscopymain as sm
import sys

if len(sys.argv) < 2:
    print("Which folder?")
    sys.exit()
else:
    folder = "Data/" + sys.argv[1]

peakfilename = folder + "/peaks.npy"

N = 2
Wx = 1
Wy = 1
filter_lw = 0.2
filter_nf = 1 / 15

spectra = sm.read_all_2dspectra(folder)
spectra = sm.iterative_smooth(spectra, N=N, Wx=Wx, Wy=Wy)

peaks = sm.find_peaks(spectra, linewidth=filter_lw, noisefraction=filter_nf)
np.save(peakfilename, peaks)
