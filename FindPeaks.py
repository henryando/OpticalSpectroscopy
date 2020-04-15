import numpy as np
import matplotlib.pyplot as plt
import spectroscopymain as sm

folders = [
    "Anthemius",  # which is the good data here?
    "Gao",
    "GlyA2",
    "GlyceriusA2",
    "Jie1",
    "Jie2",
    "Kong",
    "Leo",
    "Minnesota",
    "OlybriusA",
    "PTO_DSO",
    "PTO_DSO_10_1",
    "PTO_DSO_25_1",
    "PTO_STO",
    "PTO_STO_0p1",
    "PTO_STO_25_1",
    "Petronius",
    "Petronius 2",
    "PetroniusYVO",
    "Tang",
    "Theo",
    "Theodesius",
    "Wai",
    "Wo",
]

folder = "Data/" + folders[1]
peakfilename = folder + "/peaks.npy"
linefilename = folder + "/lines.npy"

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
