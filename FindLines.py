import numpy as np
import matplotlib.pyplot as plt
import spectroscopymain as sm

# Give the path to the folder that contains the spectra you want.
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
N = 2
Wx = 1
Wy = 2
grouping_Wx = 1
grouping_Wy = 2.5

spectra = sm.read_all_2dspectra(folder)
spectra = sm.iterative_smooth(spectra, N=N, Wx=Wx, Wy=Wy)
peaks = np.load(peakfilename)


lines = sm.find_lines(spectra, peaks, Wx=grouping_Wx, Wy=grouping_Wy)
np.save(linefilename, lines)


sm.plot_spectra(spectra)
sm.plot_peaks(peaks)
lines = np.load(linefilename, allow_pickle=True)
sm.plot_lines(lines)
plt.show()
