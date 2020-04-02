import numpy as np
import matplotlib.pyplot as plt
import spectrumtools as st
import readdata as rd
from peakfindinggui import find_peaks
import peakfindinggui as pfg
import scipy.io as sio
from peakgroupinggui import find_lines


folder = "Data/JinD 3"
peakfilename = "JinD_30K_goodpeaks.mat"
N = 4
Wx = 3
Wy = 3
filter_lw = 0.07
filter_nf = 1 / 10
grouping_Wx = 0.3
grouping_Wy = 1.5

spectra = rd.read_all_2dspectra(folder)
spectra = st.iterative_smooth(spectra, N=N, Wx=Wx, Wy=Wy)
print("There are %d spectra in this folder." % len(spectra))
print("The temperature is %.1f." % spectra[0]["temp"])

if True:
    goodpeaks = find_peaks(spectra, linewidth=filter_lw, nf=filter_nf)
    # sio.savemat(peakfilename, goodpeaks)
elif False:
    goodpeaks = sio.loadmat(peakfilename)
    goodpeaks["ex"] = goodpeaks["ex"][0]
    goodpeaks["em"] = goodpeaks["em"][0]
    exlines, emlines = find_lines(spectra, goodpeaks)
    np.save("exlines2.npy", exlines)
    np.save("emlines2.npy", emlines)

exlines = np.load("exlines2.npy")
emlines = np.load("emlines2.npy")
fig, ax = pfg.plot2d(spectra)
pfg.plot_lines((exlines, emlines))
plt.show()
