import numpy as np
import matplotlib.pyplot as plt
import spectrumtools as st
import readdata as rd
from peakfindinggui import PeakClicker
import peakfindinggui as pfg
import scipy.io as sio


def find_peaks(spectra, N=4, W=3, linewidth=0.07, nf=1 / 12):
    goodpeaks = {"ex": [], "em": []}
    for i in range(len(spectra)):
        datadict = st.iterative_smooth(spectra[i], N, W)
        peakdict = st.find_peaks(datadict, linewidth=linewidth, noisefraction=nf)
        print("There are %d peaks." % peakdict["em"].size)
        pc = PeakClicker(datadict, peakdict)
        print("Double click the good peaks, then quit the plot.")
        plt.show()
        goodpeaks["ex"] = np.concatenate((goodpeaks["ex"], pc.good_peaks()["ex"]))
        goodpeaks["em"] = np.concatenate((goodpeaks["em"], pc.good_peaks()["em"]))
        spectra[i] = datadict
    pfg.plot2d_list(spectra)
    pfg.plot_good_peaks(goodpeaks["ex"], goodpeaks["em"])
    plt.show()
    return goodpeaks


folder = "Data/JinD 3"
peakfilename = "JinD_30K_goodpeaks.mat"
spectra = rd.read_all_2dspectra("Data/JinD 3")
print("There are %d spectra in this folder." % len(spectra))
print("The temperature is %.1f." % spectra[0]["temp"])

if False:
    goodpeaks = find_peaks(spectra)
    sio.savemat(peakfilename, goodpeaks)
else:
    goodpeaks = sio.loadmat(peakfilename)
