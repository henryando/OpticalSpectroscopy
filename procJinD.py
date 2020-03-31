import numpy as np
import matplotlib.pyplot as plt
import spectrumtools as st
import readdata as rd
from peakfindinggui import PeakClicker
import peakfindinggui as pfg
import scipy.io as sio
from grouping import PeakGrouper


def find_peaks(spectra, linewidth=0.07, nf=1 / 12):
    goodpeaks = {"ex": [], "em": []}
    for s in spectra:
        peakdict = st.find_peaks(s, linewidth=linewidth, noisefraction=nf)
        print("There are %d peaks." % peakdict["em"].size)
        pc = PeakClicker(s, peakdict)
        print("Double click the good peaks, then quit the plot.")
        plt.show()
        goodpeaks["ex"] = np.concatenate(
            (goodpeaks["ex"], pc.good_peaks()["ex"]), axis=0
        )
        goodpeaks["em"] = np.concatenate(
            (goodpeaks["em"], pc.good_peaks()["em"]), axis=0
        )
    pfg.plot2d_list(spectra)
    pfg.plot_good_peaks(goodpeaks["ex"], goodpeaks["em"])
    plt.show()
    return goodpeaks


folder = "Data/JinD 3"
peakfilename = "JinD_30K_goodpeaks.mat"
N = 3
Wx = 2
Wy = 2
filter_lw = 0.07
filter_nf = 1 / 12
grouping_Wx = 0.25
grouping_Wy = 0.5
grouping_bd = 0
grouping_nm = 1

spectra = rd.read_all_2dspectra(folder)
spectra = st.iterative_smooth(spectra, N=N, Wx=Wx, Wy=Wy)
print("There are %d spectra in this folder." % len(spectra))
print("The temperature is %.1f." % spectra[0]["temp"])

if False:
    goodpeaks = find_peaks(spectra, linewidth=filter_lw, nf=filter_nf)
    sio.savemat(peakfilename, goodpeaks)
else:
    goodpeaks = sio.loadmat(peakfilename)
    goodpeaks["ex"] = goodpeaks["ex"][0]
    goodpeaks["em"] = goodpeaks["em"][0]
    pg = PeakGrouper(spectra, goodpeaks, grouping_Wx, grouping_Wy, grouping_bd)
    pg.trim_omatrix(grouping_nm)

    pfg.plot2d_list(spectra)
    pfg.plot_peaks(goodpeaks["ex"], goodpeaks["em"])
    pg.show_groups()
    pg.show_omatrix()
    plt.show()
