import numpy as np
import matplotlib.pyplot as plt
import spectrumtools as st
import readdata as rd
import time

spectra = rd.read_all_2dspectra("JinD 3")


print("There are %d spectra in this folder." % len(spectra))
print("The temperature is %.1f." % spectra[0]["temp"])

W = 3
N = 4
plt.ion()

fig = plt.figure(figsize=(6, 6))
ax = plt.gca()

for i in range(1):
    s = spectra[i]
    peaklist = st.find_peaks(
        s, linewidth=0.1, noisefraction=1 / 8, window=W, iterations=N
    )
    s = st.iterative_smooth(s, N, W)

    # st.plot2d(s, inputfigure=fig)
    # st.plot_peaks(ax, peaklist)
    # plt.show()

    # time.sleep(1)

    print("There are %d peaks." % peaklist["em"].size)
    st.plot2d(s, inputfigure=fig)
    peaks = st.manual_peak_sort(ax, peaklist)
    print(peaks["em"])
    plt.clear()
