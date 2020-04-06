import numpy as np
import conversions as conv
import matplotlib.pyplot as plt
import spectroscopymain as sm


# folder = "Data/JinD 3"
folders = ("Data/JinD 0", "Data/JinD 3", "Data/JinD 2", "Data/JinD 1")
colors = ("g", "m", "y", "b")
N = 1
Wx = 1
Wy = 1

# site = 0
# for folder in folders[1:]:
#     peakfilename = folder + "/peaks.npy"
#     linefilename = folder + "/lines%d.npy" % site
#     figfilename = folder + "/alllines.png"

#     spectra = sm.read_all_2dspectra(folder)
#     spectra = sm.iterative_smooth(spectra, N=N, Wx=Wx, Wy=Wy)
#     peaks = np.load(peakfilename)
#     lines = np.load(linefilename, allow_pickle=True)

#     sm.plot_spectra(spectra)
#     sm.plot_peaks(peaks, size=1)
#     sm.plot_lines(lines, fmt=colors[site])

# plt.show()


# pick lines

# site 1:
lines = np.load(folders[2] + "/lines0.npy", allow_pickle=True)
emlines = np.sort(lines[1])
emlines = emlines[[0, 2, 3, 5, 6, 7]]
exlines = np.sort(lines[0])
emlines = conv.nm_to_wavenumber(emlines)
exlines = conv.nm_to_wavenumber(exlines)
np.save("Data/JinD/lines0.npy", (exlines, emlines))


# site 2:
lines = np.load(folders[2] + "/lines1.npy", allow_pickle=True)
emlines = np.sort(lines[1])
emlines = emlines[1:]
exlines = np.sort(lines[0])
emlines = conv.nm_to_wavenumber(emlines)
exlines = conv.nm_to_wavenumber(exlines)
np.save("Data/JinD/lines1.npy", (exlines, emlines))

# site 3:
lines = np.load(folders[1] + "/lines2.npy", allow_pickle=True)
exlines = np.sort(lines[0])[1:]
lines = np.load(folders[2] + "/lines2.npy", allow_pickle=True)
exlines = np.concatenate((exlines, np.sort(lines[0])[0:1]))
emlines = np.sort(lines[1])[[0, 1, 2, 3, 4, 5, 6, 8, 9, 10]]
emlines = conv.nm_to_wavenumber(emlines)
exlines = conv.nm_to_wavenumber(exlines)
np.save("Data/JinD/lines2.npy", (exlines, emlines))

# site 4:
lines = np.load(folders[3] + "/lines3.npy", allow_pickle=True)
exlines = np.sort(lines[0])
lines = np.load(folders[2] + "/lines3.npy", allow_pickle=True)
emlines = np.sort(lines[1])
emlines = conv.nm_to_wavenumber(emlines)
exlines = conv.nm_to_wavenumber(exlines)
np.save("Data/JinD/lines3.npy", (exlines, emlines))
