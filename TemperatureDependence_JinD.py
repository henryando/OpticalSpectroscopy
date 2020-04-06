import numpy as np
import matplotlib.pyplot as plt
import spectroscopymain as sm
from mpl_toolkits import mplot3d

data_folders = ["Data/JinD %d" % i for i in (0, 3, 2, 1)]
spectra = [sm.read_all_2dspectra(f) for f in data_folders]
(exlines, emlines) = np.load("Data/JinD/lines0.npy", allow_pickle=True)
temps = [s[0].temp for s in spectra]
# peaks = [sm.peak_heights_at_lines(s, exlines, emlines) for s in spectra]
emspec = [sm.ex_line_nearest_peak(s, exlines[4]) for s in spectra]
em = spectra[0][0].em
print(emspec)
# ind = sum(em < emlines[0]) + 1
# emspec = [(s - min(s)) for s in emspec]
# emspec = [s / max(s[em < 1530]) for s in emspec]
# peaks = [p / p[0, 0] for p in peaks]

plt.figure()
# ax = plt.axes(projection="3d")
# xaxis = np.ones_like(peaks[0])
# yaxis = np.ones_like(peaks[0])
# xaxis = [exlines * xaxis[:, i] for i in range(xaxis.shape[1])]
# yaxis = [emlines * yaxis[i] for i in range(yaxis.shape[0])]
# for i in range(4):
#     ax.scatter3D(xaxis, yaxis, peaks[i])
for i in range(3):
    plt.plot(em, emspec[i])
plt.legend(temps)
plt.show()
