import numpy as np
import conversions as conv
import matplotlib.pyplot as plt
import spectroscopymain as sm

colors = ("g", "m", "y", "b")
N = 1
Wx = 1
Wy = 1

folder = "Data/JinD 3"
spectra = sm.read_all_2dspectra(folder)
spectra = sm.iterative_smooth(spectra, N=N, Wx=Wx, Wy=Wy)

for site in range(4):
    lines = "Data/JinD/lines%d.npy" % site

    lines = np.load(lines, allow_pickle=True)

    sm.plot_spectra(spectra)
    sm.plot_lines(lines, fmt="r")

plt.show()
