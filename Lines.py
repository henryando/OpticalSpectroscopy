import numpy as np
import matplotlib.pyplot as plt
import spectroscopymain as sm
import sys
from lineanalysis import EnergyLevelsAssignments

site = 0
grouping_Wx = 0.5
grouping_Wy = 1
if len(sys.argv) < 2:
    print("Which folder?")
    sys.exit()
else:
    folder = "Data/" + sys.argv[1]
    print(folder)
if len(sys.argv) == 3:
    site = int(sys.argv[2])
if len(sys.argv) == 4:
    grouping_Wx = float(sys.argv[2])
    grouping_Wy = float(sys.argv[3])
if len(sys.argv) == 5:
    site = int(sys.argv[2])
    grouping_Wx = float(sys.argv[3])
    grouping_Wy = float(sys.argv[4])


peakfilename = folder + "/peaks.npy"
linefilename = folder + "/lines%d.npy" % site
N = 2
Wx = 1
Wy = 2


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
