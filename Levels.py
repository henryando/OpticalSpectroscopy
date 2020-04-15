import numpy as np
import matplotlib.pyplot as plt
import spectroscopymain as sm
import sys
from lineanalysis import EnergyLevelsAssignments

if len(sys.argv) < 2:
    print("Which folder?")
    sys.exit()
else:
    folder = "Data/" + sys.argv[1]
    print(folder)

peakfilename = folder + "/peaks.npy"
linefilename = folder + "/lines.npy"
levelname = folder + "/levels.pkl"
N = 2
Wx = 1
Wy = 2

if len(sys.argv) == 4:
    grouping_Wx = sys.argv[2]
    grouping_Wy = sys.argv[3]
else:
    grouping_Wx = 1
    grouping_Wy = 2.5


spectra = sm.read_all_2dspectra(folder)
spectra = sm.iterative_smooth(spectra, N=N, Wx=Wx, Wy=Wy)
peaks = np.load(peakfilename)
lines = np.load(linefilename, allow_pickle=True)

gslines = sm.find_lines(spectra, peaks, Wx=grouping_Wx, Wy=grouping_Wy, oldlines=lines)
ela = EnergyLevelsAssignments(gslines[0], gslines[1])
ela.save_object(levelname)

sm.plot_spectra(spectra)
sm.plot_peaks(peaks)
sm.plot_lines(lines, fmt="m")
sm.plot_lines(gslines, fmt="g")
plt.show()
