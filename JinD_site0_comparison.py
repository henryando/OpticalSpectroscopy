import numpy as np
import matplotlib.pyplot as plt
import spectroscopymain as sm
import lineanalysis as la

folders = ("Data/JinD 2", "Data/JinD 1", "Data/JinD 0")
temps = (60, 32, 13)
linefilename = "Data/JinD 0/lines0.npy"
lines = np.load(linefilename, allow_pickle=True)
ela = la.EnergyLevelsAssignments(lines[0], lines[1])
spectra = [sm.read_all_2dspectra(f) for f in folders]

ela.plot_temperature_dependence(spectra, temps)

plt.xlim((6300, 6600))

ela.add_empeak(1, 1)
ela.add_empeak(2, 1)
ela.add_empeak(3, 1)
ela.add_empeak(4, 1)

ela.add_empeak(2, 4)
ela.add_empeak(1, 3)
ela.add_empeak(3, 3)
ela.add_empeak(1, 2)

plt.show()
