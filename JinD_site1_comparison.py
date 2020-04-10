import numpy as np
import matplotlib.pyplot as plt
import spectroscopymain as sm
import lineanalysis as la

folders = ("Data/JinD 3", "Data/JinD 2", "Data/JinD 1", "Data/JinD 0")
temps = (104, 60, 32, 13)
linefilename = "Data/JinD 0/lines1.npy"
lines = np.load(linefilename, allow_pickle=True)
ela = la.EnergyLevelsAssignments(lines[0], lines[1])
spectra = [sm.read_all_2dspectra(f) for f in folders]
print(spectra)
spectra[3].reverse()

ela.plot_temperature_dependence(spectra, temps)
plt.xlim((6350, 6600))

ela.add_empeak(1, 1)
ela.add_empeak(1, 2)
ela.add_empeak(2, 1)
ela.add_empeak(2, 2)

plt.show()
