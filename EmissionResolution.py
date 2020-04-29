import numpy as np
import matplotlib.pyplot as plt
import spectroscopymain as sm
from scipy.optimize import curve_fit


folder = "Data/JinD 0"
spectra = sm.read_all_2dspectra(folder)

# sm.iterative_smooth(spectra, N=1, Wx=1, Wy=1)
# sm.plot_spectra(spectra)
# plt.show()


def gauss(x, w, a, c, b):
    return b + a * np.exp(-((x - c) ** 2) / 2 * w ** 2)


def fitwidth(x, y):
    popt, _ = curve_fit(gauss, x, y, p0=(2, 100, np.mean(x), 100))
    return popt


testpoints = (6606, 6565, 6537, 6496)
emspectra = [sm.ex_line_nearest_peak(spectra, t) for t in testpoints]
em = spectra[0].em
masks = [np.bitwise_and(em > t - 5, em < t + 8) for t in testpoints]
plotdata = [(em[masks[i]], emspectra[i][masks[i]]) for i in range(4)]
fits = [fitwidth(pd[0], pd[1]) for pd in plotdata]

print([f[0] for f in fits])
print(np.mean([f[0] for f in fits]))
print(np.std([f[0] for f in fits]))

plt.figure(figsize=(5, 2))
[plt.plot(pd[0], pd[1], "k") for pd in plotdata]
[plt.plot(plotdata[i][0], gauss(plotdata[i][0], *fits[i]), "r--") for i in range(4)]
plt.xlabel("Emission energy / cm$^{-1}$")
plt.ylabel("Counts per second")
plt.savefig("Figures/EmissionResolution.png", dpi=300, bbox_inches="tight")
plt.show()
