import numpy as np
import matplotlib.pyplot as plt
import spectroscopymain as sm

# Give the path to the folder that contains the spectra you want.
folder = "Data/JinD 3"

# Read the data in
spectra = sm.read_all_2dspectra(folder)

# Smooth the spectra for peak finding and plotting
spectra = sm.iterative_smooth(spectra)

# Identify the peaks you care about
peaks = sm.find_peaks(spectra)

# Identify the lines you care about
lines = sm.find_lines(spectra, peaks)

# These lines can be saved or used in later processing
(exlines, emlines) = lines

# Plot the peaks and lines on top of the spectrum
sm.plot2d(spectra)
sm.plot_good_peaks(peaks)
sm.plot_lines(lines)
plt.show()
