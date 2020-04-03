import numpy as np
import matplotlib.pyplot as plt
import spectroscopymain as sm

# Give the path to the folder that contains the spectra you want.
folder = "Data/JinD 2"

# Decide on parameters
N = 4  # the number of smoothing iterations
Wx = 3  # size of smoothing window in x direction
Wy = 3  # size of smoothing window in y direction
filter_lw = 0.07  # minimum peak separation for initial filtration
filter_nf = 1 / 10  # fraction of noise level that peak must surpass
grouping_Wx = 0.4  # width of bins for grouping (x direction)
grouping_Wy = 1.3  # width of bins for grouping (y direction)

# Read the data in
spectra = sm.read_all_2dspectra(folder)

# Smooth the spectra for peak finding and plotting
spectra = sm.iterative_smooth(spectra, N=N, Wx=Wx, Wy=Wy)

# Identify the peaks you care about
peaks = sm.find_peaks(spectra, linewidth=filter_lw, noisefraction=filter_nf)

# Identify the lines you care about
lines = sm.find_lines(spectra, peaks, Wx=grouping_Wx, Wy=grouping_Wy)

# These lines can be saved or used in later processing
(exlines, emlines) = lines

# Plot the peaks and lines on top of the spectrum
sm.plot_spectra(spectra)
sm.plot_peaks(peaks)
sm.plot_lines(lines)
plt.show()
