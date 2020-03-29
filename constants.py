import numpy as np

# The emission axis the spectrometer returns:
CALIBRATED_EMISSION = (
    np.linspace(1, 1024, 1024) + 1.46581267e04
) / 9.91909400e00 + 20.0

# The readout noise counts we expect from each pixel of the array detector
READOUT_NOISE_COUNTS = 700
