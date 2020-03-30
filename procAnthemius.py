import numpy as np
import matplotlib.pyplot as plt
import spectrumtools as st
import readdata as rd

spectra = rd.read_all_2dspectra("Anthemius")
print("There are %d spectra in this folder." % len(spectra))
s = spectra[0]
s["spec"] = st.iterative_smooth(s["spec"], 2, 2)
st.plot2d_sideplots(spectra[0])
plt.show()
