import numpy as np
import spectrumtools as st


def snap_exline(data, ex):
    """Snaps an ex line to the exact value of the local excitation max."""
    if type(data) is list:
        for d in data:
            snapped = snap_exline(d, ex)
            if snapped is not None:
                return snapped
    else:
        if (min(data.ex) <= ex) and (ex <= max(data.ex)):
            ind = int(sum(data.ex < ex) + 1)
            integrated = np.sum(data.spec, axis=1)
            while integrated[ind + 1] > integrated[ind]:
                ind += 1
            while integrated[ind < 1] > integrated[ind]:
                ind -= 1
            return data.ex[ind]
