import numpy as np
import spectrumtools as st
import conversions as conv
import matplotlib.pyplot as plt
import pickle
import spectroscopymain as sm


class EnergyLevelsPlot:
    def __init__(self, ela, data):
        self.ela = ela
        self.data = data

    def generate_figure(self, fmt):
        if fmt == "small":
            sm.plot_spectra(self.data, figsize=(3, 4.25))
            plt.xlabel("")
            plt.ylabel("")
        elif fmt == "large":
            sm.plot_spectra(self.data, figsize=(6, 8.5))
        self.fmt = fmt

    def plot_peaks(self, peaks, color="r"):
        sm.plot_peaks(peaks, color=color)

    def plot_exline(self, zlevel, ylevel, color=0, fmt="-", offset=0):
        if self.fmt == "small":
            self.ela.plot_exline(
                zlevel, ylevel, color=color, fmt=fmt, offset=offset, fontsize=7
            )
        if self.fmt == "large":
            self.ela.plot_exline(
                zlevel, ylevel, color=color, fmt=fmt, offset=offset, fontsize=9
            )

    def plot_emline(self, zlevel, ylevel, color=0, fmt="-", offset=0):
        if self.fmt == "small":
            self.ela.plot_emline(
                zlevel, ylevel, color=color, fmt=fmt, offset=offset, fontsize=7
            )
        if self.fmt == "large":
            self.ela.plot_emline(
                zlevel, ylevel, color=color, fmt=fmt, offset=offset, fontsize=9
            )

    def savefig(self, fname, resolution):
        if resolution == "highres":
            plt.savefig(fname, dpi=800, bbox_inches="tight")
        else:
            plt.savefig(fname, dpi=300, bbox_inches="tight")


class EnergyLevelsAssignments:
    def __init__(self, exlines, emlines):
        self.z1y1 = min(exlines)
        self.y1z1 = max(emlines)
        self.zlevels = np.sort(max(emlines) - emlines)
        self.ylevels = np.sort(exlines) - min(exlines)

    def line_energy_ex(self, zlevel, ylevel):
        return self.z1y1 + self.ylevels[ylevel - 1] - self.zlevels[zlevel - 1]

    def line_energy_em(self, zlevel, ylevel):
        return self.z1y1 + self.ylevels[ylevel - 1] - self.zlevels[zlevel - 1]

    def plot_exline(self, zlevel, ylevel, color=0, fmt="-", offset=0, fontsize=9):
        if color == 0:
            color = "g"
        elif color == 1:
            color = "m"
        ylim = plt.ylim()
        # xlim = plt.xlim()
        sf = 30
        plt.plot(
            [self.line_energy_ex(zlevel, ylevel), self.line_energy_ex(zlevel, ylevel)],
            ylim,
            fmt,
            color=color,
            linewidth=1,
        )
        plt.text(
            self.line_energy_ex(zlevel, ylevel),
            ylim[1] + (0.4 + 1 * offset) * (ylim[1] - ylim[0]) / sf,
            "Z$_{%d}$Y$_{%d}$" % (zlevel, ylevel),
            fontsize=fontsize,
            color=color,
        )

    def plot_emline(self, zlevel, ylevel, color=0, fmt="-", offset=0, fontsize=9):
        if color == 0:
            color = "g"
        elif color == 1:
            color = "m"
            # ylim = plt.ylim()
        xlim = plt.xlim()
        sf = 30
        plt.plot(
            xlim,
            [self.line_energy_em(zlevel, ylevel), self.line_energy_em(zlevel, ylevel)],
            fmt,
            color=color,
            linewidth=1,
        )
        plt.text(
            xlim[1] + (0.2 + 3 * offset) * (xlim[1] - xlim[0]) / sf,
            self.line_energy_em(zlevel, ylevel),
            "Y$_{%d}$Z$_{%d}$" % (ylevel, zlevel),
            fontsize=fontsize,
            color=color,
        )

    def get_emission_spec(self, data, ex=None):
        if ex is None:
            ex = self.z1y1
        if type(data) is list:
            for s in data:
                emspec = self.get_emission_spec(s)
                if emspec is not None:
                    # print("returning em")
                    return emspec
        else:
            if (ex < min(data.ex)) or (ex > max(data.ex)):
                return None
            else:
                ind = sum(data.ex >= ex)
                d = data.spec[:, ind + 1]
                d = d - min(d)
                d = d / max(d[np.bitwise_and(data.em > ex - 5, data.em < ex + 5)])
                return d

    def plot_temperature_dependence(self, spectra, temps):
        emission = [self.get_emission_spec(s) for s in spectra]
        emaxis = spectra[0][0].em
        plt.figure(figsize=(6.5, 2.5))
        [plt.plot(emaxis, em, linewidth=1) for em in emission]
        plt.xlabel("Emission energy / cm$^{-1}$")
        plt.ylabel("Normalized emission / a.u.")
        # plt.title("Excitation at %.2f cm$^{-1}$" % self.z1y1)
        plt.legend(["%dK" % t for t in temps])
        plt.gcf().subplots_adjust(bottom=0.2)

    def add_empeak(self, zlevel, ylevel):
        plt.text(
            self.line_energy_em(zlevel, ylevel) - (plt.xlim()[1] - plt.xlim()[0]) / 150,
            1.1,
            "Y$_{%d}$Z$_{%d}$" % (ylevel, zlevel),
            fontsize=8,
            rotation=-90,
        )

    def save_object(self, fname):
        with open(fname, "wb") as output:
            pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)


def load_object(fname):
    with open(fname, "rb") as input:
        return pickle.load(input)


def get_visible_lines(exlines, emlines, temp, cutoff):
    """Return the lines where at least one point should be visible,
    where visible is defined to mean that the intensity relative to
    Z1-Y1 is more than cutoff.
    """
    z1y1 = min(exlines)
    zlevels = np.sort(max(emlines) - emlines)
    ylevels = np.sort(exlines) - min(exlines)

    numex = sum(np.exp(-zlevels / conv.kelvin_to_wavenumber(temp)) > cutoff)
    numem = sum(np.exp(-ylevels / conv.kelvin_to_wavenumber(temp)) > cutoff)

    exlines = []
    for i in range(numex):
        exlines = np.concatenate((exlines, ylevels - zlevels[i]), axis=0)

    emlines = []
    for i in range(numem):
        emlines = np.concatenate((emlines, -zlevels + ylevels[i]), axis=0)

    exlines = exlines + z1y1
    emlines = emlines + z1y1

    return exlines, emlines
