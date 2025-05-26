import numpy as np
import scipy.signal as sg
import matplotlib.pyplot as plt
import cmocean
from scipy.special import hermite
from math import factorial
from utils import extr2minth, extr2maxth, spec_chirp, dens_chirp, spec_double_chirp, dens_double_chirp

def Chirp(a, b, t):
    return np.exp(2 * 1j * np.pi * (a + b * t) * t)

class Signal:
    def __init__(self, par, duration=128, base=128, T0=None, T_length=16, fs=1):
        self.par = par
        self.duration = duration
        self.base = base
        self.T_length = T_length
        self.fs = fs
        if type(par["SNR"]) is list:
            self.type = "Double_chirp"
            if T0 is None:
                self.T0 = 0
        elif "k" in par:
            self.type = "Hermite"
            if T0 is None:
                self.T0 = -T_length/2
        else:
            self.type = "Single_chirp"
            if T0 is None:
                self.T0 = 0

        N = 4 * base
        win_length = N

        t = np.linspace(self.T0, self.T0 + T_length, N)

        if self.type == "Single_chirp":
            self.signal = np.sqrt(par["SNR"]) * Chirp(par["a"], par["b"], t)
        elif self.type == "Double_chirp":
            self.signal = np.sqrt(par["SNR"][0]) * Chirp(par["a"][0], par["b"], t) + np.sqrt(par["SNR"][1]) * Chirp(par["a"][1], par["b"], t)
        elif self.type == "Hermite":
            self.signal = np.sqrt(par["SNR"]) * hermite(par["k"])(t * np.sqrt(2 * np.pi)) * np.exp(-np.pi * t**2) / (np.pi**(1/4) * 2**(par["k"]/2) * np.sqrt(factorial(par["k"])))

        w = np.sqrt(N / T_length) * (np.random.randn(N) + 1j*np.random.randn(N)) / np.sqrt(2)
        self.noise = w
        g = sg.windows.gaussian(win_length, N/(T_length * np.sqrt(2*np.pi)))

        x0 = 2**(1/4) * (self.signal + w)
        f_spec, t_spec, spectro = sg.stft(x0, fs=fs, window=g, nperseg=win_length, noverlap=win_length-1, return_onesided=False)
        self.spectrogram = abs(spectro[np.argsort(f_spec),])**2
        self.fmax = max(f_spec) * N / T_length

        y_zero, x_zero = extr2minth(self.spectrogram)
        y_max, x_max = extr2maxth(self.spectrogram)
        f_spec_sorted = f_spec[np.argsort(f_spec)]
        self.zeros = np.vstack((self.T0 + self.T_length * x_zero/N, 2 * self.fmax * f_spec_sorted[y_zero])).T
        self.max = np.vstack((self.T0 + self.T_length * x_max/N, 2 * self.fmax * f_spec_sorted[y_max])).T

    @property
    def s_b(self):
        return np.sqrt(2 / (1 + 4 * self.par["b"]**2))

    @property
    def spread(self):
        return self.s_b * (self.par["a"][1] - self.par["a"][0]) / np.sqrt(2)

    @property
    def T_max(self):
        return self.T0 + self.T_length

    def Make_grid(self):
        x_grid = np.linspace(self.T0, self.T_max, self.spectrogram.shape[1])
        y_grid = np.linspace(-self.fmax, self.fmax, self.spectrogram.shape[0])
        return np.meshgrid(x_grid, y_grid)

    def Compute_density(self):
        X, Y = self.Make_grid()
        if self.type == "Double_chirp":
            r = self.s_b * (Y - 2 * self.par["b"] * X - self.par["a"][0]) / np.sqrt(2)
            s = self.s_b * (X + 2 * self.par["b"] * Y - sum(self.par["a"]) * self.par["b"]) / np.sqrt(2)

            return dens_double_chirp(r, self.par["SNR"][0], self.par["SNR"][1], self.spread, 2 * np.pi * self.spread * s, self.s_b)

        if self.type == "Hermite":
            r = np.pi * (X**2 + Y**2)
            mysterio = 1 / np.sqrt(np.pi)

            step1 = self.par["SNR"] * mysterio * r**(self.par["k"]-1) * np.exp(-r) / factorial(self.par["k"])
            return (1 + step1 * (self.par["k"]-r)**2) * np.exp(- step1 * r)

        if self.type == "Single_chirp":
            r = self.s_b * (Y - 2 * self.par["b"] * X - self.par["a"]) / np.sqrt(2)

            return dens_chirp(r, self.par["SNR"], self.s_b)

    def Compute_theo(self):
        X, Y = self.Make_grid()
        if self.type == "Double_chirp":
            r = self.s_b * (Y - 2 * self.par["b"] * X - self.par["a"][0]) / np.sqrt(2)
            s = self.s_b * (X + 2 * self.par["b"] * Y - sum(self.par["a"]) * self.par["b"]) / np.sqrt(2)

            return spec_double_chirp(r, self.par["SNR"][0], self.par["SNR"][1], self.spread, self.s_b, 2 * np.pi * self.spread * s)

        if self.type == "Hermite":
            r = np.pi * (X**2 + Y**2)

            return self.par["SNR"] * r**self.par["k"] * np.exp(-r) / factorial(self.par["k"])

        if self.type == "Single_chirp":
            r = self.s_b * (Y - 2 * self.par["b"] * X - self.par["a"]) / np.sqrt(2)

            return spec_chirp(r, self.par["SNR"], self.s_b)

    def plot(self, show_max=False, plot_type=["spectro", "intensity"]):
        fig, ax = plt.subplots(1, len(plot_type), figsize=(10, 10))
        labsize = 17

        if type(ax) is not np.ndarray:
            ax = [ax]

        for axes, _type in zip(ax, plot_type):
            extent = [self.T0, self.T_max, -self.fmax, self.fmax]

            if _type == "spectro":
                axes.imshow(np.log10(self.spectrogram), origin='lower', extent=extent, cmap=cmocean.cm.deep)
                axes.plot(self.zeros[:, 0], self.zeros[:, 1], 'o', color="white")
                if show_max:
                    axes.plot(self.max[:, 0], self.max[:, 1], 'o', color="red")

            if _type == "intensity":
                axes.imshow(np.minimum(self.Compute_density(), 5), origin='lower', extent=extent, cmap=cmocean.cm.deep)

            if _type == "theo":
                if self.type == "Hermite":
                    axes.imshow(np.maximum(np.log(1 + self.Compute_theo()),-5), origin='lower', extent=extent, cmap=cmocean.cm.deep)
                else:
                    axes.imshow(np.maximum(np.log10(self.Compute_theo()),-5), origin='lower', extent=extent, cmap=cmocean.cm.deep)

            if self.type == "Double_chirp":
                axes.plot([self.T0, self.T_max], [self.par["a"][0], self.par["a"][0] + 2 * self.par["b"] * self.T_length], color="red")
                axes.plot([self.T0, self.T_max], [self.par["a"][1], self.par["a"][1] + 2 * self.par["b"] * self.T_length], color="red")
                axes.set_xlim(3, self.T_length - 3)
                axes.set_xticks([3, 5, 7, 9, 11, 13])
                axes.set_ylim(0, 10)

            if self.type == "Single_chirp":
                axes.plot([self.T0, self.T_max], [self.par["a"], self.par["a"] + 2 * self.par["b"] * self.T_length], color="red")
                axes.set_xlim(3, self.T_length - 3)
                axes.set_xticks([3, 5, 7, 9, 11, 13])
                axes.set_ylim(0, 10)

            if self.type == "Hermite":
                axes.set_xlim(-5, 5)
                axes.set_ylim(-5, 5)

            axes.tick_params(axis='both', which='major', labelsize=labsize)
            axes.set_aspect("equal")
            return fig, ax

    def add_yellow_box(self, ax, N0):
        a1 = self.par["a"][0]
        a2 = self.par["a"][1]
        b = self.par["b"]
        aa = self.spread
        s_b = self.s_b
        x1 = np.linspace((-(a2 - a1) * b + (N0 - 1) * np.sqrt(2) / (-aa * s_b)) / (1 + 4 * b**2), (-(a2 - a1) * b + N0 * np.sqrt(2) / (-aa * s_b)) / (1 + 4 * b**2), 500)
        x2 = np.linspace((-(a2 - a1) * b + (-2 * aa * b - (N0 - 2) / aa) * np.sqrt(2) / s_b) / (1 + 4 * b**2), (-(a2 - a1) * b + (-2 * aa * b - (N0 - 1) / aa) * np.sqrt(2) / s_b) / (1 + 4 * b**2), 500)
        y1 = (-(np.linspace(0, aa, 500) + 2 * N0 * b / aa) * np.sqrt(2) / s_b + 2 * (a1 + a2) * b**2)/(1 + 4 * b**2)
        y2 = (-(np.linspace(0, aa, 500) + 2 * (N0 - 1) * b / aa) * np.sqrt(2)/s_b + 2 * (a1 + a2) * b**2)/(1 + 4 * b**2)

        ax.plot(x1, 2 * b * x1 + a2, color="yellow")
        ax.plot(x2, 2 * b * x2 + a2  -aa * np.sqrt(2) / s_b, color="yellow")
        ax.plot(-2 * b * y1 + (a1 + a2) * b - N0 * np.sqrt(2) / (aa * s_b), y1, color="yellow")
        ax.plot(-2 * b * y2 + (a1 + a2) * b - (N0 - 1) * np.sqrt(2) / (aa * s_b), y2, color="yellow")


# sig = Signal({"SNR":100, "k":1})
# sig.plot(plot_type=["theo"])
# sig.plot(plot_type=["spectro"])
# sig.plot(plot_type=["intensity"])

# sig = Signal({"SNR":100, "a":-1, "b":0.4})
# sig.plot(plot_type=["theo"])
# sig.plot(plot_type=["spectro"])
# sig.plot(plot_type=["intensity"])

# sig = Signal({"SNR":[100, 40], "a":[-1, 0], "b":0.4}, T_length=15)
# sig.plot(plot_type=["theo"])
# sig.plot(plot_type=["spectro"])
# sig.plot(plot_type=["intensity"])







