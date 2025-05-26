import numpy as np
import matplotlib.pyplot as plt
import cmocean
import scipy.signal as sg
from scipy.integrate import cumtrapz
from scipy.interpolate import interp1d
from scipy.spatial import distance_matrix

def dens_chirp(x, SNR, s_b):
    return ( 1 + 4 * SNR * np.pi * s_b * x**2 * np.exp(-2*np.pi*x**2) ) * np.exp( - SNR * s_b * np.exp(-2*np.pi*x**2) )

def spec_chirp(x, SNR, s_b):
    return SNR * s_b * np.exp(-2*np.pi*x**2)

def spec_double_chirp(x, SNR1, SNR2, a, s_b, theta):
    return s_b * (SNR1 * np.exp(-2 * np.pi * x**2) + SNR2 * np.exp(-2 * np.pi * (x-a)**2) + 2*np.sqrt(SNR1*SNR2)*np.exp(- np.pi * x**2 - np.pi * (x-a)**2)*np.cos(theta))

def dens_double_chirp(x, SNR1, SNR2, a, theta, s_b):
    return (1 + 4 * np.pi * s_b *(SNR1 * x**2*np.exp(-2 * np.pi * x**2) + SNR2 * (x-a)**2*np.exp(-2 * np.pi * (x-a)**2) + 2 * np.sqrt(SNR1*SNR2) * x*(x-a)*np.exp(- np.pi * x**2 - np.pi * (x-a)**2)*np.cos(theta)))*np.exp(-s_b * (SNR1 * np.exp(-2 * np.pi * x**2) + SNR2 * np.exp(-2 * np.pi * (x-a)**2) + 2*np.sqrt(SNR1*SNR2) * np.exp(- np.pi * x**2 - np.pi * (x-a)**2)*np.cos(theta)))

def extr2minth(M):
    central = M[1:-1, 1:-1]
    Bool = central <= M[2:, 1:-1]
    Bool = np.logical_and(Bool, central <= M[:-2, 1:-1])
    Bool = np.logical_and(Bool, central <= M[1:-1, 2:])
    Bool = np.logical_and(Bool, central <= M[1:-1, :-2])
    Bool = np.logical_and(Bool, central <= M[:-2, :-2])
    Bool = np.logical_and(Bool, central <= M[:-2, 2:])
    Bool = np.logical_and(Bool, central <= M[2:, 2:])
    Bool = np.logical_and(Bool, central <= M[2:, :-2])

    x, y = np.where(Bool)
    return x+1, y+1

def extr2maxth(M, Value=False):
    central = M[1:-1, 1:-1]
    Bool = central >= M[2:, 1:-1]
    Bool = np.logical_and(Bool, central >= M[:-2, 1:-1])
    Bool = np.logical_and(Bool, central >= M[1:-1, 2:])
    Bool = np.logical_and(Bool, central >= M[1:-1, :-2])
    Bool = np.logical_and(Bool, central >= M[:-2, :-2])
    Bool = np.logical_and(Bool, central >= M[:-2, 2:])
    Bool = np.logical_and(Bool, central >= M[2:, 2:])
    Bool = np.logical_and(Bool, central >= M[2:, :-2])

    x, y = np.where(Bool)

    if Value:
        val = M[x+1, y+1]
        return x+1, y+1, val
    return x+1, y+1

def pairCorrPlanarGaf(r, L):
    a = 0.5*L*r**2
    num = (np.sinh(a)**2+L**2/4*r**4)*np.cosh(a)-L*r**2*np.sinh(a)
    den = np.sinh(a)**3
    rho = num/den

    if r[0] == 0:
        rho[0] = 0
    return rho

def CrossCorrPlanarGaf(r):
    dat = np.load("Dataset/Cross-Correlation.npy")
    f = interp1d(dat[0,:],dat[1,:], bounds_error = False, fill_value = (0,1/3))
    return f(r)

def LmaxPlanarGaf(r):
    dat = np.load("Dataset/Max-LCorrelation.npy")
    f = interp1d(dat[0,:],dat[1,:])
    return f(r)

def FzeroPlanarGaf(r):
    dat = np.load("Dataset/F-zero.npy")
    f = interp1d(dat[0,:],dat[1,:])
    return f(r)

def GzeroPlanarGaf(r):
    dat = np.load("Dataset/G-zero.npy")
    f = interp1d(dat[0,:],dat[1,:])
    return f(r)

def GmaxPlanarGaf(r):
    dat = np.load("Dataset/G-max.npy")
    f = interp1d(dat[0,:],dat[1,:])
    return f(r)

def GcrossPlanarGaf(r):
    dat = np.load("Dataset/G-cross.npy")
    f = interp1d(dat[0,:],dat[1,:])
    return f(r)

def GdotmaxPlanarGaf(r):
    dat = np.load("Dataset/G-dot-max.npy")
    f = interp1d(dat[0,:],dat[1,:])
    return f(r)

def Kfunction(r, rho):
    K = np.zeros(len(rho))
    K[1:] = 2*np.pi*cumtrapz(r*rho, r)
    return K

def ginibreGaf(r, c):
    rho = 1-np.exp(-c*r**2)
    return rho

# Remi/Julien's code
def demoSpectrogramWhiteNoise(base):
    N = 4 * base
    Nfft = 2 * base

    # Noise only
    w = ( np.random.randn(N)+1j*np.random.randn(N) )/ np.sqrt(2)
    # window
    g = sg.gaussian(Nfft, np.sqrt((Nfft)/2/np.pi))
    g = g/g.sum()
    _, _, stft = sg.stft(w, window=g, nperseg=Nfft, noverlap=Nfft-1, return_onesided=False)

    Sww_t = np.abs(stft)**2
    print("STFT computed")

    tmin = base
    tmax = 3*base

    Sww = Sww_t[:, tmin:tmax+1]

    # detection
    y, x = extr2minth(Sww)   #Extraction de zéros

    u = (np.array(x))/np.sqrt(2*base)
    v = (np.array(y))/np.sqrt(2*base)

    pos_zero = np.zeros((len(x), 2))
    pos_zero[:, 0] = u
    pos_zero[:, 1] = v

    print("zeros extracted")

    y, x = extr2maxth(Sww)   #Extraction de maxima

    u = (np.array(x))/np.sqrt(2*base)
    v = (np.array(y))/np.sqrt(2*base)

    pos_max = np.zeros((len(x), 2))
    pos_max[:, 0] = u
    pos_max[:, 1] = v

    print("local maxima extracted")

    return pos_zero, pos_max, [Sww, x, y]

def Remove_Max(pos_max, base, val_th_list=[8,10,13,18,20,25,50,100,150], show=False):
    pos_max_corr = np.copy(pos_max)

    for val_th in val_th_list:

        D = distance_matrix(pos_max_corr, pos_max_corr)
        np.fill_diagonal(D, 1000000)
        Ind = np.argmin(D, axis=0)
        Val = np.min(D, axis=0)

        th = np.sqrt(val_th)/np.sqrt(2*4096)

        L1 = []
        L2 = []
        for a in range(len(Val)):
            if Val[a] <= th and a < Ind[a] :
                L1 += [a]
                L2 += [Ind[a]]

        for i in range(len(L2)):
            pos_max_corr[L2[i], :] = (pos_max_corr[L1[i], :] + pos_max_corr[L2[i], :])/2

        pos_max_corr = np.delete(pos_max_corr, L1, axis=0)
        if show:
            print(str(len(L1))+"pts deleted")

    return pos_max_corr

def LoadDataset(SNR, i, duration="full"):
    if SNR not in [0, 50, 12.5, 3, 1.5]:
        raise ValueError("SNR should be 0, 1.5, 3, 12.5 or 50.")
    if duration not in ["half", "full"]:
        raise ValueError("Duration should be half or full")

    if duration == "full":
        if SNR == 0:
            if i >= 10000:
                raise ValueError("Only 10000 datasets available for SNR=0.")
            dat = np.load("Dataset/White_Noise/WN_128_"+str(i)+".npz")
        elif SNR == 3:
            if i >= 2000:
                raise ValueError("Only 500 datasets available for SNR=3.")
            dat = np.load("Dataset/Chirp/full_128_SNR=3_"+str(i)+".npz")
        elif SNR == 12.5:
            if i >= 2000:
                raise ValueError("Only 500 datasets available for SNR=12.5.")
            dat = np.load("Dataset/Chirp/full_128_SNR=12.5_"+str(i)+".npz")
        elif SNR == 50:
            if i >= 2000:
                raise ValueError("Only 500 datasets available for SNR=50.")
            dat = np.load("Dataset/Chirp/full_128_SNR=50_"+str(i)+".npz")
        else:
            raise ValueError("No dataset for this SNR.")
    elif duration == "half":
        if SNR == 0:
            if i >= 10000:
                raise ValueError("Only 10000 datasets available for SNR=0.")
            dat = np.load("Dataset/White_Noise/WN_128_"+str(i)+".npz")
        elif SNR == 1.5:
            if i >= 2000:
                raise ValueError("Only 2000 datasets available for SNR=1.5.")
            dat = np.load("Dataset/Chirp/half_128_SNR=1.5_"+str(i)+".npz")
        elif SNR == 3:
            if i >= 2000:
                raise ValueError("Only 2000 datasets available for SNR=3.")
            dat = np.load("Dataset/Chirp/half_128_SNR=3_"+str(i)+".npz")
        elif SNR == 12.5:
            if i >= 2000:
                raise ValueError("Only 500 datasets available for SNR=12.5.")
            dat = np.load("Dataset/Chirp/half_128_SNR=12.5_"+str(i)+".npz")
        else:
            raise ValueError("No dataset for this SNR.")

    pos_zero = dat['arr_0']
    pos_max = dat['arr_1']
    return pos_zero, pos_max

def NumDataset(SNR, duration="full"):
    if SNR not in [0, 50, 12.5, 3, 1.5]:
        raise ValueError("SNR should be 0, 1.5, 3, 12.5 or 50.")
    if duration not in ["half", "full"]:
        raise ValueError("Duration should be half or full")

    if duration == "full":
        if SNR == 0:
            return 10000
        elif SNR == 3:
            return 2000
        elif SNR == 12.5:
            return 2000
        elif SNR == 50:
            return 2000
        else:
            raise ValueError("No dataset for this SNR.")
    elif duration == "half":
        if SNR == 0:
            return 10000
        elif SNR == 1.5:
            return 2000
        elif SNR == 3:
            return 2000
        elif SNR == 12.5:
            return 2000
        else:
            raise ValueError("No dataset for this SNR.")

def APF0(diag, r):
    diag0 = np.array([el[1] for el in diag if el[0] == 0])
    l = diag0[:,1]
    m = diag0[:,1]/2
    return(np.sum(l[m<r]))

def APF1(diag, r):
    diag1 = np.array([el[1] for el in diag if el[0] == 1])
    if len(diag1)==0: return(0)
    l = diag1[:,1] - diag1[:,0]
    m = (diag1[:,0] + diag1[:,1])/2
    return(np.sum(l[m<r]))

def demoSpectrogramSignal(SNR, duration, viz=False, shrink=True, cross=False):
    base = 128
    b = 150
    a = base - b

    N = 4 * base
    Nfft = 2 * base
    t = np.arange(N)

    # Noise only
    w = np.random.randn(N)
    #w = (np.random.randn(N)+1j*np.random.randn(N)) / np.sqrt(2)
    # window
    g = sg.gaussian(Nfft, np.sqrt((Nfft)/2/np.pi))
    g = g/g.sum()

    # bounds for detection (larger)
    fmin = 0
    fmax = base

    #tmin = 2*base - (base-trunc) // 2
    #tmax = 2*base + (base-trunc) // 2
    tmin = base
    tmax = 3*base

    #chirp
    duration = int(np.floor(duration))
    if duration > base:
        raise ValueError('Duration should be lesser than base')


    start_s = 2*base-duration //2
    end_s = 2*base+duration //2
    chirp = np.zeros(N)

    freq = (a + b*t[start_s:end_s]/N)*t[start_s:end_s]/N
    chirp[start_s:end_s] = sg.tukey(duration)*np.cos(2*np.pi*freq)

    x0 = np.sqrt(2*SNR)*chirp + w
    #spectro = STFT(x0, g, Nfft)
    _, _, spectro = sg.stft(x0, window=g, nperseg=Nfft, noverlap=Nfft-1, return_onesided=False)
    Sww_t = abs(spectro)**2
    #print("STFT computed")

    Sww = Sww_t[fmin:fmax+1, tmin:tmax+1]


    # detection
    y0, x0 = extr2minth(Sww)
    if shrink is True:
        # boundary conditions
        side = 110 # size of square; equivalent to trunc
        fmin_b = (max(0, (base-side)//2))
        fmax_b = (min(base, (base+side)//2))
        tmin_b = (base-side//2)
        tmax_b = (base+side//2)

        mask = (y0 > fmin_b)*(y0 < fmax_b)*(x0 > tmin_b)*(x0 < tmax_b)
        u = x0[mask]/np.sqrt(2*base)
        v = y0[mask]/np.sqrt(2*base)
    else:
        u = x0/np.sqrt(2*base)
        v = y0/np.sqrt(2*base)

    pos = np.zeros((len(u), 2))
    pos[:, 0] = u
    pos[:, 1] = v

    if cross:
        y0, x0 = extr2maxth(Sww)
        if shrink is True:
            # boundary conditions
            side = 110 # size of square; equivalent to trunc
            fmin_b = (max(0, (base-side)//2))
            fmax_b = (min(base, (base+side)//2))
            tmin_b = (base-side//2)
            tmax_b = (base+side//2)

            mask = (y0 > fmin_b)*(y0 < fmax_b)*(x0 > tmin_b)*(x0 < tmax_b)
            u = x0[mask]/np.sqrt(2*base)
            v = y0[mask]/np.sqrt(2*base)
        else:
            u = x0/np.sqrt(2*base)
            v = y0/np.sqrt(2*base)
        #u = x0/np.sqrt(2*base)
        #v = y0/np.sqrt(2*base)
    pos2 = np.zeros((len(u), 2))
    pos2[:, 0] = u
    pos2[:, 1] = v

    if viz is True:
        side = 110 # size of square; equivalent to trunc
        fmin = (max(0, (base-side)//2))/np.sqrt(2*base)
        fmax = (min(base, (base+side)//2))/np.sqrt(2*base)
        tmin = (base-side//2)/np.sqrt(2*base)
        tmax = (base+side//2)/np.sqrt(2*base)

        fig, ax = plt.subplots(figsize=(5, 5))

        ax.imshow(np.log10(Sww), origin='lower', extent=[0, (2*base)/np.sqrt(2*base), 0, (base)/np.sqrt(2*base)], cmap=cmocean.cm.deep)
        ax.scatter(pos[:, 0], pos[:, 1], color='w', s=40)
        if cross:
            ax.scatter(pos2[:, 0], pos2[:, 1], color='r', s=40)

        ax.set_xlim([tmin, tmax])
        ax.set_ylim([fmin, fmax])
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_xticks([])
        ax.set_yticks([])
        fig.tight_layout()
        fig.subplots_adjust(left=0.04, bottom=0.05)

        if cross:
            return Sww, pos, pos2, spectro, chirp
        return Sww, pos, spectro, chirp

    else:
        if cross:
            return pos, pos2
        return pos
