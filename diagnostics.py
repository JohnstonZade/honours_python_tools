'''
Code to calculate diagnostics for fluid simulations.
'''
import glob
import os
import pickle
import numpy as np
import numpy.fft as fft
import matplotlib.pyplot as plt
from pathlib import Path
from athena_read import athdf, athinput, hst
from math import pi
from matplotlib import rc
rc('text', usetex=True)  # LaTeX labels

# Path to simulation data
PATH = '/media/zade/Seagate Expansion Drive/honours_project_2020/'
# PATH = '/media/zade/STRONTIUM/honours_project_2020/'
REGIMES = ['collisionless', 'braginskii with heat fluxes', 'braginskii']
DEFAULT_PROB = 'shear_alfven'


def load_data(output_dir, n, prob=DEFAULT_PROB):
    '''Loads data from .athdf files output from Athena++.
    '''

    def f(n):
        return folder + '.out' + output_id + '.%05d' % n + '.athdf'

    # Input
    folder = PATH + output_dir + '/' + prob  # Name of output
    output_id = '2'  # Output ID (set in input file)
    filename = f(n)
    return athdf(filename)


def load_hst(output_dir, prob=DEFAULT_PROB):
    '''Loads data from .hst files output from Athena++.
    '''
    hstLoc = PATH + output_dir + '/' + prob + '.hst'
    return hst(hstLoc)


def load_athinput(athinput_path):
    '''Loads data from athinput files.
    '''
    return athinput(PATH + athinput_path)


def load_dict(output, fname=''):
    file = 'dict.pkl'
    if fname != '':
        file = fname + '_' + file

    pkl_file = open(PATH+output+file, 'rb')
    dict = pickle.load(pkl_file)
    pkl_file.close()
    return dict


def save_dict(dict, output, fname=''):
    file = 'dict.pkl'
    if fname != '':
        file = fname + '_' + file

    output = open(PATH+output+file, 'wb')
    pickle.dump(dict, output)
    output.close()


def check_dict(output, fname=''):
    file = 'dict.pkl'
    if fname != '':
        file = fname + '_' + file
    return os.path.isfile(PATH+output+file)


def make_folder(fname):
    if not os.path.exists(fname):
        path = Path(fname)
        path.mkdir(parents=True, exist_ok=True)


def get_maxn(fname):
    '''Gets the total number of simulation timesteps.
    '''
    return len(glob.glob(PATH+fname+'/*.athdf'))


# --- VECTOR FUNCTIONS --- #


def get_mag(X):
    '''Returns the magnitude of the given vector.'''
    x = np.array(X)
    return np.sqrt(np.sum(x**2, axis=1))


def get_unit(v):
    '''Calculates unit vector.'''
    v_mag = get_mag(v)
    return np.array([v[i]/v_mag[i] for i in range(len(v))])


def get_vec(v, p):
    '''Returns the vector components at a given point.'''
    tp = tuple(p)
    v1 = v[0][tp]  # x-component
    v2 = v[1][tp]  # y-component
    v3 = v[2][tp]  # z-component
    return np.array([v1, v2, v3])


def get_vol(fname, prob=DEFAULT_PROB):
    '''Returns the volume of the simulation domain.'''
    data = load_data(fname, 0, prob)
    X1 = data['RootGridX1'][1] - data['RootGridX1'][0]
    X2 = data['RootGridX2'][1] - data['RootGridX2'][0]
    X3 = data['RootGridX3'][1] - data['RootGridX3'][0]
    return abs(X1*X2*X3)  # just a check to make volume positive


# --- Δp STATISTICS --- #


def get_time_vec(fname, prob=DEFAULT_PROB):
    max_n = get_maxn(fname)
    hst_data = load_hst(fname, prob=prob)
    t = hst_data['time']
    dt = t[1] - t[0]  # hst dt
    dt_big = t[-1]/(max_n-1)  # actual data dt
    step = int(round(dt_big / dt))
    t = t[::step]
    return np.array(t)


def box_averaged_dp(fname, prob=DEFAULT_PROB):
    '''
    Averages the quantity Δp over the simulation domain.
    '''
    max_n = get_maxn(fname)
    t = get_time_vec(fname, prob=prob)
    dp, dp_std = [], []

    for n in range(max_n):
        data = load_data(fname, n, prob=prob)

        dp_data = (data['pprp'] - data['pprl'])
        dp_avg = np.mean(dp_data)  # box averaged dp
        dp_std.append(np.std(dp_data))

        B2 = data['Bcc1']**2 + data['Bcc2']**2 + data['Bcc3']**2
        B2_avg = np.mean(B2)

        dp_per_B2_avg = dp_avg / B2_avg  # not different from mean(dp/B2)
        dp.append(dp_per_B2_avg)
    dp = np.array(dp)
    dp_std = np.array(dp_std)
    return t, dp, dp_std


def dp_pdf(fname, n, prob=DEFAULT_PROB, scale_by_urms=0, do_plot=0):
    '''
    Calculates the PDF of the quantity Δp over the simulation domain
    at a given time.
    '''
    data = load_data(fname, n, prob=prob)
    # Perpendicular and parallel pressures
    dp = data['pprp'] - data['pprl']
    B2 = data['Bcc1']**2 + data['Bcc2']**2 + data['Bcc3']**2
    dp_per_B2 = dp / B2

    if scale_by_urms:
        u_mag_2 = data['vel1']**2 + data['vel2']**2 + data['vel3']**2
        u_rms = np.sqrt(np.mean(u_mag_2))
        dp_per_B2 /= u_rms

    t = str(data['Time'])

    # to stop big spikes at the limiters, modified such that it takes into
    # account scaling by u_rms
    if scale_by_urms:
        range = (np.amin(dp_per_B2) + 0.1, np.amax(dp_per_B2) - 0.1)
    else:
        range = (-0.99, 0.499)
    n, bins = np.histogram(dp_per_B2, bins=400, range=range, density=True)
    if do_plot:
        # plot histogram
        plt.plot(0.5*(bins[:-1] + bins[1:]), n)
        plt.yscale('log', nonposy='clip')
        plt.xlim(-1.05, 0.55)
        plt.title(r'PDF of $4\pi\Delta p / B^2$, $t = $ ' + t)
        plt.xlabel(r'$4\pi\Delta p / B^2$')
        plt.ylabel(r'$\mathcal{P}$')
        plt.show()
    else:
        return n, bins


def get_bbgu(data):
    '''
    Calculates the quantity bb:∇u at every point within the simulation domain.
    '''
    # Uses Fourier transform equivalent ∂_i u_j ↔ k_i 𝓕[u_j]
    def grad_u(K, u, i, j):
        return np.real(fft.ifftn(K[i]*fft.fftn(u[j])))

    def bbgu_calc(b, du, i, j):
        return b[i, :]*b[j, :] * du[i, j, :]

    def bbgu(b, du):
        sums = range(3)
        bbgu_elements = np.array([[bbgu_calc(b, du, i, j) for i in sums]
                                  for j in sums])
        return sum(bbgu_elements[i, j, :] for i in sums for j in sums)

    K = ft_grid(data, 0)[::-1]  # Kx Ky Kz
    u = (data['vel1'], data['vel2'], data['vel3'])

    sums = range(3)
    du = np.array([[grad_u(K, u, i, j) for j in sums]
                   for i in sums])

    B_vec = np.array([data['Bcc1'], data['Bcc2'], data['Bcc3']])
    B_mag = np.sqrt(data['Bcc1']**2 + data['Bcc2']**2 + data['Bcc3']**2)
    b = B_vec / B_mag

    return bbgu(b, du)


def box_averaged_bbgu(fname, prob=DEFAULT_PROB):
    '''
    Averages the quantity bb:∇u over the simulation domain.
    '''
    max_n = get_maxn(fname)
    t = get_time_vec(fname, prob=prob)
    bbgu, bbgu_std = [], []

    for n in range(max_n):
        data = load_data(fname, n, prob=prob)

        bbgu_data = get_bbgu(data)
        bbgu_avg = np.mean(bbgu_data)  # box averaged dp
        bbgu_std.append(np.std(bbgu_data))

        B2 = data['Bcc1']**2 + data['Bcc2']**2 + data['Bcc3']**2
        B2_avg = np.mean(B2)

        bbgu_per_B2_avg = bbgu_avg / B2_avg  # not different from mean(dp/B2)
        bbgu.append(bbgu_per_B2_avg)
    bbgu = np.array(bbgu)
    bbgu_std = np.array(bbgu_std)
    return t, bbgu, bbgu_std


def bbgu_pdf(fname, n, prob=DEFAULT_PROB, scale_by_urms=0, do_plot=0,
             dp_scale=0,
             scale_by_B2=1,
             regime=None,
             nu_c=None,
             p0=None):
    '''
    Calculates the PDF of the quantity bb:∇u over the simulation domain
    at a given time.
    '''
    data = load_data(fname, n, prob=prob)
    bbgu = get_bbgu(data)

    if scale_by_B2:
        B2 = data['Bcc1']**2 + data['Bcc2']**2 + data['Bcc3']**2
        bbgu_per_B2 = bbgu / B2
    else:
        bbgu_per_B2 = bbgu

    if scale_by_urms:
        u_mag_2 = data['vel1']**2 + data['vel2']**2 + data['vel3']**2
        u_rms = np.sqrt(np.mean(u_mag_2))
        bbgu_per_B2 /= u_rms
    if dp_scale:
        if 'braginskii' in regime:
            bbgu_per_B2 *= p0 / nu_c
        elif regime == 'collisionless':
            bbgu_per_B2 *= (data['pprp'] + 2*data['pprl'])/(2*np.pi)

    t = str(data['Time'])

    if dp_scale:
        range = (-0.99, 0.49)
    else:
        range = (-15.0, 15.0)
    n, bins = np.histogram(bbgu_per_B2, bins=400, range=range, density=True)
    if do_plot:
        # plot histogram
        plt.plot(0.5*(bins[:-1] + bins[1:]), n)
        plt.yscale('log', nonposy='clip')
        plt.xlim(-1.05, 0.55)
        bbgu_string = r'$4\pi\hat{\mathbf{b}}\hat{\mathbf{b}}:\nabla\mathbf{u}/B^2$'
        plt.title(r'PDF of ' + bbgu_string + ', $t = $ ' + t)
        plt.xlabel(bbgu_string)
        plt.ylabel(r'$\mathcal{P}$')
        plt.show()
    else:
        return n, bins


def avg_B(fname, n, background, prob=DEFAULT_PROB):
    data = load_data(fname, n, prob)
    return avg_B_data(data, background)


def avg_B_data(data, background):
    B_mag = np.sqrt(data['Bcc1']**2 + data['Bcc2']**2 + data['Bcc3']**2)
    B_avg = np.mean(B_mag)
    return B_avg


def get_B2(fname, n, prob=DEFAULT_PROB):
    data = load_data(fname, n, prob=prob)
    B2 = data['Bcc1']**2 + data['Bcc2']**2 + data['Bcc3']**2
    return B2


# --- DIMENSIONLESS PARAMETER CALCULATIONS --- #


def beta(p0, B0):
    '''
    Calculates the plasma beta β = thermal pressure / magnetic pressure.
    Uses Athena++ units, where factors of √4π are absorbed into B.
    '''
    return 2*p0/(B0**2)


def db_int(nu_c, omega_A, beta):
    # Assumes v_A = 1, which is true when B0 and ρ=1.
    regime = find_regime(nu_c, omega_A, beta)
    if regime == 0:  # collisionless
        return 1/np.sqrt(beta)
    else:
        return np.sqrt(nu_c/(omega_A*beta))


def find_regime(nu_c, omega_A, beta):
    eta = nu_c / omega_A
    if eta < 1:
        return 0  # collisionless
    elif eta < np.sqrt(beta):
        return 1  # Braginskii with heat fluxes/weakly collisional
    else:
        return 2  # Braginskii


def calc_omega_A(Lx):
    # Assumes background magnetic field along the x axis and v_A = 1 (B0=1,ρ=1)
    # Then k_∥ = 2π/Lx and ω_A = k_∥v_A = 2π/Lx
    return 2*np.pi / Lx


# --- FOURIER FUNCTIONS --- #


def ft_array(N):
    '''
    For given N, returns an array conforming to FT standard:
       [0 1 2 3 ... -N/2 -N/2+1 ... -1]
    '''
    return np.concatenate((np.arange(0, N//2, 1), [-N//2],
                           np.arange(-N//2+1, 0, 1)))


def ft_grid(data, k_grid):
    '''
    Creates a grid in k-space corresponding to the real grid given in data.
    k_grid is a boolean that when True calculates a regularly spaced array
    in k-space.
    '''

    # Z, Y, X
    p = (data['x3f'], data['x2f'], data['x1f'])
    Ls = [np.max(p[0]), np.max(p[1]), np.max(p[2])]
    Ns = [len(p[0])-1, len(p[1])-1, len(p[2])-1]

    K = {}
    for k in range(3):
        K[k] = 2j*pi/Ls[k]*ft_array(Ns[k])

    # Outputs Z, Y, X
    Ks = np.meshgrid(K[0], K[1], K[2], indexing='ij')
    if k_grid:
        Ks = (Ks, np.arange(0, np.max(np.imag(K[1])), 2*pi/Ls[1]))

    return Ks
