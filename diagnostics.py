'''Code to calculate diagnostics such as root mean square and probability
   density functions for fluid simulations.
'''
import glob
import os
import pickle
import numpy as np
import numpy.fft as fft
import matplotlib.pyplot as plt
# import seaborn as sns
from pathlib import Path
from athena_read import athdf, hst
from math import pi
from matplotlib import rc

rc('text', usetex=True)  # LaTeX labels
# PATH = '/media/zade/Seagate Expansion Drive/honours_project_2020/'
PATH = '/media/zade/STRONTIUM/honours_project_2020/'
DICT_PATH = PATH + 'pickle/'
REGIMES = ['collisionless', 'braginskii with heat fluxes', 'braginskii']
DEFAULT_PROB = 'shear_alfven'


def load_data(fname, n, prob=DEFAULT_PROB):
    '''Loads data from .athdf files output from Athena++, using modules
    from the athena_read code.
    '''

    def f(n):
        return folder + '.out' + output_id + '.%05d' % n + '.athdf'

    # Input
    folder = PATH + fname + '/' + prob  # Name of output
    output_id = '2'  # Output ID (set in input file)
    filename = f(n)
    return athdf(filename)


def load_hst(fname, prob=DEFAULT_PROB):
    '''Loads data from .hst files output from Athena++, using modules
    from the athena_read code.
    '''
    hstLoc = PATH + fname + '/' + prob + '.hst'
    return hst(hstLoc)


def load_dict(fname=''):
    file = 'S_dict.pkl'
    if fname != '':
        file = fname + '_' + file

    pkl_file = open(DICT_PATH+file, 'rb')
    dict = pickle.load(pkl_file)
    pkl_file.close()
    return dict


def save_dict(dict, fname=''):
    file = 'S_dict.pkl'
    if fname != '':
        file = fname + '_' + file

    output = open(DICT_PATH+file, 'wb')
    pickle.dump(dict, output)
    output.close()


def check_dict(fname):
    return os.path.isfile(DICT_PATH+fname+'_S_dict.pkl')


def make_folder(fname):
    if not os.path.exists(fname):
        path = Path(fname)
        path.mkdir(parents=True, exist_ok=True)


def get_maxn(fname):
    return len(glob.glob(PATH+fname+'/*.athdf'))


# --- VECTOR FUNCTIONS --- #


def pnts(p):
    points = []
    for zz in range(p[0]):
        for yy in range(p[1]):
            for xx in range(p[2]):
                points.append((zz, yy, xx))
    return np.array(points)


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
    data = load_data(fname, 0, prob)
    X1 = data['RootGridX1'][1] - data['RootGridX1'][0]
    X2 = data['RootGridX2'][1] - data['RootGridX2'][0]
    X3 = data['RootGridX3'][1] - data['RootGridX3'][0]
    return abs(X1*X2*X3)  # just a check to make volume positive


# --- RMS FUNCTIONS  --- #


# FIX: Find a way to get perpendicular vector components
# Assumes that B_0 is always in the x-direction
def get_rms(fname, n, do_mhd):
    data = load_data(fname, n)
    t = data['Time']

    # Get perpendicular direction from initial conditions?
    vy, vz = data['vel2'], data['vel3']
    vel_perp_2 = (vy-np.mean(vy))**2 + (vz-np.mean(vz))**2
    if do_mhd:
        By, Bz = data['Bcc2'], data['Bcc3']
        B_perp_2 = (By-np.mean(By))**2 + (Bz-np.mean(Bz))**2
    else:
        B_perp_2 = [0]

    return t, np.sqrt(np.mean(vel_perp_2)), np.sqrt(np.mean(B_perp_2))


def plot_rms(fname, do_mhd=1):
    n_max = get_maxn(fname)

    T, V, B = [], [], []
    for i in range(n_max):
        t, v_rms, b_rms = get_rms(fname, i, do_mhd)
        T.append(t)
        V.append(v_rms)
        B.append(b_rms)

    plt.plot(T, V, T, B)
    plt.plot(T, np.mean(V)*np.ones(len(T)), ':',
             T, np.mean(B)*np.ones(len(T)), ':')
    plt.title('Time Evolution of Perpendicular Fluctuations (RMS)')
    plt.xlabel('Time (s)')
    plt.ylabel(r'$\delta_{\perp \textrm{, rms}}$')
    plt.legend([r'$\delta u_{\perp}$', r'$\delta B_{\perp}$',
                r'$\delta u_{\perp}$ mean', r'$\delta B_{\perp}$ mean'])
    plt.savefig(PATH + fname + '/' + fname.split()[-1] + '_fluc.pdf')
    plt.clf()


# --- DIMENSIONLESS PARAMETER CALCULATIONS --- #


def It_Brag(fname, nu_c, prob=DEFAULT_PROB):
    '''Code to calculate the It_Brag parameter as described in Jono's
    magneto-immutable turbulence paper. ρ and v_A are assumed to be 1.
    Relavant regimes:
        It_Brag < 1 ⟹ Δp ~ B^2
        It_Brag > 1 ⟹ Δp ≪ B^2
    '''
    data = load_data(fname, 0, prob)
    # length of box ∥ to B-field, assumed to be in x-direction
    l_prl = abs(data['RootGridX1'][1] - data['RootGridX1'][0])
    # length of box ⟂ to B-field, y- and z- lengths assumed to be the same.
    l_prp = abs(data['RootGridX2'][1] - data['RootGridX2'][0])

    # pprp and pprl are set to p_0 initially over the box
    p_0 = data['pprp'][0, 0, 0]

    # (δB⟂/B0) = (L⟂/L∥)
    It = (l_prl*nu_c / p_0) * (l_prp / l_prl)**(-2)
    return It


def beta(p0, B0):
    return 2*p0/(B0**2)


def db_int(nu_c, omega_A, beta):
    # Assumes v_A = 1, which is true when B0 and ρ=1.
    regime = find_regime(nu_c, omega_A, beta)
    if regime == 0:  # collisionless
        return 1/np.sqrt(beta)
    else:
        return np.sqrt(nu_c/(omega_A*beta))


def find_regime(nu_c, omega_A, beta):
    eta =  nu_c / omega_A
    if eta < 1:
        return 0  # collisionless
    elif eta < np.sqrt(beta):
        return 1  # brag with heat fluxes
    else:
        return 2  # brag


def calc_omega_A(Lx):
    # Assumes background magnetic field along the x axis and v_A = 1 (B0=1,ρ=1)
    # Then k_∥ = 2π/Lx and ω_A = k_∥v_A = 2π/Lx
    return 2*np.pi / Lx


def avg_B(fname, n, background, prob=DEFAULT_PROB):
    data = load_data(fname, n, prob)
    return avg_B_data(data, background)


def avg_B_data(data, background):
    grid = data['RootGridSize'][::-1]
    points = pnts(grid)

    if background:
        # magnitude squared of the background magnetic field
        # assuming perturbation is in the z direction
        # this is true for Alfven waves with k and B0 in the xy plane
        # as perturbation must be in the k×B0 = z direction
        B_data = (data['Bcc1'], data['Bcc2'], np.zeros_like(data['Bcc3']))
    else:
        # average value of B2 over the whole region
        B_data = (data['Bcc1'], data['Bcc2'], data['Bcc3'])
    B_vec = np.array([get_vec(B_data, p) for p in points])
    B_avg = np.mean(get_mag(B_vec))  # spatial average of B
    return B_avg

# --- bb:∇u PDF FUNCTIONS --- #


def prlshear_pdf(fname, n, plot_title='', do_plot=1):
    # pprp and pprl are set to p_0 initially over the box
    p_0 = load_data(fname, 0)['pprp'][0, 0, 0]
    # Convention: putting ν_c (nuc) value at end of file name
    x = fname.rfind('nuc') + 3
    ν_c = float(fname[x:])

    # get delta p, nu_c and p_0 to calculate bb:grad u
    data = load_data(fname, n)
    # Perpendicular and parallel pressures
    pprp = data['pprp']
    pprl = data['pprl']
    Δp = (pprp - pprl).flatten()

    μ_Brag = p_0/ν_c
    prlshear = Δp/μ_Brag  # bb:grad u

    if do_plot:
        n, bins, patches = plt.hist(prlshear, 100, density=True)
        plt.clf()
        # plot histogram
        plt.plot(0.5*(bins[:-1] + bins[1:]), n)
        plt.yscale('log', nonposy='clip')
        plt.title(r'PDF of $\mathbf{\hat{b}\hat{b}}:\nabla\mathbf{u}$ with '
                  + plot_title)
        plt.xlabel(r'$4\pi \mathbf{\hat{b}\hat{b}}:\nabla\mathbf{u}/\langle B^2\rangle$')
        plt.ylabel(r'$\mathcal{P}$')
        plt.savefig(PATH + fname + '/' + fname.split()[-1] + '_dp_pdfnum.pdf')
        plt.clf()
    else:
        return prlshear


def prlshearft_pdf(fname, n, plot_title='', do_plot=1):

    def dv(i, j):
        '''Calculates gradient of the components of vel_data using FFT.'''
        return np.real(fft.ifftn(K[i]*fft.fftn(vel_data[j])))

    # ft stuff
    data = load_data(fname, n)
    points = pnts(data['RootGridSize'][::-1])

    K = ft_grid(data, 0)[::-1]  # gets KX, KY, KZ
    vel_data = (data['vel1'], data['vel2'], data['vel3'])
    B_data = (data['Bcc1'], data['Bcc2'], data['Bcc3'])
    B_vec = np.array([get_vec(B_data, p) for p in points])
    b = get_unit(B_vec)
    B2avg = np.mean(get_mag(B_vec)**2)  # spatial average of B^2

    # Calculate bb:∇u = b_i b_j ∂_i u_j
    # with ∂_i u_j = ifft(K_i fft(u_j))
    # over all space
    prlshear = 0
    for i in range(3):
        for j in range(3):
            prlshear += b[:, i]*b[:, j]*dv(i, j).flatten()

    if do_plot:
        n, bins, patches = plt.hist(prlshear, 100, density=True)
        plt.yscale('log', nonposy='clip')
        plt.title(r'PDF (using Fourier method) of '
                  + r'$\mathbf{\hat{b}\hat{b}}:\nabla\mathbf{u}$ with '
                  + plot_title)
        plt.xlabel(r'$4\pi \mathbf{\hat{b}\hat{b}}:\nabla\mathbf{u}/\langle B^2\rangle$')
        plt.ylabel(r'$\mathcal{P}$')
        plt.savefig(PATH + fname + '/' + fname.split()[-1] + '_dp_pdf.pdf')
        plt.clf()
    else:
        return prlshear, B2avg


# --- FOURIER FUNCTIONS --- #


def ft_array(N):
    '''For given N, returns an array conforming to FT standard:
       [0 1 2 3 ... -N/2 -N/2+1 ... -1]
    '''
    return np.concatenate((np.arange(1, N//2, 1), [-N//2],
                           np.arange(-N//2+1, 0, 1)))


def ft_grid(data, k_grid):
    '''Creates a grid in k-space corresponding to the real grid given in data.
    k_grid is a boolean that when True calculates a regularly spaced array
    in k-space.
    '''
    # Z, Y, X
    p = (data['x3f'], data['x2f'], data['x1f'])
    Ls = [np.max(p[0]), np.max(p[1]), np.max(p[2])]
    Ns = [len(p[0]), len(p[1]), len(p[2])]

    K = {}
    for k in range(3):
        K[k] = 2j*pi/Ls[k]*ft_array(Ns[k])

    # Outputs Z, Y, X
    to_ret = np.meshgrid(K[0], K[1], K[2], indexing='ij')
    if k_grid:
        to_ret = (to_ret, np.arange(0, np.max(np.imag(K[1])), 2*pi/Ls[1]))

    return to_ret
