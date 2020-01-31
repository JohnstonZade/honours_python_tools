'''Code to calculate diagnostics such as root mean square and probability
   density functions for fluid simulations.
'''
import glob
import pickle
import numpy as np
import numpy.fft as fft
import matplotlib.pyplot as plt
from athena_read import athdf, hst
from math import pi
from matplotlib import rc
rc('text', usetex=True)  # LaTeX labels


def load_data(fname, n):
    '''Loads data from .athdf files output from Athena++, using modules
    from the athena_read code.
    '''

    def f(n):
        return folder + '.out' + output_id + '.%05d' % n + '.athdf'

    # Seagate
    folder = '/media/zade/Seagate Expansion Drive/Summer_Project_2019/'

    # Thunderbird
    # folder = '/data/johza721/output/MHDTurb/'

    # Input
    folder += fname + '/Turb'  # Name of output
    output_id = '2'  # Output ID (set in input file)
    filename = f(n)
    return athdf(filename)


def load_hst(fname):
    '''Loads data from .hst files output from Athena++, using modules
    from the athena_read code.
    '''
    # Seagate
    folder = '/media/zade/Seagate Expansion Drive/Summer_Project_2019/'

    hstLoc = folder + fname + '/Turb.hst'
    return hst(hstLoc)


def load_dict(fname=''):
    file = 'S_dict.pkl'
    if fname != '':
        file = fname + '_' + file

    pkl_file = open(file, 'rb')
    dict = pickle.load(pkl_file)
    pkl_file.close()
    return dict


def save_dict(dict, fname=''):
    file = 'S_dict.pkl'
    if fname != '':
        file = fname + '_' + file

    output = open(file, 'wb')
    pickle.dump(dict, output)
    output.close()


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
    folder = '/media/zade/Seagate Expansion Drive/Summer_Project_2019/'
    filename = folder + fname
    n_max = len(glob.glob(filename+'/*.athdf'))

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
    plt.savefig('animate/' + fname + '/' + fname + '_fluc.png')
    plt.clf()


# TODO: implement saving plot
def plot_energy_evo(fname, plot_title='test'):

    def get_vol(fname):
        data = load_data(fname, 0)
        X1 = data['RootGridX1'][1] - data['RootGridX1'][0]
        X2 = data['RootGridX2'][1] - data['RootGridX2'][0]
        X3 = data['RootGridX3'][1] - data['RootGridX3'][0]
        return abs(X1*X2*X3)  # just a check to make volume positive

    hstData = load_hst(fname)
    vol = get_vol(fname)
    t = hstData['time']
    # tau = get_turnover_time(fname, 0)[0]
    # t /= tau

    KE = (hstData['2-KE'] + hstData['3-KE']) / vol  # kinetic energy
    ME = (hstData['2-ME'] + hstData['3-ME']) / vol  # magnetic energy
    TE = hstData['tot-E'] / vol                     # thermal energy
    norm = 1
    plt.semilogy(t, KE/norm, t, ME/norm, t, TE/norm)


# TODO: implement calculation
def It_Brag(fname):
    It = 1
    return It


# --- bb:∇u PDF FUNCTIONS --- #


# TODO: focus on getting working
def prlshear_pdf(fname, n):
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

    n, bins, patches = plt.hist(prlshear, 100, density=True)
    # plot histogram


# TODO: make independent function
def prlshearft_pdf(fname, n):

    def dv(i, j):
        '''Calculates gradient of the components of vel_data using FFT.'''
        return np.real(fft.ifftn(K[i]*fft.fftn(vel_data[j])))

    def pnts(p):
        points = []
        for zz in range(p[0]):
            for yy in range(p[1]):
                for xx in range(p[2]):
                    points.append((zz, yy, xx))
        return np.array(points)

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
    prlshear *= 4*pi/B2avg

    n, bins, patches = plt.hist(prlshear, 100, density=True)
    plt.yscale('log', nonposy='clip')
    plt.xlabel(r'$4\pi \mathbf{\hat{b}\hat{b}}:\nabla\mathbf{u}/\langle B^2\rangle$')
    plt.ylabel(r'$\mathcal{P}$')
    # plt.title('PDF for bb:Gu with [CALCULATE AND ADD PARAMETERS HERE]')

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
