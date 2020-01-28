'''Code to calculate diagnostics such as root mean square and probability
   density functions for fluid simulations.
'''
import glob
import numpy as np
import numpy.fft as fft
import matplotlib.pyplot as plt
from athena_read import athdf
from math import pi
from matplotlib import rc
rc('text', usetex=True)  # LaTeX labels


def load_data(fname, n):

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
    path = '/media/zade/Seagate Expansion Drive/Summer_Project_2019/'
    filename = path + fname
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


# --- bb:∇u PDF FUNCTIONS --- #


def prlshear_pdf(fname, n, do_ft):
    # pprp and pprl are set to p_0 initially
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

    if do_ft:
        # use ft for gradient to find bb:nabla u
        # to compare with the above way of calculating
        ft_n, ft_bins, ft_patches = prlshearft_pdf(data)

    n, bins, patches = plt.hist(prlshear, 100, density=True)
    # plot histogram


def prlshearft_pdf(data):

    def dv(i, j):
        return np.real(fft.ifftn(K[i]*fft.fftn(vel_data[j])))

    def pnts(p):
        points = []
        for zz in range(p[0]):
            for yy in range(p[1]):
                for xx in range(p[2]):
                    points.append((zz, yy, xx))
        return np.array(points)

    # ft stuff
    points = pnts(data['RootGridSize'][::-1])
    K = ft_grid(data, 0)
    vel_data = (data['vel1'], data['vel2'], data['vel3'])
    B_data = (data['Bcc1'], data['Bcc2'], data['Bcc3'])
    B_vec = np.array([get_vec(B_data, p) for p in points])
    b = get_unit(B_vec)
    B2avg = np.mean(get_mag(B_vec)**2)  # spatial average of B^2

    # Calculate b_i b_j ∂_i u_j
    # with ∂_i u_j = ifft(K_i fft(u_j))
    # over all space
    prlshear = 0
    for i in range(3):
        for j in range(3):
            prlshear += b[:, i]*b[:, j]*dv(i, j).flatten()
    prlshear *= 4*pi/B2avg
    n, bins, patches = plt.hist(prlshear, 100, density=True)
    plt.yscale('log', nonposy='clip')


# --- FOURIER FUNCTIONS --- #


def ft_array(N):
    '''For given N, returns an array conforming to FT standard:
       [0 1 2 3 ... -N/2 -N/2+1 ... -1]
    '''
    return np.concatenate((np.arange(1, N//2, 1), [-N//2],
                           np.arange(-N//2+1, 0, 1)))


def ft_grid(data, k_grid):
    p = (data['x3f'], data['x2f'], data['x1f'])
    Ls = [np.max(p[0]), np.max(p[1]), np.max(p[2])]
    Ns = [len(p[0]), len(p[1]), len(p[2])]

    K = {}
    for k in range(3):
        K[k] = 2j*pi/Ls[k]*ft_array(Ns[k])

    to_ret = np.meshgrid(K[0], K[1], K[2], indexing='ij')
    if k_grid:
        to_ret = (to_ret, np.arange(0, np.max(np.imag(K[1])), 2*pi/Ls[1]))

    return to_ret
