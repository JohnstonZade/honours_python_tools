import glob
import numpy as np
import matplotlib.pyplot as plt
from math import pi
from matplotlib import rc
rc('text', usetex=True)  # LaTeX labels

# Change to athena_read dir
import sys
sys.path.insert(1, '/home/zade/athena/vis/python')
# Thunderbird
# sys.path.insert(1, '/nfs_physics/users/stud/johza721/athena/vis/python')
from athena_read import athdf  # no h5py on thunderbird



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


def get_unit(v):
    '''Calculates unit vector.'''
    v_mag = get_mag(v)
    return np.array([v[i]/v_mag[i] for i in range(len(v))])


def get_mag(X):
    '''Returns the magnitude of the given vector.'''
    x = np.array(X)
    return np.sqrt(np.sum(x**2, axis=1))


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


# --- PDF --- #


def shearstrain_pdf(fname, n, do_ft):
    # get delta p, nu_c and p_0 to calculate bb:nabla u

    if do_ft:
        # use ft for gradient to find bb:nabla u
        ft_pdf()

    # plot/return histogram
    return 1


def ft_pdf():
    # ft stuff
    # plot/return histogram
    return 1


# --- FOURIER --- #


def ft_array(N):
    '''For given N, returns an array conforming to FT standard:
       [0 1 2 3 ... -N/2 -N/2+1 ... -1]
    '''
    return np.concatenate((np.arange(1, N//2, 1), [-N//2],
                           np.arange(-N//2+1, 0, 1)))


# Very tentative, don't think this is right for a Fourier transform
# Will work for spectrum but am wanting to try get a velocity Fourier transform
# to get Kx from ux etc
def ft(v1, v2, v3, k_grid):
    Ls = [np.max(v1), np.max(v2), np.max(v3)]
    Ns = [len(v1), len(v2), len(v3)]

    K = {}
    for k in range(3):
        K[k] = 2j*pi/Ls[k]*ft_array(Ns[k])

    to_ret = np.meshgrid(K[0], K[1], K[2])
    if k_grid:
        to_ret = (to_ret, np.arange(0, np.max(np.imag(K[1])), 2*pi/Ls[1]))

    return to_ret
