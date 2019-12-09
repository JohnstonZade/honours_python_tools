'''Code to calculate structure function of a fluid simulation.'''
# Change to athena_read dir
import sys
sys.path.insert(1, '/home/zade/athena/vis/python')
# Thunderbird
# sys.path.insert(1, '/nfs_physics/users/stud/johza721/athena/vis/python')
import os
from athena_read import athdf  # no h5py on thunderbird
import numpy as np
from numpy.random import randint, random
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)  # LaTeX labels


def get_mag(X):
    '''Returns the magnitude of the given vector.'''
    x = np.array(X)
    return np.sqrt(np.sum(x**2, axis=1))


def generate_points(grid, N):
    '''Generate N pairs of random points on the grid.'''
    N_p = int(N)
    L1 = randint(0, grid, size=(N_p, len(grid)))
    # Generate second set of random points to bias closer to points in L1
    r_pow = 6
    L2 = L1 + np.ceil(grid*random(size=(N_p, len(grid)))**r_pow).astype(int)
    # Ensure points are in the grid
    L2 = np.mod(L2, grid)
    return L1, L2


def get_vec(v, p):
    '''Returns the vector components at a given point.'''
    tp = tuple(p)
    v1 = v[0][tp]  # x-component
    v2 = v[1][tp]  # y-component
    v3 = v[2][tp]  # z-component
    return np.array([v1, v2, v3])


def select_y(x, x_bin, y, i, mask=[], use_mask=0, return_mask=0):
    '''Selects and averages the y vector values with the condition that
    the corresponding x values satisfy x_bin[i] <= x <= x_bin[i+1].

    If a mask is provided, it is used on both x and y so that the order of
    the elements is unchanged when selecting based on x_bin.
    '''

    if use_mask:
        x = x[mask]
        y = y[mask]
    assert i < len(x_bin)-1, 'Can\'t access element at len(x_bin)'
    x_mask = np.logical_and(x_bin[i] <= x, x < x_bin[i+1])

    return x_mask if return_mask else y[x_mask]


def get_mean(x, x_bin, y, i, mask=[], use_mask=0):
    '''Calculate the mean of the selection.'''
    y_sel = select_y(x, x_bin, y, i, mask, use_mask)

    # Stops error for mean of empty slice
    # Might be a better way to handle this
    if len(y_sel) == 0:
        return float('nan')

    return np.mean(y_sel)


def get_unit(v):
    '''Calculates unit vector.'''
    v_mag = get_mag(v)
    return np.array([v[i]/v_mag[i] for i in range(len(v))])


def get_l_perp(L1, L2, l, B):
    # Calculate average B field between point pairs
    B1_vec = np.array([get_vec(B, p) for p in L1])
    B2_vec = np.array([get_vec(B, p) for p in L2])
    B_mean = 0.5*(B1_vec + B2_vec)

    # Dot product of unit vectors to get cos(θ)
    cθ = abs(np.sum(get_unit(B_mean)*get_unit(l), axis=1))
    θ_data = np.arccos(cθ)
    θ = np.linspace(0, 90, 4, endpoint=True)  # 7 default
    θlen = len(θ) - 1
    θ_rad = (np.pi/180)*θ

    l_mask = [select_y(θ_data, θ_rad, get_mag(l), i, return_mask=1)
              for i in range(θlen)]
    return θ, l_mask


def calc_struct(L1, L2, v, l_mag, L_max, mask=[], use_mask=0):
    # For each pair of position vectors x1, x2
    # Get velocity/Bfield vectors v1, v2 at each point
    # Calculate Δv2 = abs(v1 - v2)**2
    # We now have a mapping of l to Δv2 <- structure function
    v1_vec = np.array([get_vec(v, p) for p in L1])
    v2_vec = np.array([get_vec(v, p) for p in L2])
    Δv_vec = v1_vec - v2_vec
    Δv_mag2 = get_mag(Δv_vec)**2

    # Bin and plot structure function
    # Plot in the middle of the grid points otherwise the size of arrays
    # won't match up.
    N_l = 30
    l_bin = np.logspace(np.log10(2*L_max/N_l), np.log10(L_max/2), N_l+1)
    l_grid = 0.5*(l_bin[:-1] + l_bin[1:])
    Δv_avg = np.array([get_mean(l_mag, l_bin, Δv_mag2, i, mask, use_mask)
                       for i in range(N_l)])

    print('Structure Function Calculated')
    return l_grid, Δv_avg


def plot_MHD(l, t, titles, vels, Bs, fname):
    # Check whether folder to save exists
    # Seagate
    filename = '/media/zade/Seagate Expansion Drive/Summer_Project_2019/'\
                + 'figs/' + fname

    if not os.path.exists(filename):
        os.makedirs(filename)

    for i in range(len(titles)):
        # plt.subplot(3, 2, i+1)
        plt.loglog(l, vels[i], l, Bs[i], l, l**(2/3))
        plt.title(r'$S_2(l)$ with ' + titles[i])
        plt.xlabel(r'log($l$)')
        plt.ylabel(r'log($S_2(l)$))')
        plt.legend(['Vel Structure Function', 'B-field Structure Function',
                    r'$l^{2/3}$'])
        plt.savefig(filename + '/t' + t + '_' + str(i) + '.png')
        plt.clf()
    print('Plotted MHD')


def plot_struct(l_grid, v_avg, t, fname):
    filename = '/media/zade/Seagate Expansion Drive/Summer_Project_2019/'\
                + 'figs/' + fname

    if not os.path.exists(filename):
        os.makedirs(filename)

    plt.loglog(l_grid, v_avg, l_grid, l_grid**(2/3))
    plt.title(r'$S_2(l)$ at $t=$ ' + t)
    plt.xlabel(r'log($l$)')
    plt.ylabel(r'log($S_2(l)$))')
    plt.legend(['Structure Function', r'$l^{2/3}$'])
    plt.savefig(filename + '/struct_t' + t + '.png')
    plt.clf()


def structure_function(fname, n, do_mhd=0, N=1e6, do_ldist=0):
    '''Calculates and plots structure function.
    Takes about 20s for a 128 cube grid file with 1 million random pairs.
    '''

    def get_length():
        '''Returns the dimensions of the simulation.'''
        names = ['RootGridX3', 'RootGridX2', 'RootGridX1']
        return np.array([data[name][1]-data[name][0] for name in names])

    def get_point(p):
        '''Returns coordinate of given grid point.'''
        tp = tuple(p)
        return np.array([zz[tp], yy[tp], xx[tp]])

    def f(n):
        return folder + '.out' + output_id + '.%05d' % n + '.athdf'

    # Read in data and set up grid

    # Seagate
    folder = '/media/zade/Seagate Expansion Drive/Summer_Project_2019/'
    if do_mhd:
        folder += 'MHD/'

    # Thunderbird
    # folder = '/data/johza721/output/MHDTurb/'
    # Input
    folder += fname + '/Turb'  # Name of output
    output_id = '2'  # Output ID (set in input file)
    filename = f(n)
    data = athdf(filename)

    # Following (z, y, x) convention from athena_read
    grid = data['RootGridSize'][::-1]
    t = '{:.1f}'.format(data['Time']) + ' s'

    zz, yy, xx = np.meshgrid(data['x3f'], data['x2f'], data['x1f'],
                             indexing='ij')
    vel_data = (data['vel1'], data['vel2'], data['vel3'])
    if do_mhd:
        B_data = (data['Bcc1'], data['Bcc2'], data['Bcc3'])

    L1, L2 = generate_points(grid, N)
    print('Generated points')

    # Get actual position vector components for each pair of grid points
    # and difference vector between them
    x1_vec = np.array([get_point(p) for p in L1])
    x2_vec = np.array([get_point(p) for p in L2])
    l_vec = x1_vec - x2_vec
    # Find distance between each pair of points
    l_mag = get_mag(l_vec)
    L = np.max(get_length())
    print('Lengths calculated')

    # Check distribution of l's
    if do_ldist:
        n, bins, patches = plt.hist(l_mag, 100)
        plt.title(r'Distribution of $l$ at $t =$ ' + t)
        plt.xlabel('Distance between points')
        plt.ylabel('Counts')
        plt.show()

    # MHD: Get magnetic field components
    # Get l and B unit vectors then get
    # cos(theta) = l dot b
    # bin with theta = [0, 15, 30, 45, 60, 75, 90]
    # 0-15 is l_parallel
    # 45-90 is l_perp
    # display structure function for each bin (eg 0-15)
    # using both average velocity at each pair of points
    # and average B field at each pair of points (12 total)
    if do_mhd:
        θ, l_mask = get_l_perp(L1, L2, l_vec, B_data)
        titles, vels, Bs = [], [], []
        l_grid = calc_struct(L1, L2, vel_data, l_mag, L)[0]

        for i, l_m in enumerate(l_mask):
            titles.append(str(θ[i]) + r'$^\circ$ $\leq \theta <$ '
                          + str(θ[i+1]) + r'$^\circ$' + ' at t = ' + t)
            vels.append(calc_struct(L1, L2, vel_data, l_mag, L, l_m, 1)[1])
            Bs.append(calc_struct(L1, L2, B_data, l_mag, L, l_m, 1)[1])
            print('{:.2f}'.format((i+1)/len(l_mask)*100) + '% done')

        print('Plotting MHD')
        plot_MHD(l_grid, t, titles, vels, Bs, fname)

    print('Calculating full velocity structure function')
    l_grid, Δv_avg = calc_struct(L1, L2, vel_data, l_mag, L)
    print('Plotting Struct')
    plot_struct(l_grid, Δv_avg, t, fname)
