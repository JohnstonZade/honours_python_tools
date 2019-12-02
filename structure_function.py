'''Code to calculate structure function of a fluid simulation.
Not sure if it works on Thunderbird as it seems to be running
Python 2.
'''
# Change to athena_read dir
import sys
sys.path.insert(1, '/home/zade/athena/vis/python')
from athena_read import athdf
from numpy import array, ceil, linspace, log10, logspace, logical_and, max
from numpy import mean, meshgrid, mod, sqrt, sum, where, zeros
from numpy.random import randint, random
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)  # LaTeX labels


def get_mag(X):
    '''Returns the magnitude of the given vector.'''
    x = array(X)
    return sqrt(sum(x**2, axis=1))


# Check todos
def structure_function(filename, n, N_p=1e6, do_ldist=0, do_plot=0):
    '''Calculates and plots structure function
    assuming isotropy of turbulence.

    Takes about 20s for a 128 cube grid file with 1 million random pairs.
    '''
    def get_mean(i):
        '''Selects and averages the difference in velocity values squared at
        between two points with distance l between them.
        Selection depends on whether l_bin[i] <= l < l_bin[i+1].
        '''
        l_logic = where(logical_and(l_bin[i] <= l_mag, l_mag < l_bin[i+1]))

        # Stops error for mean of empty slice
        # Might be a better way to handle this
        if len(Δu2_mag[l_logic]) == 0:
            return float('nan')

        return mean(Δu2_mag[l_logic])

    def get_length():
        '''Returns the dimensions of the simulation.'''
        names = ['RootGridX3', 'RootGridX2', 'RootGridX1']
        return array([data[name][1]-data[name][0] for name in names])

    def get_vel(p):
        '''Returns the velocity vector at a given point.'''
        tp = tuple(p)
        v1 = vel_data[0][tp]  # x-component
        v2 = vel_data[1][tp]  # y-component
        v3 = vel_data[2][tp]  # z-component
        return array([v1, v2, v3])

    def get_point(p):
        '''Returns coordinate of given grid point.'''
        tp = tuple(p)
        return array([zz[tp], yy[tp], xx[tp]])

    def fname(n):
        return folder + '.out' + output_id + '.%05d' % n + '.athdf'

    # Read in data and set up grid
    # filename = '/media/zade/Seagate Expansion Drive/Summer_Project_2019/'
    # filename += 'hydro_cont_turb_128/Turb.out2.00030.athdf'
    folder = '/media/zade/Seagate Expansion Drive/Summer_Project_2019/'
    folder += filename + '/Turb'  # Name of output
    output_id = '2'  # Output ID (set in input file)

    filename = fname(n)
    data = athdf(filename)

    # Following (z, y, x) convention from athena_read
    grid = data['RootGridSize'][::-1]
    zz, yy, xx = meshgrid(data['x3f'], data['x2f'], data['x1f'], indexing='ij')
    t = '{:.3f}'.format(data['Time']) + ' s'

    # Get N_p pairs of random points on the grid
    # Assuming turbulence is isotropic
    N_p = int(N_p)
    L1 = randint(0, grid, size=(N_p, len(grid)))
    # Generate second set of random points to bias closer to points in L1
    r_pow = 6
    L2 = L1 + ceil(grid*random(size=(N_p, len(grid)))**r_pow).astype(int)
    L2 = mod(L2, grid)
    print("Points generated")

    # Get actual position vector components for each pair of grid points
    # and difference vector between them
    x1_vec = array([get_point(p) for p in L1])
    x2_vec = array([get_point(p) for p in L2])
    l_vec = x1_vec - x2_vec
    # Find distance between each pair of points
    l_mag = get_mag(l_vec)
    print('Distances calculated')

    # Check distribution of l's
    if do_ldist:
        n, bins, patches = plt.hist(l_mag, 100)
        plt.title(r'Distribution of $l$ at $t =$ ' + t)
        plt.xlabel('Distance between points')
        plt.ylabel('Counts')
        plt.show()

    # For each pair of position vectors x1, x2
    # Get velocity vectors u1, u2 at each point
    # Calculate Δu2 = abs(u1 - u2)**2
    # We now have a mapping of l to Δu2! <- structure function!
    vel_data = (data['vel1'], data['vel2'], data['vel3'])
    u1_vec = array([get_vel(p) for p in L1])
    u2_vec = array([get_vel(p) for p in L2])
    Δu2_vec = u1_vec - u2_vec
    Δu2_mag = get_mag(Δu2_vec)**2
    print('Velocities calculated')

    # Bin and plot structure function
    # Plot in the middle of the grid points otherwise the size of arrays
    # won't match up.
    N_l = 30
    L = max(get_length())
    l_bin = logspace(log10(2*L/N_l), log10(L/2), N_l+1)
    l_grid = 0.5*(l_bin[:-1] + l_bin[1:])
    Δu_avg = array([get_mean(i) for i in range(N_l)])

    # Update labels
    if do_plot:
        plt.loglog(l_grid, Δu_avg, l_grid, l_grid**(2/3))
        plt.title(r'$S_2(l)$ at $t=$ ' + t)
        plt.xlabel(r'log($l$)')
        plt.ylabel('log Structure Function')
        plt.legend(['Structure Function', r'$l^{2/3}$'])
        plt.show()
    else:
        return l_grid, Δu_avg, t

# structure_function('hydro_cont_turb_128', 1, do_plot=1)
