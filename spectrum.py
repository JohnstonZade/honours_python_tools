# -- IMPORTS -- #

# -- for testing on laptop -- #
import sys
sys.path.insert(1, '/home/zade/athena/vis/python')
# -- for testing on laptop -- #
from athena_read import athdf, vtk
import numpy as np
from math import pi
# plotting
# fft

# TODO: Add unit tests


def spectrum(folder_name, plot_title, do_mhd=1, do_full_calc=1):
    # Add file names

    folder = '/data/johza721/output/MHDTurb/' + folder_name  # Folder of output
    folder += 'Turb'  # Name of output
    output_id = '2'  # Output ID (set in input file)
    is_cont = 'cont' in folder_name

    # Getting turnover time and converting to file number
    τ = get_turnover_time(folder_name, is_cont)
    τ_file = 10*int(τ)

    # File numbers to average over
    # Figure out more adaptive way for decay
    nums = range(τ_file, 6*τ_file+1) if is_cont else range(τ_file, 301)
    return nums

    if do_full_calc:
        # read hdf5 files
        fname = lambda n: folder + '.out' + output_id + '.%05d' % n + '.athdf'
        read_f = lambda f: athdf(f)  # put athena read file here

        # create grid of K from first time step
        # filename = '/media/zade/Seagate Expansion Drive/Summer_Project_2019/'
        # filename += 'hydro_cont_turb_32/Turb.out2.00128.athdf'

        D = read_f(fname(nums[0]))
        x = D['x1f']
        y = D['x2f']
        z = D['x3f']
        dx = x[1] - x[0]
        dy = y[1] - y[0]
        dz = z[1] - z[0]
        Ls = [2, 3, 4] # [np.max(x), np.max(y), np.max(z)]
        Ns = [32, 32, 32] # [len(x), len(y), len(z)]

        # create grid in k-space
        # not sure if best way to get FT standard [0 1 2 3...-N/2 -N/2+1...-1]
        K = {}
        for k in range(3):
            K[k] = 2j*pi/Ls[k]*ft_array(Ns[k])

        KX, KY, KZ = np.meshgrid(K[0], K[1], K[2])
        Kprl = np.abs(KX)
        Kperp = np.sqrt(np.abs(KY)**2 + np.abs(KZ)**2)
        Kmag = np.sqrt(Kprl**2+Kperp**2)
        Kspec = Kmag

        kgrid = np.arange(0, np.max(np.imag(K[1])), 2*pi/Ls[1])
        # TODO: make a structure/dict to hold spectrum?

        # normalize modes

        # average over snapshots in nums

        # save spectrum data somehow
    else:
        # load spectrum data
        plpl = 2
    # plot spectrum
    # save figure
    # plot energy time evolution
    # save figure


def ft_array(N):
    '''For given N, returns an array conforming to FT standard:
       [0 1 2 3 ... -N/2 -N/2+1 ... -1]
    '''
    return np.concatenate((np.arange(1, N//2, 1), [-N//2], np.arange(-N//2+1, 0, 1)))


def ath_plot_hst(folder_name, plot_title):
    # Add file names
    # Import data from .hst files, use athena_read.py
    # Get turnover time
    return None


def get_turnover_time(folder_name, is_cont):
    # Add file names
    # Import data from .hst files

    # if is_cont
    # average turnover time in saturated state

    # else
    # calculate initial turnover time

    # return turnover time
    return 12.56
