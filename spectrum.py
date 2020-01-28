'''Code to calculate the energy spectrum of turbulence
   in a fluid simulation.
'''
import numpy as np
from diagnostics import ft_grid, load_data
from math import pi
# plotting
# fft

# TODO: Add unit tests


def spectrum(fname, plot_title, do_mhd=1, do_full_calc=1):
    is_cont = 'cont' in fname

    # Getting turnover time and converting to file number
    tau = get_turnover_time(fname, is_cont)
    tau_file = 10*int(tau)

    # File numbers to average over
    # Figure out more adaptive way for decay
    nums = range(tau_file, 6*tau_file+1) if is_cont else range(tau_file, 301)
    return nums

    if do_full_calc:
        # create grid of K from first time step
        # filename = '/media/zade/Seagate Expansion Drive/Summer_Project_2019/'
        # filename += 'hydro_cont_turb_32/Turb.out2.00128.athdf'

        data = load_data(fname, tau_file)
        KX, KY, KZ, kgrid = ft_grid(data, 1)
        Kprl = np.abs(KX)
        Kperp = np.sqrt(np.abs(KY)**2 + np.abs(KZ)**2)
        Kmag = np.sqrt(Kprl**2+Kperp**2)
        Kspec = Kmag

        # normalize modes

        # average over snapshots in nums

        # save spectrum data somehow
    else:
        # load spectrum data
        return 2
    # plot spectrum
    # save figure
    # plot energy time evolution
    # save figure


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
