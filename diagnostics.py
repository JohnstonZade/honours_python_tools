import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from structure_function import load_data
# rc('text', usetex=True)  # LaTeX labels


# FIX:
# Assumes that B_0 is always in the x-direction
def get_rms(fname, n, do_mhd):
    data = load_data(fname, n)
    t = data['Time']
    # Get perpendicullar direction from initial conditions?
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
    # plt.plot(T, 0.5*np.ones(len(T)), ':')
    plt.plot(T, np.mean(B)*np.ones(len(T)), ':', T, np.mean(V)*np.ones(len(T)), ':')
    plt.title('Time Evolution of Perpendicular Fluctuations (RMS)')
    plt.xlabel('Time (s)')
    plt.ylabel(r'$\delta_{\perp \textrm{, rms}}$')
    plt.legend([r'$\delta u_{\perp}$', r'$\delta B_{\perp}$',
                r'$\delta u_{\perp}$ mean', r'$\delta B_{\perp}$ mean'])
    plt.show()


fname = 'MHD/cgl_cont_turb_12864'
plot_rms(fname)
