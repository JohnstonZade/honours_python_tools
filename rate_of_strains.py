import diagnostics as diag
import numpy as np
import numpy.fft as fft
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)  # LaTeX labels


def identity(i, j):
    return 1.0 if i == j else 0.0


def grad_u(K, u, i, j):
    return np.real(fft.ifftn(K[i]*fft.fftn(u[j])))


def bbgu_calc(b, du, i, j):
    return b[:, i]*b[:, j] * du[i, j, :]


def bbgu(b, du, no_points):
    sums = range(3)
    bbgu_elements = np.array([[bbgu_calc(b, du, i, j) for i in sums]
                              for j in sums])
    return np.array([np.sum(bbgu_elements[:, :, p]) for p in range(no_points)])


def Gprp_uprp_calc(b, du, i, j):
    return (identity(i, j) - b[:, i]*b[:, j]) * du[i, j, :]


def Gprp_uprp(b, du, no_points):
    sums = range(3)
    Gprp_uprp_elements = np.array([[Gprp_uprp_calc(b, du, i, j) for i in sums]
                                   for j in sums])
    return np.array([np.sum(Gprp_uprp_elements[:, :, p])
                     for p in range(no_points)])


def plot_strains(output_dir, save_path, fname, prob):
    bin_size = 400
    n_max = diag.get_maxn(output_dir)  # 4 alfven periods
    n_max = int(n_max / 2)  # 2 alfven periods
    n_bins = [np.zeros(bin_size),
              np.zeros(bin_size), np.zeros(bin_size), np.zeros(bin_size)]
    bins = [np.zeros(bin_size+1), np.zeros(bin_size+1),
            np.zeros(bin_size+1), np.zeros(bin_size+1)]
    labels = [r'$|\nabla\mathbf{u}|/u_\textrm{rms}$',
              r'$|\nabla_{\|}\mathbf{u}_{\|}|/u_\textrm{rms}$',
              r'$|\nabla_{\perp}\mathbf{u}_{\|}|/u_\textrm{rms}$',
              r'$|\nabla_{\|}\mathbf{u}_{\perp}|/u_\textrm{rms}$']

    for n in range(n_max):
        data = diag.load_data(output_dir, n, prob=prob)
        points = diag.pnts(data['RootGridSize'][::-1])
        no_points = points.shape[0]
        K = diag.ft_grid(data, 0)[::-1]  # Kx Ky Kz

        u = (data['vel1'], data['vel2'], data['vel3'])
        u_mag = diag.get_mag(u)
        u_rms = np.sqrt(np.mean(u_mag**2))

        sums = range(3)
        du = np.array([[grad_u(K, u, i, j).flatten() for j in sums]
                       for i in sums])
        du_mag_2 = du[0, 0]**2 + du[1, 0]**2 + du[2, 0]**2\
                 + du[0, 1]**2 + du[1, 1]**2 + du[2, 1]**2\
                 + du[0, 2]**2 + du[1, 2]**2 + du[2, 2]**2

        B_data = (data['Bcc1'], data['Bcc2'], data['Bcc3'])
        B_vec = np.array([diag.get_vec(B_data, p) for p in points])
        b = diag.get_unit(B_vec)

        dprluprl = bbgu(b, du, no_points)**2  # / du_mag_2
        # dprpuprp = Gprp_uprp(b, du, no_points)**2   # / du_mag_2

        dprpuprl = (du[0, 0, :]*b[:, 0] + du[0, 1, :]*b[:, 1] + du[0, 2, :]*b[:, 2])**2\
                 + (du[1, 0, :]*b[:, 0] + du[1, 1, :]*b[:, 1] + du[1, 2, :]*b[:, 2])**2\
                 + (du[2, 0, :]*b[:, 0] + du[2, 1, :]*b[:, 1] + du[2, 2, :]*b[:, 2])**2\
                 - dprluprl

        dprluprp = (b[:, 0]*du[0, 0, :] + b[:, 1]*du[1, 0, :] + b[:, 2]*du[2, 0, :])**2\
                 + (b[:, 0]*du[0, 1, :] + b[:, 1]*du[1, 1, :] + b[:, 2]*du[2, 1, :])**2\
                 + (b[:, 0]*du[0, 2, :] + b[:, 1]*du[1, 2, :] + b[:, 2]*du[2, 2, :])**2\
                 - dprluprl

        gradients = [du_mag_2, dprluprl, dprpuprl, dprluprp]

        for i, grad in enumerate(gradients):
            grad /= u_rms**2
            n_bins_temp, bins_temp = np.histogram(np.sqrt(grad),
                                                  bin_size, density=True)
            n_bins[i] += n_bins_temp
            bins[i] += bins_temp
        print('Done ' + str(n))

    for i, n_bin in enumerate(n_bins):
        n_bin /= n_max
        bins[i] /= n_max
        plt.semilogy(0.5*(bins[i][:-1] + bins[i][1:]), n_bin, label=labels[i])
    plt.legend()
    title_string = 'MHD' if fname.split('_')[-1] == 'mhd' else ''
    plt.xlabel(r'$|\nabla \mathbf{u}|$')
    plt.ylabel(r'$\mathcal{P}$')
    plt.title(title_string)
    plt.savefig(save_path + '/' + fname + '_strainrate.pdf')
    plt.savefig(save_path + '/' + fname + '_strainrate.png')
    plt.close()
