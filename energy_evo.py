import matplotlib.pyplot as plt
import diagnostics as diag
from matplotlib import rc
rc('text', usetex=True)  # LaTeX labels
default_prob = diag.DEFAULT_PROB


def plot_energy_evo(output_dir, fname, nu_c, p0,
                    theta=None,
                    B0=1.0,
                    Lx=1.0,
                    plot_tot_only=1,
                    prob=default_prob):
    t, KE, ME = get_energy_data(output_dir, fname, prob, B0=B0, Lx=Lx)
    tot_E = KE + ME  # ignoring constant background magnetic field

    omega_A = diag.calc_omega_A(Lx)
    beta = diag.beta(p0, B0)
    regime = diag.REGIMES[diag.find_regime(nu_c, omega_A, beta)]
    title_str = 'Evolution of energy in the ' + regime + ' regime. '
    if theta is not None:
        title_str += r'$\theta =$ ' + str(round(theta, 2))
    norm = 1

    if plot_tot_only:
        plt.semilogy(t, tot_E/norm)
    else:
        plt.semilogy(t, KE/norm, t, ME/norm, t, tot_E/norm)
        plt.legend([r'$E_K$: Kinetic Energy',
                    r'$E_{B_\perp}$: Perpendicular Magnetic Energy',
                    r'$E_\textrm{tot} = E_K + E_{B_\perp}$'])
    plt.title(title_str)
    plt.xlabel(r'$t/\tau_A$')
    plt.ylabel('Energy Density')
    plt.savefig(diag.PATH + output_dir + '/' + fname + '_energy_evo.pdf')
    plt.grid()
    plt.close()


def get_energy_data(output_dir,
                    prob=default_prob,
                    volume=None,
                    B0=1.0,
                    Lx=1.0):
    hstData = diag.load_hst(output_dir, prob)
    vol = diag.get_vol(output_dir, prob) if volume is None else volume
    t_A = Lx  # true when v_A = 1 and B0 along the x axis
    t = hstData['time'] / t_A  # scale by Alfven period
    KE = (hstData['1-KE'] + hstData['2-KE'] + hstData['3-KE']) / vol  # kinetic
    ME = hstData['3-ME'] / vol  # perp magnetic energy
    return t, KE, ME
