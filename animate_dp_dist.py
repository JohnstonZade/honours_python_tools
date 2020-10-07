import diagnostics as diag
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import rc
rc('text', usetex=True)  # LaTeX labels


def plot_dpdist_all_time(output_dir, a_save, fname,
                         B0=1.0,
                         Lx=1.0,
                         prob='shear_alfven'):
    t_dp = diag.get_time_vec(output_dir, prob=prob)  # time vector
    colormap = plt.cm.viridis
    norm = mpl.colors.Normalize(vmin=t_dp[0], vmax=t_dp[-1])

    for n in [5, 20, 40, 80]:  # t/τ_A = 0.25, 1.0, 2.0, 4.0
        t = '{:.2f}'.format(t_dp[n])  # scale by Alfven time
        n_hist, bins = diag.dp_pdf(output_dir, n, prob=prob)
        plt.semilogy(0.5*(bins[:-1] + bins[1:]), n_hist, color=colormap(norm(t_dp[n])),
                     label=r'$t/\tau_A= \ $' + t)
    plt.xlabel(r'$4\pi\Delta p / \langle B^2 \rangle$')
    plt.ylabel(r'$\mathcal{P}$')
    plt.xlim(-1.1, 0.6)
    plt.ylim(1e-5, 1e3)
    plt.legend()
    plt.grid()
    plt.axvline(-1, ls='--', c='r')  # firehose
    plt.axvline(0.5, ls='--', c='k')  # mirror

    plt.savefig(a_save + fname + '_dp_pdf_evo.pdf')
    plt.savefig(a_save + fname + '_dp_pdf_evo.png')
    plt.close()


def plot_dpdist_and_bbgudist(output_dir, a_save, fname, regime, nu_c, p0,
                             prob='shear_alfven', show_plot=0):
    t_dp = diag.get_time_vec(output_dir, prob=prob)  # time vector
    colours = ['b', 'k']

    for i, n in enumerate([20, 80]):  # t/τ_A = 1.0, 4.0
        t = '{:.2f}'.format(t_dp[n])  # scale by Alfven time
        dp_hist, bins = diag.dp_pdf(output_dir, n, prob=prob)
        bbgu_hist = diag.bbgu_pdf(output_dir, n, prob=prob,
                                           dp_scale=1,
                                           regime=regime,
                                           nu_c=nu_c,
                                           p0=p0)[0]

        x = 0.5*(bins[:-1] + bins[1:])

        plt.semilogy(x, dp_hist, color=colours[i],
                     label=r'$t/\tau_A= \ $' + t)
        plt.semilogy(x, bbgu_hist, color=colours[i], ls=':')

    plt.xlabel(r'$4\pi\Delta p / \langle B^2 \rangle$')
    plt.ylabel(r'$\mathcal{P}$')
    plt.xlim(-1.1, 0.6)
    plt.ylim(1e-5, 1e3)
    # plt.legend()
    plt.grid()
    plt.axvline(-1, ls='--', c='r')  # firehose
    plt.axvline(0.5, ls='--', c='k')  # mirror

    if show_plot:
        plt.show()
    else:
        plt.savefig(a_save + fname + '_dpbbgu_pdf_evo.pdf')
        plt.savefig(a_save + fname + '_dpbbgu_pdf_evo.png')
        plt.close()


def avg_dppdf_oneperiod(output_dir, a_save, fname,
                        B0=1.0,
                        Lx=1.0,
                        prob='shear_alfven'):

    avg = 20
    n_hist, bins = diag.dp_pdf(output_dir, 20, prob=prob)
    for n in range(21, 41):
        n_hist_temp = diag.dp_pdf(output_dir, n, prob=prob)[0]
        n_hist += n_hist_temp
    n_hist /= avg
    plt.semilogy(0.5*(bins[:-1] + bins[1:]), n_hist)
    plt.xlabel(r'$4\pi\Delta p / \langle B^2 \rangle$')
    plt.ylabel(r'$\mathcal{P}$')
    plt.xlim(-1.1, 0.6)
    plt.ylim(1e-5, 1e3)
    plt.grid()
    plt.axvline(-1, ls='--', c='r')  # firehose
    plt.axvline(0.5, ls='--', c='k')  # mirror
    plt.legend()

    plt.savefig(a_save + fname + '_dp_pdf_1tau.pdf')
    plt.savefig(a_save + fname + '_dp_pdf_1tau.png')
    plt.close()


def plot_bbgudist_all_time(output_dir, a_save, fname,
                           prob='shear_alfven'):
    t_dp = diag.get_time_vec(output_dir, prob=prob)  # time vector
    colormap = plt.cm.viridis
    norm = mpl.colors.Normalize(vmin=t_dp[0], vmax=t_dp[-1])

    for n in [5, 20, 40, 80]:  # t/τ_A = 0.25, 1.0, 2.0, 4.0
        print('Doing ' + str(n))
        t = '{:.2f}'.format(t_dp[n])  # scale by Alfven time
        n_hist, bins = diag.bbgu_pdf(output_dir, n, prob=prob)
        plt.semilogy(0.5*(bins[:-1] + bins[1:]), n_hist, color=colormap(norm(t_dp[n])),
                     label=r'$t/\tau_A= \ $' + t)
    plt.xlabel(r'$4\pi\hat{\mathbf{b}} \hat{\mathbf{b}} : \nabla \mathbf{u} / \langle B^2 \rangle$')
    plt.ylabel(r'$\mathcal{P}$')
    plt.xlim(-15, 15)
    plt.ylim(1e-4, 1e1)
    plt.grid()
    plt.legend()

    plt.savefig(a_save + fname + '_bbgu_pdf_evo.pdf')
    plt.savefig(a_save + fname + '_bbgu_pdf_evo.png')
    plt.close()


def do_animation(output_dir, a_save, fname, nu_c, p0, amp,
                 B0=1.0,
                 Lx=1.0,
                 prob='shear_alfven'):
    max_n = diag.get_maxn(output_dir)  # to get the number of frames
    t_A = Lx  # true when v_A = 1 and B0 along the x axis

    beta = diag.beta(p0, B0)
    omega_A = diag.calc_omega_A(Lx)
    db_int = diag.db_int(nu_c, omega_A, beta)
    db = '{:.2f}'.format(amp/db_int)
    regime = diag.REGIMES[diag.find_regime(nu_c, omega_A, beta)]

    t_dp, dp, dp_std = diag.box_averaged_dp(output_dir, prob=prob)  # mean dp

    fig = plt.figure()
    ax = plt.axes(xlim=(-1.1, 0.6), ylim=(1e-5, 1e3))
    ax.set_yscale('log')
    line, = ax.plot([], [], lw=2)
    dp_mean_line = ax.axvline(x=0.0, ls=':', color='g')
    title = ax.set_title('')
    ax.axvline(-1, ls='--', c='r')  # firehose
    ax.axvline(0.5, ls='--', c='k')  # mirror
    ax.grid()

    def init():
        line.set_data([], [])
        dp_mean_line.set_xdata(0.0)
        title.set_text(regime + r', $\delta b / \delta b_\textrm{int} = \ $'
                       + db + r', $t/\tau_A= \ $')
        return line, dp_mean_line, title

    def animate(n):
        print(str(n))  # remove
        data = diag.load_data(output_dir, n, prob)
        t = '{:.2f}'.format(data['Time']/t_A)  # scale by Alfven time

        n_hist, bins = diag.dp_pdf(output_dir, n, prob=prob)
        line.set_data(0.5*(bins[:-1] + bins[1:]), n_hist)
        ax.set(xlabel=r'$4\pi\Delta p / B^2$', ylabel=r'$\mathcal{P}$')
        dp_mean = dp[n]
        dp_mean_line.set_xdata(dp_mean)  # mean dp
        title.set_text(regime + r', $\delta b / \delta b_\textrm{int} = \ $'
                       + db + r', $t/\tau_A= \ $' + t)
        return line, dp_mean_line, title

    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=max_n, interval=80, blit=False)

    # a_save = diag.PATH + a_save + fname
    # a_save += output_dir.split('/')[-1] + '/' + fname
    a_save += fname
    anim.save(a_save + '_dp_pdf.mp4', fps=10,
              extra_args=['-vcodec', 'libx264'])
    print('Saved ' + fname)
    plt.close()

    plot_dp_time_evo(output_dir, a_save, fname, regime, db, t_dp, dp, dp_std, prob=prob)


def plot_dp_time_evo(output_dir, fig_save, fname, regime, db, t, dp, dp_std,
                     prob='shear_alfven'):
    # already loaded in function above, takes ages for large simulations
    # t, dp, dp_std = diag.box_averaged_dp(output_dir, prob=prob)
    plt.plot(t, dp)

    # stop stdev going over limiters
    dp_below = dp - dp_std
    dp_below[dp_below < -1.0] = -1.0
    dp_above = dp + dp_std
    dp_above[dp_above > 0.5] = 0.5

    plt.fill_between(t, dp_below, dp_above, alpha=0.25, color='k')
    plt.axhline(-1.0, ls='--', c='r')
    plt.axhline(0.5, ls='--', c='k')
    plt.xlabel(r'$t/\tau_A$')
    plt.ylabel(r'$\langle 4\pi\Delta p / B^2 \rangle$')
    plt.title(regime + r', $\delta b / \delta b_\textrm{int} = \ $' + db)
    plt.grid()

    fig_save += fname
    fig_save += '_dp_mean'
    plt.savefig(fig_save + '.pdf')
    plt.savefig(fig_save + '.png')
    plt.close()
