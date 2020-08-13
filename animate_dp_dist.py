import diagnostics as diag
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import rc
rc('text', usetex=True)  # LaTeX labels


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

    dp = diag.box_averaged_dp(output_dir, prob=prob)[1]  # mean dp

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

    a_save += output_dir.split('/')[-1] + '/' + fname
    anim.save(a_save + '_dp_pdf.mp4', fps=10,
              extra_args=['-vcodec', 'libx264'])
    print('Saved ' + fname)
    plt.close()

    plot_dp_time_evo(output_dir, a_save, regime, db, prob=prob)


def plot_dp_time_evo(output_dir, fig_save, regime, db,
                     prob='shear_alfven'):
    t, dp, dp_std = diag.box_averaged_dp(output_dir, prob=prob)
    plt.plot(t, dp)
    plt.fill_between(t, dp - dp_std, dp+dp_std, alpha=0.25, color='k')
    plt.axhline(-1.0, ls='--', c='r')
    plt.axhline(0.5, ls='--', c='k')
    plt.xlabel(r'$t/\tau_A$')
    plt.ylabel(r'$\langle 4\pi\Delta p / B^2 \rangle$')
    plt.title(regime + r', $\delta b / \delta b_\textrm{int} = \ $' + db)
    plt.grid()

    fig_save += '_dp_mean'
    plt.savefig(fig_save + '.pdf')
    plt.savefig(fig_save + '.png')
    plt.close()
