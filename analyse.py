import diagnostics as diag
import matplotlib.pyplot as plt
import seaborn as sns
from math import pi
from spectrum import calc_spectrum
from structure_function import structure_function


def analyse(fname, do_title=1):
    if do_title:
        plot_title = gen_plot_title(fname)
    else:
        plot_title = 'test'

    fig_path = diag.FIG_PATH + fname
    diag.make_folder(fig_path)

    # Default analyse at (total sim time)/2 as am focusing on
    # continuously forced turbulence at the moment
    max_n = diag.get_maxn(fname)
    tau = int(max_n/2)

    # Spectrum
    print('Doing Spectrum')
    calc_spectrum(fname, plot_title)
    # PDF
    print('Doing PDF')
    diag.prlshearft_pdf(fname, tau, plot_title)
    diag.prlshear_pdf(fname, tau, plot_title)
    # Structure function
    print('Doing Structure Function')
    structure_function(fname, tau)
    # Energy evolution
    print('Doing Energy Evolution')
    diag.plot_energy_evo(fname, plot_title)
    # Output


def gen_plot_title(fname):
    # do some parsing of file name?
    # or use It_Brag parameter
    return 'It Brag = ' + str(diag.It_Brag(fname))


# This is a mess, only made for quickly generating some plots
# Not intended to be kept
def plot_pdfs():
    fnames = ['cgl_31Jan_320160_nuc1',
              'cgl_31Jan_320160_nuc4',
              'cgl_31Jan_320160_nuc10',
              'cgl_10Feb_320160_nuc1e8']

    ItBrags = [diag.It_Brag(fname) for fname in fnames]
    X, prlshears, B2avgs, Δps = [], [], [], []
    for fname in fnames:
        max_n = diag.get_maxn(fname)
        tau = int(max_n/2)

        prlshear, B2avg = diag.prlshearft_pdf(fname, tau, do_plot=0)
        Δp_bb = diag.prlshear_pdf(fname, tau, do_plot=0)
        N, bins, patches = plt.hist(prlshear, 100, density=1)
        X.append(0.5*(bins[:-1] + bins[1:]))
        prlshears.append(N)
        B2avgs.append(B2avg)
        Δps.append(Δp_bb)

        print(fname + ' done')
    plt.clf()
    for i, prlshear in enumerate(prlshears):
        plt.plot(X[i], prlshear)
        plt.yscale('log')
    for i in range(len(fnames)-1):
        ItBrag = ItBrags[i]
        B2 = B2avgs[i]
        x_low = -(2/pi)*ItBrag*B2
        x_high = (1/pi)*ItBrag*B2
        plt.vlines([x_low, x_high], 1e-4, 1, colors="C{}".format(i),
                   linestyles='dashed')

    ax = plt.gca()
    plt.legend(labels=[gen_plot_title(fname) for fname in fnames])
    plt.title(r'Comparison of $\mathbf{\hat{b}\hat{b}}:\nabla\mathbf{u}$ PDFs at It Brag = 0.8, 3.2, 8, and MHD')
    plt.xlabel(r'$4\pi \mathbf{\hat{b}\hat{b}}:\nabla\mathbf{u}/\langle B^2\rangle$')
    plt.ylabel(r'$\mathcal{P}$')
    ax.set_xlim([-15, 15])
    ax.set_ylim([1e-4, 1])
    plt.show()
    plt.clf()

    # do similar plot as above
    for Δp in Δps:
        sns.kdeplot(Δp, bw=0.5)
    for i in range(len(fnames)-1):
        ItBrag = ItBrags[i]
        B2 = B2avgs[i]
        x_low = -(2/pi)*ItBrag*B2
        x_high = (1/pi)*ItBrag*B2
        plt.vlines([x_low, x_high], 1e-4, 1, colors="C{}".format(i),
                   linestyles='dashed')
    plt.legend(labels=[gen_plot_title(fname) for fname in fnames])
    plt.title(r'Comparison of $\mathbf{\hat{b}\hat{b}}:\nabla\mathbf{u}$ PDFs from $\Delta p$ at It Brag = 0.8, 3.2, 8, and MHD')
    plt.xlabel(r'$(\nu_c / p_0)\Delta p$')
    plt.ylabel(r'$\mathcal{P}$')
    plt.show()
