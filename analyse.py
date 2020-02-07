import diagnostics as diag
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
