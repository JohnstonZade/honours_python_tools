import diagnostics as diag
from spectrum import calc_spectrum
from structure_function import structure_function


# Maybe write a function to generate plot title from fname?
def analyse(fname, plot_title='test'):
    # plot_title  = gen_plot_title(fname)
    fig_path = diag.FIG_PATH + fname
    diag.make_folder(fig_path)

    # Default analyse at (total sim time)/2
    max_n = diag.get_maxn(fname)
    tau = int(max_n/2)

    # Spectrum
    print('Doing Spectrum')
    calc_spectrum(fname, plot_title)
    # PDF
    print('Doing PDF')
    diag.prlshearft_pdf(fname, tau)
    # Structure function
    print('Doing Structure Function')
    structure_function(fname, tau)
    # Energy evolution
    print('Doing Energy Evolution')
    diag.plot_energy_evo(fname, plot_title)
    # Output
