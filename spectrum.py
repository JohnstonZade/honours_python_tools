'''Code to calculate the energy spectrum of turbulence
   in a fluid simulation.
'''
import numpy as np
import numpy.fft as fft
import matplotlib.pyplot as plt
from diagnostics import ft_grid, load_data, load_hst  # , plot_energy_evo


# TODO: work on saving/loading dictionary
def calc_spectrum(fname, plot_title='test', do_mhd=1, do_full_calc=1):

    # Getting turnover time and converting to file number
    tau_file, nums = get_turnover_time(fname, 0)

    if do_full_calc:
        # create grid of K from first time step

        data = load_data(fname, tau_file)
        (KX, KY, KZ), kgrid = ft_grid(data, 1)
        Kprl = np.abs(KX)
        Kperp = np.sqrt(np.abs(KY)**2 + np.abs(KZ)**2)
        Kmag = np.sqrt(Kprl**2+Kperp**2)
        Kspec = Kmag

        # Dictionary to hold spectrum information
        S = {}
        S['Nk'] = len(kgrid) - 1
        S['kgrid'] = 0.5*(kgrid[:-1] + kgrid[1:])

        # Count the number of modes in each bin to normalize later -- this
        # gives a smoother result, as we want the average energy in each bin.
        # TODO: check where this is used
        oneGrid = np.ones(KX.shape)
        S['nbin'] = spect1D(oneGrid, oneGrid, Kspec, kgrid)*np.size(oneGrid)**2
        S['nnorm'] = S['nbin']/S['kgrid']**2
        S['nnorm'] /= np.mean(S['nnorm'])

        # average over snapshots in nums
        def m3(a):
            return np.mean(np.mean(np.mean(a)))

        ns = 0
        fields = ['vel1', 'vel2', 'vel3', 'Bcc1', 'Bcc2', 'Bcc3',
                  'EK', 'EM', 'B', 'rho']

        # Initializing variable fields in spectrum dictionary
        for var in fields:
            S[var] = 0

        for n in nums:
            try:
                data = load_data(fname, n)
            except IOError:
                print('Could not load file', n)
                break

            print('Doing n =', n)

            for vel in fields[:3]:
                ft = fft.fftn(data[vel])
                S[vel] += spect1D(ft, ft, Kspec, kgrid)
                S['EK'] += S[vel]  # Total spectrum is sum of each component

            if do_mhd:
                Bmag = 0
                for Bcc in fields[3:6]:
                    ft = fft.fftn(data[Bcc])
                    S[Bcc] += spect1D(ft, ft, Kspec, kgrid)
                    S['EM'] += S[Bcc]
                    Bmag += data[Bcc]**2

                Bmag = np.sqrt(Bmag)
                ft_Bmag = fft.fftn(Bmag)
                S['B'] += spect1D(ft_Bmag, ft_Bmag, Kspec, kgrid)

            ft_rho = fft.fftn(data['rho'] - m3(data['rho']))
            S['rho'] += spect1D(ft_rho, ft_rho, Kspec, kgrid)

            ns += 1

        # Average over times done
        for var in fields:
            S[var] /= ns
        S['nums'] = nums

        # save spectrum data somehow
        # try using pickle
    else:
        # load spectrum data
        return 0
    plot_spectrum(S, plot_title, do_mhd)
    # plot_energy_evo(fname, plot_title)


# TODO: work on saving figure
def plot_spectrum(S, plot_title, do_mhd=1):
    # plot spectrum
    if do_mhd:
        plt.loglog(S['kgrid'], S['EK'], S['kgrid'], S['B'],
                   S['kgrid'], S['kgrid']**(-5/3), ':')
        plt.legend([r'$E_K$', r'$E_B$', r'$k^{-5/3}$'])
    else:
        plt.loglog(S['kgrid'], S['EK'], S['kgrid'], S['kgrid']**(-5/3), ':')
        plt.legend([r'$E_K$', r'$k^{-5/3}$'])
    plt.xlabel(r'$k$')
    plt.ylabel(r'$E_K$')
    plt.title('Energy  Spectrum: ' + plot_title)
    plt.show()
    # save figure


def spect1D(v1, v2, K, kgrid):
    '''Function to find the spectrum < v1 v2 >,
    K is the kgrid associated with v1 and v2
    kgrid is the grid for spectral shell binning
    '''
    nk = len(kgrid) - 1
    out = np.zeros((nk, 1))
    NT2 = np.size(K)**2
    for k in range(nk):
        mask = np.logical_and(K < kgrid[k+1], K > kgrid[k])
        sum = np.sum(np.real(v1[mask])*np.conj(v2[mask]))
        out[k] = np.real(sum) / NT2
    return out


# TODO: work on making more general, check against MATLAB
def get_turnover_time(fname, do_decay):
    hstData = load_hst(fname)
    KEx, KEy, KEz = hstData['1-KE'], hstData['2-KE'], hstData['3-KE']

    u_L = 0.
    if do_decay:
        u_L = np.sqrt(2*(KEx[1]+KEy[1]+KEz[1]))
    else:
        n = 0
        for f in range(100, 3001):
            try:
                u_L += np.sqrt(2*(KEx[f]+KEy[f]+KEz[f]))
                n += 1
            except IndexError:
                print('No entry at', f)
                break
            u_L /= n
    u_L
    tau = 1 / u_L
    tau_file = int(tau)

    tau_file = 10

    if do_decay:
        nums = range(tau_file, 301)
    else:
        nums = range(tau_file, 6*tau_file+1)


    return tau_file, nums
