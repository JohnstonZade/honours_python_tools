import numpy as np
import diagnostics as diag


def calculate(output, fname, nu, B_prp, p_0, Ly, Lx=1.0,
              prob='shear_alfven',
              vel_pert=0):
    assert nu >= 0, 'nu should be non-negative'

    # assuming ω_A is always set to 1
    # and assuming a sinusoidal perturbation in the z-direction to a magnetic
    # field in the xy plane. Only for 2D waves

    # Average of background field initially
    # Just B0 for shear alfven case but will keep for more general turbulence
    B_0 = diag.avg_B(output, 0, 1, prob)
    beta = diag.beta(p_0, B_0)  # beta = 2p/B^2 (Gaussian units)
    db = B_prp/B_0

    omega_A = diag.calc_omega_A(Lx)
    regime = diag.find_regime(nu, omega_A, beta)
    db_int = diag.db_int(nu, omega_A, beta)
    db_int *= 2.5 if vel_pert else 1.0
    interrupt_limit = 'will' if db >= db_int else 'will not'

    # write to text file and save
    f_inputs = 'Inputs: \n Collision Frequency ν_c: ' + str(nu)\
               + '\n Alfven frequency ω_A: ' + str(omega_A)\
               + '\n ν_c / ω_A: ' + str(nu/omega_A) + '\n B_prp: '\
               + str(B_prp) + '\n Inital Pressure p_0: ' + str(p_0) + '\n \n'
    f_beta = 'β = ' + str(beta) + '\nδb = ' + str(db) + '\n√β = ' + str(np.sqrt(beta))\
             + '\nδb_int = ' + str(db_int) + '\nLy / Lx = ' + str(Ly) + '/' + str(Lx) + '\n'
    f_regime = '\nPlasma is in the ' + diag.REGIMES[regime] + ' regime \n'
    f_interrupt = 'Waves ' + interrupt_limit + ' be interrupted.'
    f_string = f_inputs + f_beta + f_regime + f_interrupt
    file_path = diag.PATH + output + '/analysis'
    diag.make_folder(file_path)
    with open(file_path + '/' + fname + '.txt', 'w') as f:
        f.write(f_string)
