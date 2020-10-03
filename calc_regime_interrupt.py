import numpy as np
import diagnostics as diag
default_prob = diag.DEFAULT_PROB


def calculate(output_dir, save_dir, fname, nu_c, dedt, amp, p_0,
              Lx=1.0,
              prob=default_prob,
              vel_pert=0):
    assert nu_c >= 0, 'nu should be non-negative'
    B_0 = 1.0  # diag.avg_B(output_dir, 0, 1, prob)
    beta = diag.beta(p_0, B_0)  # beta = 2p/B^2 (Gaussian units)
    db = amp/B_0
    omega_A = diag.calc_omega_A(Lx)
    regime = diag.find_regime(nu_c, omega_A, beta)
    db_int = diag.db_int(nu_c, omega_A, beta)
    db_int *= 2.5 if vel_pert else 1.0
    interrupt_limit = 'will' if db >= db_int else 'will not'

    # write to text file and save
    f_inputs = 'Inputs: \n Collision Frequency ν_c: ' + str(nu_c)\
               + '\n Alfven frequency ω_A: ' + str(omega_A)\
               + '\n ν_c / ω_A: ' + str(nu_c/omega_A) + '\n dedt: '\
               + str(dedt) + '\n amp = 2√dedt: ' + str(amp) + '\n Inital Pressure p_0: ' + str(p_0) + '\n \n'
    f_beta = 'β = ' + str(beta) + '\nδb = ' + str(db) + '\n√β = ' + str(np.sqrt(beta))\
             + '\nδb_int = ' + str(db_int) + '\n'
    f_regime = '\nPlasma is in the ' + diag.REGIMES[regime] + ' regime \n'
    f_interrupt = 'Waves ' + interrupt_limit + ' be interrupted.'
    f_string = f_inputs + f_beta + f_regime + f_interrupt

    # saving
    if save_dir[-1] != '/':
        save_dir += '/'

    with open(diag.PATH + save_dir + fname + '.txt', 'w') as f:
        f.write(f_string)
    save_to_dictionary(output_dir, save_dir, fname, nu_c, omega_A, amp, p_0, beta, db,
                       db_int, regime)


def save_to_dictionary(output_dir, save_dir, fname, nu_c, ω_A, amp, p_0, beta,
                       db, db_int, regime):
    params = {}
    params['collision_freq'] = nu_c
    params['alfven_freq'] = ω_A
    params['perp_perturbation'] = amp
    params['initial_press'] = p_0
    params['beta'] = beta
    params['db'] = db
    params['db_int'] = db_int
    params['regime'] = regime

    if save_dir[-1] != '/':
        save_dir += '/'

    diag.save_dict(params, save_dir, fname)


def quick_calc(nu_c, p_0, dedt, B_0=1.0, Lx=1.0):
    amp = 2*np.sqrt(dedt)  # via critical balance

    beta = diag.beta(p_0, B_0)  # beta = 2p/B^2 (Gaussian units)
    db = amp/B_0
    omega_A = diag.calc_omega_A(Lx)
    regime = diag.REGIMES[diag.find_regime(nu_c, omega_A, beta)]
    db_int = diag.db_int(nu_c, omega_A, beta)
    interrupt_limit = 'will' if db >= db_int else 'will not'

    print(regime)
    print('β = ' + str(beta))
    print('db = ' + str(db))
    print('db_int = ' + str(db_int))
    print('Waves ' + interrupt_limit + ' be interrupted.')
