import numpy as np
import diagnostics as diag


def calculate(output, fname, nu, B_prp, p_0):
    assert nu >= 0, 'nu should be non-negative'

    # assuming B_0 and omega_A are always set to 1
    # and assuming a sinusoidal perturbation along x
    # only for 2D CPAW Waves
    data = diag.load_data(output, 0, 'cpaw')
    B_2_avg = avg_B2(data)
    beta = 2*p_0 / B_2_avg  # beta = 2p/B^2
    db = B_prp    # db = B_prp/B_0

    regime = find_regime(nu, beta)
    is_interrupt = find_interrupt(regime, nu, db, beta)
    interrupt = str(1/np.sqrt(beta)) if regime=='collisionless' else str(np.sqrt(nu/beta))
    interrupt_limit = 'will' if is_interrupt else 'will not'

    # write to text file and save
    f_inputs = 'Inputs: \n Collision Frequency nu_c: ' + str(nu) + '\n B_perp: '\
               + str(B_prp) + '\n Inital Pressure p_0: ' + str(p_0) + '\n \n'
    f_beta = 'β = ' + str(beta) + '\nδb = ' + str(db) + '\n√β = ' + str(np.sqrt(beta))\
             + '\nδb_max = ' + interrupt + '\n'
    f_regime = '\nPlasma is in the ' + regime + ' regime \n'
    f_interrupt = 'Waves ' + interrupt_limit + ' be interrupted.'
    f_string = f_inputs + f_beta + f_regime + f_interrupt

    # change to fname, implement
    with open(diag.PATH + output + '/' + fname + '.txt', 'w') as f:
        f.write(f_string)


def find_regime(nu, beta):
    if nu < 1:
        return 'collisionless'
    elif nu < np.sqrt(beta):
        return 'braginskii with heat fluxes'
    else:
        return 'braginskii'


def find_interrupt(regime, nu, db, beta):
    if regime == 'collisionless':
        return db >= 1/np.sqrt(beta)
    else:
        return db >= np.sqrt(nu/beta)


def avg_B2(data):
    def pnts(p):
        points = []
        for zz in range(p[0]):
            for yy in range(p[1]):
                for xx in range(p[2]):
                    points.append((zz, yy, xx))
        return np.array(points)

    grid = data['RootGridSize'][::-1]
    points = pnts(grid)

    B_data = (data['Bcc1'], data['Bcc2'], data['Bcc3'])
    B_vec = np.array([diag.get_vec(B_data, p) for p in points])
    B2_avg = np.mean(diag.get_mag(B_vec)**2)  # spatial average of B^2
    return B2_avg
