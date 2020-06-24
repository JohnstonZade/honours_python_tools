from numpy import sqrt


def calculate(nu, B_prp, p_0):
    assert nu >= 0, 'nu should be positive'

    # assuming B_0 and omega_A are always set to 1
    beta = 2*p_0  # beta = 2p/(B_0)^2
    db = B_prp    # db = B_prp/B_0

    regime = find_regime(nu, beta)
    is_interrupt = find_interrupt(regime, nu, db, beta)
    interrupt_limit = 'will' if is_interrupt else 'will not'

    # write to text file and save
    f_inputs = 'Inputs: \n Collision Frequency nu_c: ' + str(nu) + '\n B_perp: '\
               + str(B_prp) + '\n Inital Pressure p_0: ' + str(p_0) + '\n \n'
    f_regime = 'Plasma is in the ' + regime + ' regime \n'
    f_interrupt = 'Waves ' + interrupt_limit + ' be interrupted.'
    f_string = f_inputs + f_regime + f_interrupt
    print(f_string)

    # change to fname, implement
    # with open('testfile.txt', 'w') as f:
    #     f.write(f_string)

def find_regime(nu, beta):
    if nu < 1:
        return 'collisionless'
    elif nu < sqrt(beta):
        return 'braginskii with heat fluxes'
    else:
        return 'braginskii'


def find_interrupt(regime, nu, db, beta):
    if regime == 'collisionless':
        return db >= 1/sqrt(beta)
    else:
        return db >= sqrt(nu/beta)


nu = float(input('Collision Frequency: '))
Bprp = float(input('B_perp: '))
p0 = float(input('Initial Pressure: '))
calculate(nu, Bprp, p0)
