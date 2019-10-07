"""
Create odes strings for stan
"""
from bcodes.ratevector import subs_id_by_value

def create_transdict(params2estimate, params_dict, id_sp):
    """
    """
    params = params_dict.copy()

    # NOTE stan starts counting form 1 not zero
    # Substitute p[] for parameters to estimate
    params2estimate_dict = dict(zip(
        params2estimate,
        ['p[{}]'.format(i + 1) for i, p in enumerate(params2estimate)]
        ))
    params.update(params2estimate_dict)

    # Create id_sp substitution dictionary
    id_sp_dict = dict(zip(
        id_sp, ['y[{}]'.format(i + 1) for i, sp in enumerate(id_sp)]
        ))

    return dict(**params, **id_sp_dict)

def create_substituted_rates_dict(rates_dict, trans_dict):
    rates = {}
    for rxn in rates_dict:
        rates[rxn] = subs_id_by_value(rates_dict[rxn], trans_dict)
    return rates

def create_odes_str(id_sp, id_rs, rates, mass_balances):
    odes = """
    real[] odes(real t, real[] y, real[] p, real[] x_r, int[] x_i)
        {{
        real dydt[{}];\n""".format(len(id_sp))

    for i, sp in enumerate(id_sp):
        dydt = '\t\tdydt[{}] = '.format(i + 1)  # NOTE stan counts from 1
        for rxn in mass_balances[sp]:
            if mass_balances[sp][rxn] > 0:
                sign = '+'
            else:
                sign = '-'
            dydt += '{0} {1} * {2}'.format(
                    sign, abs(mass_balances[sp][rxn]), rates[rxn])
        dydt += ';\n'
        odes += dydt
    odes += '\t\treturn dydt;\n'
    odes += '\t}'
    return odes

def create_stan_odes_str(
        id_sp, id_rs, rates, mass_balances, params, params2estimate
        ):
    """
    Creates a string describing a system of differential equations from bcodes,
    to be used in stan.
    ACCEPTS
    id_sp [list of str] species ids
    id_rs [list of str] reacion ids
    rates [dict] {rxn_id: rxn_equation}
    mass_balances [dict] {species id: {rxn_id: stoic coefficient}}
    params [dict} {parameter id: numerical value}
    params2estimate [list] ids of the parameters to estimate
    """
    trans_dict = create_transdict(params2estimate, params, id_sp)
    rates_ = create_substituted_rates_dict(rates, trans_dict)
    odes_str = create_odes_str(id_sp, id_rs, rates_, mass_balances)
    return odes_str, trans_dict, rates


