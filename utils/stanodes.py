"""
Create odes strings for stan
"""
from bcodes.ratevector import subs_id_by_value


def create_transdict(params2estimate, params_dict, id_sp, params2tune=None):
    """
    Create translation dictionary for string substitution. Parameters to
    estimate are assigned elements of a vector "p". Parameters to tune are
    assigned elements of a vector "x_r" (recognized by Stan). species ids are
    assinged elements of a vector "y". Each of the key-value assignments just
    described are used to update the params_dict, the dictionary of paramter
    values.

    NOTE. vector elements are enumerated from 1 onwards, so that Stan
    understands them.
    Parameters
    ----------
    params2estimate : list
        list of strings. Parameter ids to estimate.
    param_dict : dict
        {parameter_id : parameter_value} Dictionay of paramter values
    id_sp : list
        list of strings with species ids
    params2tune : list, optional
        list of strings. Parameter ids that can be tuned, or given as inputs

    Returns
    -------
     : dict
        {ids: values or vector elements} Dictionary to be used to subsititute
        species and parameter ids for vector elements, or values.
    """
    params = params_dict.copy()

    params2estimate_dict = dict(
        zip(
            params2estimate,
            ["p[{}]".format(i + 1) for i, p in enumerate(params2estimate)],
        )
    )
    params.update(params2estimate_dict)
    # Substitute x_r[] for tunable parameters
    if params2tune:
        params2tune_dict = dict(
            zip(
                params2tune,
                ["x_r[{}]".format(i + 1) for i, p in enumerate(params2tune)],
            )
        )
        params.update(params2tune_dict)

    # Create id_sp substitution dictionary
    id_sp_dict = dict(zip(id_sp, ["y[{}]".format(i + 1) for i, sp in enumerate(id_sp)]))

    return dict(**params, **id_sp_dict)


def create_substituted_rates_dict(rates_dict, trans_dict):
    """
    Substitute the keys in trans_dict for their corresponding values in the
    values (strings) of rates_dict. In addition the python power operator ("**")
    is substituted for the operator used in stan ("^").

    Parameters
    ----------
    rates_dict : dict
        {rxn_ids : rate_law} Dictionary of rate ids and their corresponding rate
        laws as strings.
    trans_dict : dict
        {ids : values or vector elements} Dictionary specifying which strings
        (ids) should be substituted by vector elements or numerical values.

    Returns
    -------
    rates : dict
        {rxn_ds : rate_law} rate ides and their corresponding rate laws as
        strings, with the ids specified in trans_dict substituted by their
        corresponging numerical values of vector elements.
    """
    rates = {}
    for rxn in rates_dict:
        rates[rxn] = subs_id_by_value(rates_dict[rxn], trans_dict)
        rates[rxn] = rates[rxn].replace("**", "^")
    return rates


def create_odes_str(id_sp, id_rs, rates, mass_balances, name="odes"):
    """
    Create odes string for substituting inside a Stan code string.

    Parameters
    ----------
    id_sp : list
        list of strings. species ids
    id_rs : list
        list of strings. reaction ids
    rates : dict
        {rxn_ids : rate_law} Dictionary of rate ids and their corresponding rate
        laws as strings.
    mass_balances : dict of dicts
        {sp_id : {rxn_id: stoichiometric coefficient, ...}}. Specifies the
        mass balances using the stoichiometic coefficients of the reactions.
    name: str, optional
        name to assing to the system of odes fucntion in Stan.

    Returns
    -------
    odes : str
        string encoding the system of ordinary differential equations for a
        kineitc model, formatted as a function for Stan.
    """
    odes = """
    real[] {name}(real t, real[] y, real[] p, real[] x_r, int[] x_i)
        {{
        real dydt[{len_id_sp}];\n""".format(
        len_id_sp=len(id_sp), name=name
    )

    for i, sp in enumerate(id_sp):
        dydt = "\t\tdydt[{}] = ".format(i + 1)  # NOTE stan counts from 1
        for rxn in mass_balances[sp]:
            if mass_balances[sp][rxn] > 0:
                sign = "+"
            else:
                sign = "-"
            dydt += "{0} {1} * {2}".format(
                sign, abs(mass_balances[sp][rxn]), rates[rxn]
            )
        dydt += ";\n"
        odes += dydt
    odes += "\t\treturn dydt;\n"
    odes += "\t}"
    return odes


def create_stan_odes_str(
    id_sp,
    id_rs,
    rates,
    mass_balances,
    params,
    params2estimate,
    params2tune=None,
    name="odes",
):
    """
    Create odes string for substituting inside a Stan. First it creates a
    translation dictionary translation dictionary for string substitutions in
    the values of the "rates" dictionary (which are strings), using
    create_transdict(). Parameters to estimate are assigned elements of a vector
    "p". Parameters to tune are assigned elements of a vector "x_r" (recognized
    by Stan). species ids are assinged elements of a vector "y". Each of the
    key-value assignments just described are used to update the params_dict, the
    dictionary of paramter values. Then creates to odes string using
    create_odes_str()

    NOTE. vector elements are enumerated from 1 onwards, so that Stan
    understands them.

    Parameters
    ----------
    id_sp : list
        list of strings. species ids
    id_rs : list
        list of strings. reaction ids
    rates : dict
        {rxn_ids : rate_law} Dictionary of rate ids and their corresponding rate
        laws as strings.
    mass_balances : dict of dicts
        {sp_id : {rxn_id: stoichiometric coefficient, ...}}. Specifies the
        mass balances using the stoichiometic coefficients of the reactions.
    params : dict
        {parameter_id : parameter_value} Dictionay of paramter values
    name: str, optional
        name to assing to the system of odes fucntion in Stan.
    params2estimate : list
        list of strings. Parameter ids to estimate.
    params2tune : list, optional
        list of strings. Parameter ids that can be tuned, or given as inputs
    name: str, optional
        name to assing to the system of odes fucntion in Stan.

    Returns
    -------
    odes : str
        string encoding the system of ordinary differential equations for a
        kineitc model, formatted as a function for Stan.
    """
    trans_dict = create_transdict(params2estimate, params, id_sp, params2tune)
    rates_ = create_substituted_rates_dict(rates, trans_dict)
    odes_str = create_odes_str(id_sp, id_rs, rates_, mass_balances, name=name)
    return odes_str, trans_dict, rates
