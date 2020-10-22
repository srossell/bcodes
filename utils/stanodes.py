"""
Create odes strings for stan
"""
import csv

import numpy as np
import pandas as pd

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


def create_odes_str(id_sp, id_rs, rates, mass_balances, name):
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
    name: str
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


# For reading optimization resutls
def decoment(file_io):
    """
    Reads a file and skips files starting with "#"

    Parameters
    ----------
    file_io : _ioTextIOWraper
        obtained when using open(myfile.txt)
    """
    for row in file_io:
        raw = row.split("#")[0].strip()
        if raw:
            yield raw


def read_csv_ignoring_comments(myfile):
    """
    reads a csv file ignoring rows starting with "#".

    Parameters
    ----------
    myfile : str
        path to csv file

    Returns
    -------
    df : pd.DataFrame
        csv table as a pandas dataframe.
    """
    with open(myfile, "r") as f:
        reader = csv.reader(decoment(f))
        table = [row for row in reader]
    df = pd.DataFrame(np.array(table).T)
    return df

def read_stanoptimize_result(opt_output_file, params2estimate, t_obs, id_sp, name):
    """
    Build estimated parameters dictionary and fitted species in time dataframe
    from the output of stan's optimize (a csv file).

    Parameters
    ----------
    opt_output_file : str
        path to stan's optimize csv file
    params2estimate : list of str
        id of the parameters estimated
    t_obs : list or 1d-array of floats
        times to which the paramters were fitted
    id_sp : list of str
        id of the species fitted
    name : str
        name given to the model (and parameters) in the stan model.

    Returns
    -------
    opt_p_dict : dict {parameter_id: value}
        Optimized parameters dictionary
    y_hat : pd.DataFrame
        Fitted species in time
    """
    opt_df = read_csv_ignoring_comments(opt_output_file)

    # Get optimized parameters
    p_mask = opt_df[0].str.contains(f"p_{name}")
    opt_p_dict = dict(zip(params2estimate, opt_df.loc[p_mask, 1].astype(float)))

    #########
    # Get y_hat
    yh_mask = opt_df[0].str.contains(f"y_hat_{name}")

    # build indexes for time and species from Stan's result name
    # y_hat_{name}_{index_time}_{index_species}
    yh_indexes = (
        opt_df[yh_mask][0]
        .str.split(".", expand=True)[[1, 2]]
        .rename(columns={1: "ind_t", 2: "ind_sp"})
    )
    # Get fitted species in time values
    y_hat_long = pd.merge(
        yh_indexes, opt_df[yh_mask][[1]], left_index=True, right_index=True
    ).rename(columns={1: "value"})

    y_hat_long["sp"] = (
        y_hat_long["ind_sp"]
        .astype(int)
        .map(dict(zip(range(1, len(id_sp) + 1), id_sp)))
    )
    y_hat_long["t"] = (
        y_hat_long["ind_t"].astype(int).map(dict(zip(range(1, len(t_obs) + 1), t_obs)))
    )

    y_hat = y_hat_long.pivot(columns="sp", values="value", index="t").astype(float)

    return opt_p_dict, y_hat

