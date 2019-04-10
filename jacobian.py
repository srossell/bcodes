"""
Build the jacobian for system of odes from metabolic models.
"""

from sympy import Matrix
from sympy import Symbol
import numpy as np

from bcodes.ratevector import subs_id_by_value

def odes_rhs2str(id_sp, id_rs, rates, mass_balances):
    """
    Creates a string representing the rhs of a system of odes

    ACCEPTS
    id_sp [list of str] species ids
    id_rs [list of str] reacion ids
    rates [dict] {rxn_id: rxn_equation}
    mass_balances [dict] {species id: {rxn_id: stoic coefficient}}

    RETURNS
    rhs [str] rhs of the system of odes
    """
    rhs = ''  # rhs of the mass balances
    for sp in id_sp:
        sp_bal = ''
        for rxn in mass_balances[sp]:
            if mass_balances[sp][rxn] > 0:
                sign = '+'
            else:
                sign = '-'
            sp_bal += '{0} {1} * {2}'.format(
                    sign, abs(mass_balances[sp][rxn]), rates[rxn])
        sp_bal += ',\n'
        rhs += sp_bal
    return rhs


def odes_rhs2mat(id_sp, id_rs, rates, mass_balances, param_list):
    """
    Build the code to construct a sympy Matrix for the rhs of a
    system of odes.

    ACCEPTS
    id_sp [list of str] species ids
    id_rs [list of str] reacion ids
    rates [dict] {rxn_id: rxn_equation}
    mass_balances [dict] {species id: {rxn_id: stoic coefficient}}
    param_list[list of str] parameter ids

    RETURNS
    rhs [str] code for a sympy Matrix with the rhs of the system of odes

    """
    rhs = odes_rhs2str(id_sp, id_rs, rates, mass_balances)
    rhs_mat = ''
    for sp in id_sp:
        rhs_mat += "{0} = Symbol('{0}')\n".format(sp)
    for p in param_list:
        rhs_mat += "{0} = Symbol('{0}')\n".format(p)
    rhs_mat += 'rhs = Matrix([{0}])'.format(rhs)
    return rhs_mat


def odes_lhs2mat(id_sp):
    """
    Build the code to construct a sympy martix for the lhs of a system of odes.

    ACCEPTS
    id_sp [list of str] species ids
    id_rs [list of str] reacion ids
    rates [dict] {rxn_id: rxn_equation}
    mass_balances [dict] {species id: {rxn_id: stoic coefficient}}
    param_list[list of str] parameter ids

    RETURNS
    sp_mat [str] code for a sympy Matrix with the rhs of the system of odes.
    """
    sp_str = ''
    for sp in id_sp:
        sp_str += '{0}, '.format(sp)
    sp_str = sp_str[:-1]  # remove last comma
    sp_mat = 'lhs = Matrix([{0}])'.format(sp_str)
    return sp_mat


def create_sym_jac(id_sp, id_rs, rates, mass_balances, param_list):
    """
    Compute symbolically the Jacobian of a system of odes from a metabolic
    model.

    ACCEPTS
    id_sp [list of str] species ids
    id_rs [list of str] reacion ids
    rates [dict] {rxn_id: rxn_equation}
    mass_balances [dict] {species id: {rxn_id: stoic coefficient}}
    param_list[list of str] parameter ids

    RETURNS
    jac [Sympy Matrix] symbolic Jacobian for the system of odes.

    """
    # Getting the necessary sympy classes within the scope of the exec
    # statement.
    d = {'Symbol': Symbol, 'Matrix': Matrix}

    rhs = odes_rhs2mat(id_sp, id_rs, rates, mass_balances, param_list)
    lhs = odes_lhs2mat(id_sp)
    exec(rhs, d)
    exec(lhs, d)
    # Use sympy to calculate the jacobian
    jac = d['rhs'].jacobian(d['lhs'])
    return jac

def func_from_str(func_str):
    return lambda t, y, p: eval(func_str)

def sym_jac2jac_func(sym_jac, trans_dict):
    """
    Convert a symbolic jacobian (sympy.Matrix) into a numpy array function.

    ACCEPTS
    sym_jac [Sympy Matrix] symbolic Jacobian
    trans_dict [dict] dictionary to translate parameters and species into
        vector elements.

    RETURNS
    jac [2d array of lambdas] jacobian function
    """
    j_str = np.array(sym_jac).astype(str).ravel()
    j_str_subs = [
            subs_id_by_value(expression, trans_dict) for expression in j_str]
    jac = lambda t, y, p: np.array([i(t, y, p) for i in
            [
                func_from_str(j_elem)
                for j_elem in j_str_subs
                ]
            ]).reshape(*sym_jac.shape)
    return jac

def build_jacobian(id_sp, id_rs, rates, mass_balances, param_list, trans_dict):
    """
    Build a jacobian function for a metabolic model.

    ACCEPTS
    id_sp [list of str] species ids
    id_rs [list of str] reacion ids
    rates [dict] {rxn_id: rxn_equation}
    mass_balances [dict] {species id: {rxn_id: stoic coefficient}}
    param_list[list of str] parameter ids
    trans_dict [dict] dictionary to translate parameters and species into
        vector elements.

    RETURNS
    jac [2d array of lambdas] jacobian function
    """
    symbolic_jacobian = create_sym_jac(
            id_sp, id_rs, rates, mass_balances, param_list
            )
    jacobian_function = sym_jac2jac_func(symbolic_jacobian, trans_dict)
    return jacobian_function

