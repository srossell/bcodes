"""
Functions to build rate vectors. Rate vectors are list of lambda fucntions that
accept species values and, optionally, time and parameters as inputs.
"""

import re

def subs_id_by_value(rate_eq, trans_dict):
    """
    Takes a rate equation (rate_eq), which is a string, and substitutes the
    keys in trans_dict for their values.
    ACCEPTS
    rate_eq [str] rate equation
    trans_dict [dict] {key : value}

    RETURNS
    s [str] rate_eq with keys substituted by values
    """
    s = rate_eq[:]
    for key in trans_dict:
        if type(trans_dict[key]) != str:
            s = re.sub(r'\b%s\b' %  key, str(trans_dict[key]), s)
        else:
            s = re.sub(r'\b%s\b' %  key, trans_dict[key], s)
    return s

def rate_func_from_eq(rate_eq, trans_dict, eval_scope=None, t_is_arg=False,
        p_is_arg=False):
    """
    Takes a rate equation as a string (rate_eq) and substitutes parameter
    and species names by their values or vector variables in trans_dict
    and creates a lambda function that dependens on the vector variables.

    Uses the function subs_id_by_value to substitute the keys in trans_dict
    for their values.

    ACCEPTS
    rate_eq [str] rate equation
    trans_dict [dict] {paramId : paramValue, speciesId : x[..]}
    eval_scope[dict] Dictionary defining the scope where eval (below) functions
    t_is_input [bool] whether t (e.g. time) may be an argument to the function
    p_is_input [bool] whether p (e.g. parameter vector to fit) may be an
        argument of the fuction returned

    RETURNS [lambda] function that depends on the vector x
    """
    if not eval_scope:
        eval_scope = {}
    if t_is_arg and p_is_arg:
        return lambda x, t, p: eval(subs_id_by_value(rate_eq, trans_dict),
                dict({'x' : x, 't': t, 'p': p}, **eval_scope))
    elif t_is_arg and not p_is_arg:
        return lambda x, t: eval(subs_id_by_value(rate_eq, trans_dict),
                dict({'x' : x, 't': t}, **eval_scope))
    elif p_is_arg and not t_is_arg:
        return lambda x, p: eval(subs_id_by_value(rate_eq, trans_dict),
                dict({'x' : x, 'p': p}, **eval_scope))
    else:
        return lambda x: eval(subs_id_by_value(rate_eq, trans_dict),
                                dict({'x' : x}, **eval_scope))



def create_rate_vector(id_sp, id_rs, rate_eq_dict, param_dict, eval_scope=None,
        t_is_arg=False, p_is_arg=False):
    """
    Creates a rate vector v = f(x) where x is the concentration vector
    of species. The vector conserves the order in id_rs

    Uses the function rate_func_from_eq to produce a lambda function
    for each rate equation if rate_eq_dict

    ACCEPTS
    id_sp [list] ordered species ids
    id_rs [list] ordered reaction ids
    rate_eq_dict [dict] {rxnId : rate equation (as a string)
    param_dict [dict] {paramId : paramValue}
    eval_scope[dict] Dictionary defining the scope where the 'eval' statement,
        within rate_func_from_eq, functions.
    t_is_input [bool] whether t (e.g. time) may be an argument to the function
    p_is_input [bool] whether p (e.g. parameter vector to fit) may be an
        argument of the fuction returned

    RETURNS [list of lambda functions] v = f(x)
    """
    # species names to concentration vector 'x'
    id_sp_dict = dict(zip(id_sp, ['x[%i]' % i for i, sp in enumerate(id_sp)]))
    subs_dict = param_dict.copy()
    subs_dict.update(id_sp_dict)
    if t_is_arg and p_is_arg:
        return lambda x, t, p: [i(x, t, p) for i in [
            rate_func_from_eq(rate_eq_dict[r], subs_dict, eval_scope, t_is_arg,
                p_is_arg) for r in id_rs]]
    elif t_is_arg and not p_is_arg:
        return lambda x, t: [i(x, t) for i in [
            rate_func_from_eq(rate_eq_dict[r], subs_dict, eval_scope, t_is_arg,
                p_is_arg) for r in id_rs]]
    elif p_is_arg and not t_is_arg:
        return lambda x, p: [i(x, p) for i in [
            rate_func_from_eq(rate_eq_dict[r], subs_dict, eval_scope, t_is_arg,
                p_is_arg) for r in id_rs]]
    else:
        return lambda x: [i(x) for i in [
            rate_func_from_eq(rate_eq_dict[r], subs_dict, eval_scope, t_is_arg,
                p_is_arg) for r in id_rs]]
