"""
Utility functions to create a tellurim model string using a model defined in 
bcodes as input
"""


# Tellurium model template
model_str_template = """
model {model_name}()

# mass balances
{mass_balances}

# rate equations
{rates}

# functions
{functions}

# parameters
{params}

# initial conditions
{species}
{init}

# events
{events}

end
"""

###############################################################################
# FUNCTIONS

# For rates and functions
def create_func_str(func_dict):
    func_str = ""
    for key in func_dict:
        func_str += "{0}:= {1};\n".format(key, func_dict[key])
        func_str = func_str.replace("**", "^")
    return func_str


# For parameters and init
def create_assignment_str(assign_dict):
    assign_str = ""
    for key in assign_dict:
        assign_str += "{0} = {1};\n".format(key, assign_dict[key])
    return assign_str


def create_massbalance_str(mb_dict):
    mb_str = ""
    for sp in mb_dict:
        s = "{0}' = ".format(sp)
        for rxn in mb_dict[sp]:
            if mb_dict[sp][rxn] > 0:
                sign = "+"
            else:
                sign = "-"
            s += "{0} {1} * {2}".format(sign, abs(mb_dict[sp][rxn]), rxn)
        s += "; \n"
        mb_str += s
    return mb_str


def create_sp_str(id_sp):
    sp_str = "species "
    for sp in id_sp:
        sp_str += "{0}, ".format(sp)
    # remove last comma
    sp_str = sp_str[:-2]
    return sp_str


def create_event_str():
    pass
