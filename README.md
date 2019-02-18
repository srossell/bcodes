# Write ODEs for biochemical simulations

Functions to build systems of ordinary differential equations for biochemical
simulations. The building blocks for these odes are: a list of reaction
identifiers (id_rs), a list of species identifiers (id_sp), a dictionary of
rate equations (with id_rs as keys), a dictionary of parameter values and a
dictionary of mass balancaes, with (id_sp as keys).


Initial commit from
[moduleForOdeModeling.py](https://bitbucket.org/sergiorossell/odekineticmodels)
with minor modifications.

