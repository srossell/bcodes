from scipy.sparse import lil_matrix

def build_stoichiometry_matrix(id_sp, id_rs, mass_balances):
    """
    Construct a stoichiometry matrix.

    ACCEPTS
    id_sp [list] species identifiers
    id_rs [list] reaction identifiers
    mass_balances [nested dict] {metab: {rxn_id: stoichiometry}}

    RETURNS
    stoichiometry matrix [numpy.ndarray] with rows ordered by id_sp and columns
        ordered by id_rs
    """
    sparse_s = lil_matrix((len(id_sp), len(id_rs)))
    for i, sp in enumerate(id_sp):
        for rxn in mass_balances[sp]:
            ri = id_rs.index(rxn)
            sparse_s[i, ri] = mass_balances[sp][rxn]
    return sparse_s.toarray()

