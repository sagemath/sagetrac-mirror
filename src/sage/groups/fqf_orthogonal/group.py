from sage.groups.abelian_gps.abelian_aut import AbelianGroupAutomorphismGroup_subgroup
from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
from sage.groups.fqf_orthogonal.lift import _gens, _reverse_homogeneous_blocks
from sage.misc.cachefunc import cached_method
from sage.rings.all import mod, ZZ, IntegerModRing, Zp
from sage.matrix.all import matrix

class FqfOrthogonalGroup_ambient(AbelianGroupAutomorphismGroup_subgroup):
    r"""
    INPUT:

        - T a non degenerate torsion quadratic module.
    """
    def __init__(self, T):
        r"""
        """
        # TODO: check normal form

        invs = T.invariants()
        A = AbelianGroupGap(invs)
        ambient = A.aut()
        gens = _compute_gens(T)
        gens = [ambient(g) for g in gens]
        AbelianGroupAutomorphismGroup_subgroup.__init__(self, ambient, gens)

def _compute_gens(T):
    r"""
    INPUT:

    - ``T`` -- torsion orthogonal module in normal form.
    """
    # normal form gens
    invs = T.invariants()
    if len(invs) == 0:
        return []
    blocks = []
    gensT_orders = [t.order() for t in T.gens()]
    n = len(T.gens())
    P = T.annihilator().gen().prime_divisors()
    for p in P:
        indices = []
        for k in range(len(gensT_orders)):
            if mod(gensT_orders[k], p) == 0:
                indices.append(k)
        blocks.append([p,tuple(indices)])

    # transformation matrices between generators
    to_smith = matrix([t.vector() for t in T.gens()])
    R = IntegerModRing(invs[-1])
    E = matrix.identity(R,len(invs))
    # view T as a submodule of (ZZ/nZZ)^(len(invs))
    B = to_smith.change_ring(R)
    for k in range(B.ncols()):
        B[:,k] *= invs[-1]//invs[k]
        E[:,k] *= invs[-1]//invs[k]
    # unfortunatly solve_left does not work with matrices
    to_normal = matrix([B.solve_left(e) for e in E.rows()])
    to_normal = to_normal.change_ring(ZZ)

    # check that we did not mess up
    tmp = to_normal*to_smith
    E = matrix.identity(R,len(invs))
    for k in range(E.ncols()):
        assert 0 == (tmp[:,k]-E[:,k]).change_ring(ZZ).mod(invs[k])

    # compute generators of the orthogonal groups
    gens = []
    for bl in blocks:
        # compute the generators of the p-primary part
        # the whole group is the direct product of the p-primary parts
        p = bl[0]
        indices = bl[1]
        q_p = T.gram_matrix_quadratic()[indices,indices]
        b = invs[-1].valuation(p)
        G_p = q_p*(p**b)
        R = Zp(p, type='fixed-mod', prec=b+5, show_prec=False)
        G_p = G_p.change_ring(R)

        # there is an unfortunate convention in _gens
        # the blocks of G are descending
        # our blocks are ascendings
        U = _reverse_homogeneous_blocks(G_p)
        Uinv = (~U).change_ring(ZZ)
        gens_mat = _gens(U.T*G_p*U, b)
        gens_mat = [U * g * Uinv for g in gens_mat]

        # extend as identity on the orthogonal complement
        E1 = matrix.identity(indices[0])
        E2 = matrix.identity(n - indices[-1] - 1)
        for g in gens_mat:
            g = g.change_ring(ZZ)
            g = matrix.block_diagonal([E1,g,E2])
            g = to_normal * g * to_smith
            gens.append(g)
    return gens

def _kernel_order(G):
    pass
