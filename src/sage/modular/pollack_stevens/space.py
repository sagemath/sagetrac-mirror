from sage.modules.module import Module
from sage.structure.factory import UniqueFactory
from distributions import Distributions
from sage.modular.dirichlet import DirichletCharacter
from sage.modular.arithgroup.all import Gamma0
from sage.rings.integer import Integer
from sage.rings.rational_field import QQ
from modsym import PSModularSymbolElement, PSModSymAction
from fund_domain import ManinRelations
from sage.rings.padics.precision_error import PrecisionError
from sage.rings.infinity import infinity as oo

class PSModularSymbols_constructor(UniqueFactory):
    def create_key(self, group, weight=None, sign=0, base_ring=None, p=None, prec_cap=None, coefficients=None):
        if isinstance(group, (int, Integer)):
            group = Gamma0(group)
        if base_ring is None and p is None:
            base_ring = QQ
        if coefficients is None:
            if p is not None and prec_cap is None:
                prec_cap = 20
            if isinstance(group, DirichletCharacter):
                character = group.minimize_base_ring()
                group = Gamma0(character.modulus())
                character = (character, None)
            else:
                character = None
            if weight is None: raise ValueError("you must specify a weight or coefficient module")
            k = weight - 2
            coefficients = Distributions(k, p, prec_cap, base_ring, character)
        else:
            # TODO: require other stuff to be None
            pass
        return (group, coefficients, sign)

    def create_object(self, version, key):
        return PSModularSymbolSpace(*key)

PSModularSymbols = PSModularSymbols_constructor('PSModularSymbols')

class PSModularSymbolSpace(Module):
    """
    A class for spaces of modular symbols that use Glenn Stevens'
    conventions.

    There are two main differences between the modular symbols in this
    directory and the ones in sage.modular.modsym.

    - There is a shift in the weight: weight `k=0` here corresponds to
      weight `k=2` there.

    - There is a duality: these modular symbols are functions from
      `Div^0(P^1(\QQ))`, the others are formal linear combinations of
      such elements.

    INPUT:

    - ``V`` -- the coefficient module, which should have a right action of `M_2(\ZZ)`

    - ``domain`` -- a set or None, giving the domain
    """
    Element = PSModularSymbolElement
    def __init__(self, group, coefficients, sign):
        self._group = group
        self._coefficients = coefficients
        self._sign = sign
        self._manin_relations = ManinRelations(group.level()) # should distingish between Gamma0 and Gamma1...
        act = PSModSymAction(self)
        self._populate_coercion_lists_(action_list=[act])

    def source(self):
        return self._manin_relations

    def coefficient_module(self):
        r"""
        Returns the coefficient module of self.

        EXAMPLES:

        ::

        """
        return self._coefficients

    def group(self):
        r"""
        Returns the congruence subgroup of self.

        EXAMPLES:

        ::

        """
        return self._group

    def ngens(self):
        r"""
        Returns the number of generators defining self.

        EXAMPLES:

        ::
        sage: E = EllipticCurve('11a')
        sage: from sage.modular.overconvergent.pollack.modsym_symk import form_modsym_from_elliptic_curve
        sage: phi = form_modsym_from_elliptic_curve(E); phi
        [-1/5, 3/2, -1/2]
        sage: phi.ngens()
        3
        sage: E = EllipticCurve('37a')
        sage: from sage.modular.overconvergent.pollack.modsym_symk import form_modsym_from_elliptic_curve
        sage: phi = form_modsym_from_elliptic_curve(E); phi
        [0, 1, 0, 0, 0, -1, 1, 0, 0]
        sage: phi.ngens()
        9
        """
        return len(self._manin_relations.indices())

    def ncoset_reps(self):
        r"""
        Returns the number of coset representatives defining the full_data of self

        OUTPUT:
        The number of coset representatives stored in the manin relations. (Just the size
        of P^1(Z/NZ))

        EXAMPLES:

        ::

        sage: E = EllipticCurve('11a')
        sage: from sage.modular.overconvergent.pollack.modsym_symk import form_modsym_from_elliptic_curve
        sage: phi = form_modsym_from_elliptic_curve(E); phi
        [-1/5, 3/2, -1/2]
        sage: phi.ncoset_reps()
        12
        sage: E = EllipticCurve('37a')
        sage: phi = form_modsym_from_elliptic_curve(E); phi
        [0, 1, 0, 0, 0, -1, 1, 0, 0]
        sage: phi.ncoset_reps()
        38
        """
        return len(self._manin_relations.reps())

    def level(self):
        r"""
        Returns the level `N` when self is of level `Gamma_0(N)`.

        INPUT:
            none

        OUTPUT:

        The level `N` when self is of level `Gamma_0(N)`

        EXAMPLES:

        ::

        sage: E = EllipticCurve('11a')
        sage: from sage.modular.overconvergent.pollack.modsym_symk import form_modsym_from_elliptic_curve
        sage: phi = form_modsym_from_elliptic_curve(E); phi
        [-1/5, 3/2, -1/2]
        sage: phi.level()
        11
        """

        return self._manin_relations.level()

    def _grab_relations(self):

        # Should this return a dictionary as opposed to a list?
        v = []
        for r in range(len(self._manin_relations.reps())):
            for j in range(self._manin_relations.reps()):
                R = self._manin_relations[j]
                if (len(R) == 1) and (R[0][2] == self._manin_relations.indices(r)):
                    if R[0][0] <> -1 or R[0][1] <> Id:
                        v = v + [R]
        return v

    def zero_elt(self):
        r"""
        Returns the zero element of the space where self takes values.

        INPUT:
            none

        OUTPUT:

        The zero element in the space where self takes values

        EXAMPLES:

        ::

        sage: E = EllipticCurve('11a')
        sage: from sage.modular.overconvergent.pollack.modsym_symk import form_modsym_from_elliptic_curve
        sage: phi = form_modsym_from_elliptic_curve(E); phi
        [-1/5, 3/2, -1/2]
        sage: z=phi.zero_elt(); z
        0
        sage: parent(z)
        <class 'sage.modular.overconvergent.pollack.symk.symk'>
        sage: z.weight
        0
        sage: z.poly
        0
        """
        return self.coefficient_module().zero()

    def zero(self):
        r"""
        Returns the modular symbol all of whose values are zero.

        INPUT:
            none

        OUTPUT:

        The zero modular symbol of self.

        EXAMPLES:

        ::

        sage: E = EllipticCurve('11a')
        sage: from sage.modular.overconvergent.pollack.modsym_symk import form_modsym_from_elliptic_curve
        sage: phi = form_modsym_from_elliptic_curve(E); phi
        [-1/5, 3/2, -1/2]
        sage: zero_sym = phi.zero(); zero_sym
        [0, 0, 0]
        sage: parent(zero_sym)
        <class 'sage.modular.overconvergent.pollack.modsym_symk.modsym_symk'>

        """

        dd = {}
        for rep in self._manin_relations.reps():
            dd[rep] = self.zero_elt()
        #v = [self.zero_elt() for i in range(0, self.ngens())]
        #return C(v, self._manin_relations)
        return PSModularSymbolElement(dd, self)

    def precision_cap(self):
        r"""
        Returns the number of moments of each value of self

        EXAMPLES:

        ::

        """
        return self.coefficient_module()._prec_cap

    def weight(self):
        r"""
        Returns the weight of self

        EXAMPLES:

        ::

        """
        return self.coefficient_module()._k

    def prime(self):
        r"""
        Returns the prime of self

        EXAMPLES:

        ::

        """
        return self.coefficient_module()._p

    def random_element(self, M):
        r"""
        Returns a random OMS with tame level `N`, prime `p`, weight `k`, and
        `M` moments --- requires no `2` or `3`-torsion
        INPUT:

        - M: the number of moments

        OUTPUT:

        - A random element of self

        EXAMPLES:

        ::


        """

        if M > self.precision_cap():
            raise PrecisionError ("Too many moments for self.")

        # Somebody who understands the details should take a careful look.
        N = self.level()
        p = self.prime()
        manin = ManinRelations(N * p)
        dd = {}
        for g in manin.gens()[1:]:
            mu = self.coefficient_module().random_element(M)
            if g in manin.reps_with_two_torsion:
                rj = manin.twotor_index[g]
                gam = manin.twotorrels[rj]
                dd[g] = mu.act_right(gam)._sub_(mu)._lmul_(ZZ(1) / ZZ(2))
            else:
                dd[g] = mu
        t = self.zero()
        for j in range(2, len(manin.coset_relations())):
            R = manin.coset_relations(j)
            if len(R) == 1:
                if R[0][0] == 1:
                    rj = manin.gens(j - 1)
                    t = t + dd[rj]
                    # Should t do something?
                else:
                    index = R[0][2]
                    rj = manin.gens(index - 1)
                    mu = dd[rj]
                    t = t + mu.act_right(R[0][1])._lmul_(R[0][0])
                    # Should t do something?
        dd[manin.gens(0)] =  mu._lmul_(-1)
        return modsym_dist(dd, self)

def cusps_from_mat(g):
    a, b, c, d = g.list()
    if c: ac = a/c
    else: ac = oo
    if d: bd = b/d
    else: bd = oo
    return ac, bd

def form_modsym_from_elliptic_curve(E):
    N = E.conductor()
    V = PSModularSymbols(Gamma0(N), 2)
    D = V.coefficient_module()
    manin = V.source()
    plus_sym = E.modular_symbol(sign = 1)
    minus_sym = E.modular_symbol(sign = -1)
    val = {}
    for g in manin.gens():
        ac, bd = cusps_from_mat(g)
        val[g] = D([plus_sym(ac) + minus_sym(ac) - plus_sym(bd) - plus_sym(bd)])
    return V(val)

def form_modsym_from_decomposition(A):
    """
    """
    M = A.ambient_module()
    w = A.dual_eigenvector()
    K = w.base_ring()
    V = PSModularSymbols(A.group(), A.weight(), base=K) # should eventually add sign as well.
    D = V.coefficient_module()
    N = V.level()
    k = V.weight() # = A.weight() - 2
    manin = V.source()
    val = {}
    for g in manin.gens():
        ac, bd = cusps_from_mat(g)
        v = []
        for j in range(k+1):
            # The following might be backward: it should be the coefficient of X^j Y^(k-j)
            v.append(w.dot_product(M.modular_symbol([j, ac, bd]).element()) * binomial(k, j))
        val[g] = D(v)
    return V(val)
