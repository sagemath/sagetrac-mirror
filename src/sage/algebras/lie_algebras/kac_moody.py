"""
Kac-Moody Algebras

AUTHORS:

- Travis Scrimshaw (2013-05-03): Initial version
"""

#*****************************************************************************
#  Copyright (C) 2013 Travis Scrimshaw <tscrim@ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty
#    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#  See the GNU General Public License for more details; the full text
#  is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element import RingElement
from sage.categories.algebras import Algebras
from sage.categories.lie_algebras import LieAlgebras
from sage.misc.cachefunc import cached_method
from sage.misc.indexed_generators import IndexedGenerators
from sage.misc.misc import repr_lincomb

from sage.algebras.free_algebra import FreeAlgebra
from sage.algebras.lie_algebras.lie_algebra import FinitelyGeneratedLieAlgebra
from sage.algebras.lie_algebras.finitely_presented import FinitelyPresentedLieAlgebra
from sage.algebras.lie_algebras.lie_algebra_element import LieAlgebraElement
from sage.algebras.lie_algebras.free_lie_algebra import FreeLieAlgebra
from sage.combinat.root_system.cartan_matrix import CartanMatrix
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.free_module import CombinatorialFreeModule
from sage.rings.arith import binomial
from sage.sets.family import Family

class KacMoodyAlgebra(FinitelyPresentedLieAlgebra):
    r"""
    A Kac-Moody algebra.

    A Kac-Moody algebra over a ring `R` is given by the following data,
    known as Cartan datum:

    - a `n \times n` (generalized) Cartan matrix `A = (a_{ij})` of rank `r`,

    - a vector space `\mathfrak{h}` over `R` of dimension `2n - r`,

    - a set of `n` linearly independent simple coroots `\alpha_i^{\vee} \in
      \mathfrak{h}` and roots `\alpha_i \in \mathfrak{h}^*` such that
      `\alpha_i(\alpha_i^{\vee}) = a_{ij}`.

    The Kac-Moody algerba is the Lie algebra generated by `e_i, f_i`, and
    `h \in \mathfrak{h}` which satisfy the following relations:

    .. MATH::

        \begin{aligned}
        [h, h^{\prime}] &= 0,
        \\ [h, e_i] &= \alpha_i(h) e_i,
        \\ [h, f_i] &= -\alpha_i(h) f_i,
        \\ [e_i, f_j] &= \delta_{ij} \alpha_i^{\vee},
        \\ \mathrm{ad}_{e_i}^{1-a_{ij}} e_j &= 0,
        \\ \mathrm{ad}_{f_i}^{1-a_{ij}} f_j &= 0.
        \end{aligned}

    where `\mathrm{ad}_x y = [x, y]`, the adjoint operator of `\mathfrak{g}`.
    """
    @staticmethod
    def __classcall_private__(cls, R, cm, index_set=None):
        """
        Normalize arguments to ensure a unique representation.

        TESTS::

            sage: KM1 = KacMoodyAlgebra(QQ, ['A', 3])
            sage: KM2 = KacMoodyAlgebra(QQ, 'A3')
            sage: KM3 = KacMoodyAlgebra(QQ, CartanType(['A', 3]))
            sage: KM1 is KM2 and KM2 is KM3
            True
        """
        cm = CartanMatrix(cm, index_set=index_set) # index_set not currently supported
        return super(KacMoodyAlgebra, cls).__classcall__(cls, R, cm)

    def __init__(self, R, cm):
        """
        Initialize ``self``.

        TESTS::

            sage: TestSuite(KacMoodyAlgebra(QQ, ['B', 3])).run()
            sage: TestSuite(KacMoodyAlgebra(QQ, ['G', 2])).run()
            sage: TestSuite(KacMoodyAlgebra(QQ, ['A', 2, 1])).run()
            sage: TestSuite(KacMoodyAlgebra(QQ, ['E', 6, 2])).run()
        """
        index_set = cm.index_set()
        n = len(index_set)
        self._n = n
        self._cartan_matrix = cm

        # Construct the base free Lie algebra
        names  = ['e%s'%i for i in index_set]
        names += ['f%s'%i for i in index_set]
        names += ['h%s'%i for i in index_set]
        names += ['d%s'%i for i in range(n-cm.rank())]
        F = FreeLieAlgebra(R, names).Lyndon()
        e = F.gens()[:n]
        f = F.gens()[n:2*n]
        h = F.gens()[2*n:3*n]
        d = F.gens()[3*n:]

        # Construct the basic relations
        P_check = h+d
        rels = [F.bracket(x, y) for i,x in enumerate(P_check) for y in P_check[i+1:]]
        rels += [F.bracket(x, y) - cm[i,j]*x for i,x in enumerate(e) for j,y in enumerate(h)]
        rels += [F.bracket(x, y) for x in e for y in d]
        rels += [F.bracket(x, y) + cm[i,j]*x for i,x in enumerate(f) for j,y in enumerate(h)]
        rels += [F.bracket(x, y) for x in f for y in d]
        for i,x in enumerate(e):
            for j,y in enumerate(f):
                if i == j:
                    rels.append(F.bracket(x, y) - h[i])
                else:
                    rels.append(F.bracket(x, y))
        # Construct the Serre relations
        for i in range(n):
            for j in range(n):
                if i != j:
                    cur_e = e[i]
                    cur_f = f[i]
                    for k in range(1-cm[i,j]):
                        cur_e = F.bracket(e[j], cur_e)
                        cur_f = F.bracket(f[j], cur_f)
                    rels.append(cur_e)
                    rels.append(cur_f)

        category = LieAlgebras(R).WithBasis()
        if cm.is_finite():
            category = category.FiniteDimensional()
        FinitelyPresentedLieAlgebra.__init__(self, F, rels, category=category)

        # Define the actual generators
        elt = lambda x: self.element_class(self, x)
        self._e = map(elt, e)
        self._f = map(elt, f)
        self._h = map(elt, h)
        self._d = map(elt, d)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: KacMoodyAlgebra(QQ, ['B',2])
            Kac-Moody algebra of type ['B', 2]
            sage: cm = CartanMatrix([[2,-5],[-4,2]])
            sage: KacMoodyAlgebra(QQ, cm)
            Kac-Moody algebra with Cartan matrix:
            [ 2 -5]
            [-4  2]
        """
        ct = self._cartan_matrix.cartan_type()
        if ct is not self._cartan_matrix:
            return "Kac-Moody algebra of type {}".format(ct)
        return "Kac-Moody algebra with Cartan matrix:\n{}".format(self._cartan_matrix)

    def _basis_cmp(self, x, y):
        """
        Return the ``cmp`` of two basis elements ``x`` and ``y``.
        """
        return cmp(x,y)
        #return cmp(x.value, y.value)

    def _construct_UEA(self):
        """
        Return the universal enveloping algebra of ``self``.

        EXAMPLES::

            sage: K = KacMoodyAlgebra(QQ, ['B',2])
            sage: UEA = K.universal_enveloping_algebra() # indirect doctest
            sage: UEA
            Poincare-Birkhoff-Witt basis corresponding to Kac-Moody algebra of type ['B', 2]
            sage: TestSuite(UEA).run()
        """
        indices = map(lambda x: x.value.leading_support(), self.basis())
        return self.pbw_basis(indices)

    def cartan_type(self):
        """
        Return the Cartan type of ``self``.

        EXAMPLES::

            sage: KM = KacMoodyAlgebra(QQ, ['B',2])
            sage: KM.cartan_type()
            ['B', 2]
            sage: cm = CartanMatrix([[2,-5],[-4,2]])
            sage: KM = KacMoodyAlgebra(QQ, cm)
            sage: KM.cartan_type()
            [ 2 -5]
            [-4  2]
        """
        return self._cartan_matrix.cartan_type()

    def cartan_matrix(self):
        """
        Return the defining Cartan matrix of ``self``.

        EXAMPLES::

            sage: KM = KacMoodyAlgebra(QQ, ['B',2])
            sage: KM.cartan_matrix()
            [ 2 -1]
            [-2  2]
            sage: cm = CartanMatrix([[2,-5],[-4,2]])
            sage: KM = KacMoodyAlgebra(QQ, cm)
            sage: KM.cartan_matrix()
            [ 2 -5]
            [-4  2]
        """
        return self._cartan_matrix

    @cached_method
    def basis(self):
        """
        Return a basis of ``self``.

        EXAMPLES::

            sage: KM = KacMoodyAlgebra(QQ, ['B',2])
            sage: sorted(KM.basis(), cmp=lambda x,y: cmp(x.value, y.value))
            [e1,
             -[e1, [e1, e2]],
             [e1, [e1, e2]],
             -[e1, e2],
             [e1, e2],
             e2,
             f1,
             -[f1, [f1, f2]],
             [f1, [f1, f2]],
             -[f1, f2],
             [f1, f2],
             f2,
             h1,
             h2]
        """
        if not self._cartan_matrix.is_finite():
            raise NotImplementedError

        todo_e = list(self._e)
        todo_f = list(self._f)
        basis = set(self.gens())
        while len(todo_e) > 0:
            next_e = todo_e.pop()
            next_f = todo_f.pop()
            for i,e in enumerate(self._e):
                f = self._f[i]
                print e, next_e
                b = self.bracket(e, next_e)
                if b != self.zero() and b not in basis:
                    print "===== new basis element ===="
                    print b
                    print "----------------------------"
                    basis.add(b)
                    todo_e.append(b)
                    b = self.bracket(f, next_f)
                    basis.add(b)
                    todo_f.append(b)

        return Family(basis)

    def alpha(self, i, x):
        r"""
        The `i`-th simple root `\alpha_i` applied to ``x``.
        """
        for j in range(self._n):
            if x == self._h[j]:
                return self._cartan_matrix[i][j]
        return 0

    def e(self, i):
        r"""
        Return the generator `e_i`.
        """
        return self._e[i]

    def f(self, i):
        r"""
        Return the generator `f_i`.
        """
        return self._f[i]

    def h(self, i):
        """
        Return the generator `h_i`.
        """
        return self._h[i]

    def d(self, i):
        """
        Return the generator `d_i`.
        """
        return self._d[i]

# This will need work
class NoSerreKacMoodyAlgebra(KacMoodyAlgebra):
    """
    This is the Lie algebra generated by `e_i`, `f_i`, and `h_i` subject
    with all relations as for Kac-Moody algebras except for the Serre
    relations.
    """
    @staticmethod
    def __classcall_private__(cls, R, cm, index_set=None):
        """
        Normalize arguments to ensure a unique representation.

        TESTS::

            sage: KM1 = NoSerreKacMoodyAlgebra(QQ, ['A', 3])
            sage: KM2 = NoSerreKacMoodyAlgebra(QQ, 'A3')
            sage: KM3 = NoSerreKacMoodyAlgebra(QQ, CartanType(['A', 3]))
            sage: KM1 is KM2 and KM2 is KM3
            True
        """
        cm = CartanMatrix(cm, index_set=index_set) # index_set not currently supported
        return super(KacMoodyAlgebra, cls).__classcall__(cls, R, cm)

    def __init__(self, R, cm):
        """
        Initialize ``self``.

        TESTS::
        """
        index_set = cm.index_set()
        n = len(index_set)
        self._n = n
        self._cartan_matrix = cm

        # Construct the base free Lie algebra
        names  = ['e%s'%i for i in index_set]
        names += ['f%s'%i for i in index_set]
        names += ['h%s'%i for i in index_set]
        names += ['d%s'%i for i in range(n-cm.rank())]
        F = FreeLieAlgebra(R, names).Lyndon()
        e = F.gens()[:n]
        f = F.gens()[n:2*n]
        h = F.gens()[2*n:3*n]
        d = F.gens()[3*n:]

        # Construct the basic relations
        P_check = h+d
        rels = [F.bracket(x, y) for i,x in enumerate(P_check) for y in P_check[i+1:]]
        rels += [F.bracket(x, y) - cm[i,j]*x for i,x in enumerate(e) for j,y in enumerate(h)]
        rels += [F.bracket(x, y) for x in e for y in d]
        rels += [F.bracket(x, y) + cm[i,j]*x for i,x in enumerate(f) for j,y in enumerate(h)]
        rels += [F.bracket(x, y) for x in f for y in d]
        for i,x in enumerate(e):
            for j,y in enumerate(f):
                if i == j:
                    rels.append(F.bracket(x, y) - h[i])
                else:
                    rels.append(F.bracket(x, y))

        category = LieAlgebras(R).WithBasis()
        if cm.is_finite():
            category = category.FiniteDimensional()
        FinitelyPresentedLieAlgebra.__init__(self, F, rels, category=category)

        # Define the actual generators
        elt = lambda x: self.element_class(self, x)
        self._e = map(elt, e)
        self._f = map(elt, f)
        self._h = map(elt, h)
        self._d = map(elt, d)

class LieAlgebraChevalleyBasis(FinitelyGeneratedLieAlgebra, IndexedGenerators):
    r"""
    A simple Lie algebra in the Chevalley basis.

    Let `L` be a simple complex Lie algebra with roots `\Phi`, then the
    Chevalley basis is given by `e_{\alpha}` for all `\alpha \in \Phi` and
    `h_{\alpha_i} := h_i` where `\alpha_i` is a simple root subject. These
    generators are subject to the relations:

    .. MATH::

        \begin{aligned}
        [h_i, h_j] & = 0
        \\ [h_i, e_{\beta}] & = A_{\alpha_i, \beta} e_{\beta}
        \\ [e_{\beta}, e_{-\beta}] & = \sum_i A_{\beta, \alpha_i} h_i
        \\ [e_{\beta}, e_{\gamma}] & = \begin{cases}
        N_{\beta,\gamma} e_{\beta + \gamma} & \beta + \gamma \in \Phi \\
        0 & \text{otherwise.} \end{cases}
        \end{aligned}

    where `A_{\alpha, \beta} = \frac{2 (\alpha, \beta)}{(\alpha, \alpha)}` and
    `N_{\alpha, \beta}` is the maximum such that
    `\alpha - N_{\alpha, \beta} \beta \in \Phi`.

    For computing the signs of the coefficients, see Section 3 of [CMT]_.

    REFERNCES:

    .. [CMT] A. M. Cohen, S. H. Murray, D. E. Talyor. *Groups of Lie type*.
       http://www.win.tue.nl/~amc/pub/papers/cmt.pdf
    """
    @staticmethod
    def __classcall_private__(cls, R, cartan_type):
        """
        Normalize ``self`` to ensure a unique represntation.

        TESTS::

            sage: L1 = LieAlgebra(QQ, cartan_type=['A', 2])
            sage: L2 = LieAlgebra(QQ, CartanType(['A', 2]))
            sage: L3 = LieAlgebra(QQ, cartan_type=CartanMatrix(['A', 2]))
            sage: L1 is L2 and L2 is L3
            True
        """
        return super(LieAlgebraChevalleyBasis, cls).__classcall__(
            cls, R, CartanType(cartan_type))

    def __init__(self, R, cartan_type):
        r"""
        Initialize ``self``.

        TESTS::

            sage: TestSuite(LieAlgebra(QQ, cartan_type=['A',2])).run()
        """
        self._cartan_type = cartan_type
        RL = cartan_type.root_system().root_lattice()
        alpha = RL.simple_roots()
        p_roots = list(RL.positive_roots_by_height())
        n_roots = map(lambda x: -x, p_roots)
        alphacheck = RL.simple_coroots()
        roots = RL.roots()
        num_sroots = len(alpha)

        names = ['h%s'%i for i in range(1, num_sroots+1)]
        e_names = ['e%s'%i for i in range(1, num_sroots+1)]
        f_names = ['f%s'%i for i in range(1, num_sroots+1)]

        # Determine the signs for the structure coefficients from the root system
        # We first create the special roots
        sp_sign = {}
        for i,a in enumerate(p_roots):
            for b in p_roots[i+1:]:
                if a + b not in p_roots:
                    continue

                # Compute the sign for the extra special pair
                x, y = (a + b).extraspecial_pair()

                if (x, y) == (a, b): # If it already is an extra special pair
                    sp_sign[(x, y)] = 1
                    sp_sign[(y, x)] = -1
                    continue

                if b - x in roots:
                    t1 = (b-x).norm_squared() / b.norm_squared() * sp_sign[(x, b-x)] * sp_sign[(a, y-a)]
                else:
                    t1 = 0
                if a - x in roots:
                    t2 = (a-x).norm_squared() / a.norm_squared() * sp_sign[(x, a-x)] * sp_sign[(b, y-b)]
                else:
                    t2 = 0

                if t1 - t2 > 0:
                    sp_sign[(a,b)] = 1
                else:
                    sp_sign[(a,b)] = -1
                sp_sign[(b,a)] = -sp_sign[(a,b)]

        # Function to construct the structure coefficients (up to sign)
        def e_coeff(r, s):
            p = 1
            while r - p*s in roots:
                p += 1
            return p

        # Now we can compute all necessary structure coefficients
        coeffs = {}
        for i,r in enumerate(p_roots):
            # [e_r, h_i] and [h_i, f_r]
            for ac in alphacheck:
                c = r.scalar(ac)
                if c == 0:
                    continue
                coeffs[(r, ac)] = {r: -c}
                coeffs[(ac, -r)] = {-r: -c}

            # [e_r, f_r]
            h_sum = {}
            for j, c in r.associated_coroot():
                h_sum[alphacheck[j]] = c
            coeffs[(r, -r)] = h_sum

            # [e_r, e_s] and [e_r, f_s] with r != +/-s
            for j, s in enumerate(p_roots[i+1:]):
                j += i+1
                # Since h(s) > h(r), we will always have s - r > 0 (if it is a root)
                # [e_r, f_s]
                if s - r in p_roots:
                    c = e_coeff(r, -s)
                    a,b = s-r, r
                    if p_roots.index(a) + 1 > p_roots.index(b): # Note a != b
                        c = -c * sp_sign[(b, a)]
                    else:
                        c *= sp_sign[(a, b)]
                    coeffs[(-r, s)] = {a: -c}
                    coeffs[(r, -s)] = {a: c}

                # [e_r, e_s]
                a = r + s
                if a in p_roots:
                    # (r, s) is a special pair
                    c = e_coeff(r, s) * sp_sign[(r, s)]
                    index = p_roots.index(r+s) + 1
                    coeffs[(r, s)] = {a: c}
                    coeffs[(-r, -s)] = {-a: -c}

        # Lastly, make sure a < b for all (a, b) in the coefficients and flip if necessary
        for a,b in coeffs.keys():
            if self._basis_cmp(a, b) > 0:
                coeffs[(b,a)] = {k:-v for k,v in coeffs[(a, b)].items()}
                del coeffs[(a, b)]
        self._s_coeff = coeffs

        names = e_names + f_names + names
        category = LieAlgebras(R).FiniteDimensional().WithBasis()
        FinitelyGeneratedLieAlgebra.__init__(self, R, names, category=category)
        IndexedGenerators.__init__(self, p_roots + n_roots + list(alphacheck),
                                   prefix='E', monomial_cmp=self._basis_cmp)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: LieAlgebra(QQ, cartan_type=['A', 2])
            Lie algebra of ['A', 2] in the Chevalley basis
        """
        return "Lie algebra of %s in the Chevalley basis"%self._cartan_type

    def _repr_generator(self, m):
        """
        Return a string representation of the basis element indexed by ``m``.
        """
        if m in self._cartan_type.root_system().root_lattice().simple_coroots():
            return "h{}".format(m.support()[0])
        return IndexedGenerators._repr_generator(self, m)

    def _basis_cmp(self, x, y):
        """
        Compare two basis element indices. We order the basis elements by
        positive roots, coroots, and negative roots and then according to
        height.

        OUTPUT:

        If ``x == y``, return 0. If ``x < y``, return -1. Else return 1.
        """
        if x == y:
            return 0

        RL = self._cartan_type.root_system().root_lattice()
        p_roots = list(RL.positive_roots_by_height())
        n_roots = map(lambda x: -x, p_roots)
        alphacheck = RL.simple_coroots()

        if x in p_roots:
            if y in p_roots:
                return cmp(p_roots.index(x), p_roots.index(y))
            return -1

        if x in alphacheck:
            if y in p_roots:
                return 1
            if y in alphacheck:
                return cmp(x, y)
            return -1

        # x is in n_roots
        if y not in n_roots:
            return 1
        return cmp(n_roots.index(x), n_roots.index(y))

    def dimension(self):
        """
        Return the dimension of ``self``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, cartan_type=['A', 2])
            sage: L.dimension()
            8
        """
        return len(self.basis())

    def bracket_on_basis(self, x, y):
        """
        Return the Lie bracket of ``[x, y]``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, cartan_type=['A', 2])
            sage: e1,e2,f1,f2,h1,h2 = L.gens()
            sage: L.bracket(e1, f1) # indirect doctest
            h1
            sage: L.bracket(h1, e1)
            2*E[alpha[1]]
            sage: L.bracket(h1, f1)
            -2*E[-alpha[1]]
            sage: L.bracket(e1, e1)
            0
            sage: k = L.bracket(e1, e2); k
            E[alpha[1] + alpha[2]]
            sage: L.bracket(e1, k)
            0
        """
        b = (x, y)
        if b not in self._s_coeff:
            return self.zero()
        return self.element_class(self, self._s_coeff[b])

    def cartan_type(self):
        """
        Return the Cartan type of ``self``.
        """
        return self._cartan_type

    def weyl_group(self):
        """
        Return the Weyl group of ``self``.
        """
        from sage.combinat.root_system.weyl_group import WeylGroup
        return WeylGroup(self._cartan_type)

    def affine(self, kac_moody=False):
        """
        Return the (untwisted) affine Lie algebra of ``self``.
        """
        from sage.algebras.lie_algebras.classical_lie_algebra import AffineLieAlgebra
        return AffineLieAlgebra(self, kac_moody)

    # Useful in creating the UEA
    @cached_method
    def indices_to_positive_roots_map(self):
        """
        Return the map from indices to positive roots.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, cartan_type=['A', 2])
            sage: L.indices_to_positive_roots_map()
            {1: alpha[1], 2: alpha[2], 3: alpha[1] + alpha[2]}
        """
        RL = self._cartan_type.root_system().root_lattice()
        return {i+1:r for i,r in enumerate(RL.positive_roots())}

    @cached_method
    def gens(self):
        """
        Return the generators of ``self`` in the order of `e_i`, `f_i`,
        and `h_i`.
        """
        index_set = self._cartan_type.index_set()
        RL = self._cartan_type.root_system().root_lattice()
        alpha = RL.simple_roots()
        alphacheck = RL.simple_coroots()
        B = self.basis()

        ret = []
        for i in index_set:
            ret.append(B[alpha[i]])
        for i in index_set:
            ret.append(B[-alpha[i]])
        for i in index_set:
            ret.append(B[alphacheck[i]])
        return tuple(ret)

    def gen(self, i):
        """
        Return the ``i``-th generator of ``self``.
        """
        return self.gens()[i]

    def basis(self):
        """
        Return the basis of ``self``.
        """
        one = self.base_ring().one()
        return Family(self._indices, lambda x: self.element_class(self, {x: one}))

    def algebra_generators(self):
        """
        Return the Lie algebra generators of ``self``.
        """
        RL = self._cartan_type.root_system().root_lattice()
        alpha = RL.simple_roots()
        alphacheck = RL.simple_coroots()
        d = {}
        B = self.basis()
        for i in self._cartan_type.index_set():
            d['e{}'.format(i)] = B[alpha[i]]
            d['f{}'.format(i)] = B[-alpha[i]]
            d['h{}'.format(i)] = B[alphacheck[i]]
        return Family(d)

    def highest_root_basis_elt(self, pos=True):
        r"""
        Return the basis element corresponding to the highest root `\theta`.
        If ``pos`` is ``True``, then returns `e_{\theta}`, otherwise it
        returns `f_{\theta}`.
        """
        RL = self._cartan_type.root_system().root_lattice()
        theta = RL.highest_root()
        B = self.basis()
        if pos:
            return B[theta]
        return B[-theta]

    class Element(LieAlgebraElement):
        def _repr_(self):
            r"""
            Return a string representation of ``self``.
            """
            return repr_lincomb(self._sorted_items_for_printing(),
                                scalar_mult=self.parent()._print_options['scalar_mult'],
                                repr_monomial = self.parent()._repr_generator,
                                strip_one = True)

        def _latex_(self):
            r"""
            Return a `\LaTeX` representation of ``self``.
            """
            return repr_lincomb(self._sorted_items_for_printing(),
                                scalar_mult       = self.parent()._print_options['scalar_mult'],
                                latex_scalar_mult = self.parent()._print_options['latex_scalar_mult'],
                                repr_monomial = self.parent()._latex_generator,
                                is_latex=True, strip_one = True)

