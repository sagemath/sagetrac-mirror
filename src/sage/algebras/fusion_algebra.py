r"""
Fusion algebras

AUTHORS:

- Travis Scrimshaw (2013-10-18): Initial version
"""

#*****************************************************************************
#  Copyright (C) 2013 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.misc_c import prod

from sage.categories.finite_dimensional_algebras_with_basis import FiniteDimensionalAlgebrasWithBasis
from sage.categories.commutative_algebras import CommutativeAlgebras
from sage.rings.all import ZZ
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.backtrack import TransitiveIdeal
from sage.combinat.crystals.monomial_crystals import CrystalOfNakajimaMonomials

class FusionAlgebra(CombinatorialFreeModule):
    r"""
    The fusion algebra with a given affine Cartan type.

    Abstractly a *fusion algebra* `F` is a finite dimensional commutative
    associative algebra over a ring `R` (typically `\QQ`) with a basis
    `(B_a)_{a \in I}` indexed by `I` with structure coefficients
    `N_{a,b}^c \in \ZZ_{\geq 0}` given by

    .. MATH::

        B_a \cdot B_b = \sum_{c \in I} N_{a,b}^c B_c.

    Let `D` denote the dimension of `F`. There must exist a distinguished
    index `\Omega \in I` such that the matrix
    `C_{\Omega} = (N_{a,b}^{\Omega})_{a,b \in I}` satisfies
    `C_{\Omega}^2 = I_D`. We call the matrix `C_{\Omega}` a
    *conjugation matrix*. Thus we can define an order 2 permutation
    `\sigma_{\Omega}` whose matrix is `C_{\Omega}` and
    `N_{a,b}^{\sigma_{\Omega}(c)}` is symmetric for all `a,b,c \in I`.
    The multiplication of two basis elements in a fusion algebra is known as
    a *fusion rule*.

    Using a conjugation matrix `C_{\Omega}` we can define an automorphism
    `\mathcal{C}_{\Omega}` called *conjugation* given by

    .. MATH::

        \mathcal{C}_{\Omega}(B_a) = \sum_{b \in I} C^{\Omega}_{a,b} B_b.

    We note that `\mathcal{C}_{\Omega}^2 = I_F`.

    Fusion algebras are strongly tied to WZW and conformal field theories
    and encode information about the operator product of two primary fields
    of the conformal field theory. In particular, the coefficients given from
    the operator product vanish if and only if the corresponding fusion
    product coefficients vanish. See [Fuchs1992]_ for more information.

    Let `\mathfrak{g}` denote an affine Kac-Moody algebra and the
    corresponding classical weights by `\overline{\beta}`. Let
    `\overline{\rho}` denote the (classical) Weyl vector and `\theta` denote
    the (classical) highest root of `\overline{\mathfrak{g}}`. Let `W` be the
    (affine) Weyl group of `\mathfrak{g}` and `h^{\vee}` be the dual coxeter
    number of `\overline{\mathfrak{g}}`. We can define fusion algebras using
    level `k` representations of `\mathfrak{g}` which has a basis consisting
    of level `k` dominant weights and the structure coefficients are given
    by the Kac-Waldron formula:

    .. MATH::

        N_{\lambda,\mu}^{(k) \, \nu} = \sum_{w \in W} \epsilon(w)
        M_{\overline{\lambda}}\bigl( w(\overline{\nu} + \overline{\rho})
        - \overline{\mu} - \overline{\rho} \bigr)

    where `M_{\overline{\lambda}}(\beta)` is the multiplicity of the
    weight `\beta` in the (classical) highest weight representation
    `\overline{\lambda}`, and we have

    .. MATH::

        s_0(\beta) = s_{\theta}(\beta) + (k + h^{\vee}) \theta.

    We note we only need classical weights to define the fusion algebra since
    we can describe level `k` weights only using classical weights `\beta`
    such that `\langle \beta, \theta \rangle \leq k` and
    `k - \langle \beta, \theta \rangle \equiv 0 \mod l_0` where `l_0` is the
    level of `\Lambda_0`.

    INPUT:

    - ``base_ring`` -- the base ring
    - ``cartan_type`` -- the affine Cartan type
    - ``level`` -- the level

    .. TODO::

        Implement for twisted types following [QRS2002]_.

    EXAMPLES::

        sage: FusionAlgebra(QQ, ['A',2,1], 4)
        The fusion algebra ['A', 2, 1] of level 4 over Rational Field

    REFERENCES:

    .. [Feingold04] Alex J. Feingold.
       *Fusion rules for affine Kac-Moody algebras*. Contemporary Mathematics.
       Volume **343** (2004).

    .. [Fuchs1992] Jurgen Fuchs. *Affine Lie Algebras and Quantum Groups*.
       Cambridge Monographs on Mathematical Physics.
       Cambridge University Press. (1992).

    .. [QRS2002] Thomas Quella, Ingo Runkel, and Christoph Schweigert.
       *An algorithm for twisted fusion rules*. Adv. Theor. Math. Phys.
       **6** (2002). 197-205. :arxiv:`0203133`.
    """
    @staticmethod
    def __classcall_private__(cls, base_ring, cartan_type, level):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: F1 = FusionAlgebra(QQ, ['A',2,1], 4)
            sage: F2 = FusionAlgebra(QQ, ('A',2,1), int(4))
            sage: F3 = FusionAlgebra(QQ, CartanType('A2~'), 4)
        """
        cartan_type = CartanType(cartan_type)
        if not cartan_type.is_affine():
            raise ValueError("the Cartan type must be affine")
        if not cartan_type.is_untwisted_affine():
            raise NotImplementedError("only implemented for untwisted types")
        return super(FusionAlgebra, cls).__classcall__(cls, base_ring, cartan_type, level)

    def __init__(self, base_ring, cartan_type, level):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: F = FusionAlgebra(QQ, ['A',2,1], 4)
            sage: TestSuite(F).run()
            sage: F = FusionAlgebra(QQ, ['C',2,1], 3)
            sage: TestSuite(F).run()
        """
        self._cartan_type = cartan_type
        self._level = level
        category = (FiniteDimensionalAlgebrasWithBasis(base_ring), CommutativeAlgebras(base_ring))

        P = cartan_type.root_system().weight_lattice()
        P_cl = cartan_type.classical().root_system().weight_lattice()

        # Generate all dominant classical weights (thought of as affine weights)
        #   with level <= self._level
        levels = {i:la.level() for i,la in dict(P.fundamental_weights()).items()}
        La = P_cl.fundamental_weights()
        next_weights = lambda x: [(levels[i] + x[0], x[1] + la)
                                  for i,la in dict(La).items()
                                  if levels[i] + x[0] <= level]
        basis_keys = TransitiveIdeal(next_weights, [(0, P_cl.zero())])
        basis_keys = [wt for l,wt in basis_keys if (level - l) % levels[0] == 0]

        CombinatorialFreeModule.__init__(self, base_ring, tuple(basis_keys),
                                         prefix='F', category=category)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: FusionAlgebra(QQ, ['A',2,1], 4)
            The fusion algebra ['A', 2, 1] of level 4 over Rational Field
        """
        return "The fusion algebra {} of level {} over {}".format(
                self._cartan_type, self._level, self.base_ring())

    def level(self):
        """
        Return the level of ``self``.

        EXAMPLES::

            sage: F = FusionAlgebra(QQ, ['A',2,1], 4)
            sage: F.level()
            4
        """
        return self._level

    def cartan_type(self):
        """
        Return the Cartan type of ``self``.

        EXAMPLES::

            sage: F = FusionAlgebra(QQ, ['A',2,1], 4)
            sage: F.cartan_type()
            ['A', 2, 1]
        """
        return self._cartan_type

    @cached_method
    def one_basis(self):
        """
        Return the basis index of the element `1`.

        EXAMPLES::

            sage: F = FusionAlgebra(QQ, ['A',2,1], 4)
            sage: F.one_basis()
            0
            sage: F.one_basis().parent()
            Weight lattice of the Root system of type ['A', 2]
        """
        return self._cartan_type.classical().root_system().weight_lattice().zero()

    @cached_method
    def product_on_basis(self, la, mu):
        """
        Return the fusion product of two basis elements indexed by
        ``la`` and ``mu``.

        EXAMPLES::

            sage: F = FusionAlgebra(QQ, ['A',1,1], 3)
            sage: la = RootSystem(['A',1]).weight_lattice().fundamental_weight(1)
            sage: B = [F.basis()[i*la] for i in range(4)]; B
            [F[0], F[Lambda[1]], F[2*Lambda[1]], F[3*Lambda[1]]]
            sage: F0,F1,F2,F3 = B
            sage: F2 * F3 # indirect doctest
            F[Lambda[1]]
            sage: F1 * F3
            F[2*Lambda[1]]
            sage: F1 * F2
            F[Lambda[1]] + F[3*Lambda[1]]

            sage: F = FusionAlgebra(QQ, ['A',2,1], 2)
            sage: La = RootSystem(['A',2]).weight_lattice().fundamental_weights()
            sage: B = F.basis()
            sage: B[La[1]+La[2]] * B[La[1]+La[2]]
            F[0] + F[Lambda[1] + Lambda[2]]
            sage: B[La[1]] * B[La[1]+La[2]]
            F[Lambda[1]] + F[2*Lambda[2]]
            sage: B[2*La[1]] * B[2*La[2]]
            F[0]

            sage: F = FusionAlgebra(QQ, ['B',3,1], 1)
            sage: matrix([[x*y for y in F.basis()] for x in F.basis()])
            [               F[0]        F[Lambda[3]]        F[Lambda[1]]]
            [       F[Lambda[3]] F[0] + F[Lambda[1]]        F[Lambda[3]]]
            [       F[Lambda[1]]        F[Lambda[3]]                F[0]]

            sage: F = FusionAlgebra(QQ, ['C',2,1], 1)
            sage: matrix([[x*y for y in F.basis()] for x in F.basis()])
            [               F[0]        F[Lambda[2]]        F[Lambda[1]]]
            [       F[Lambda[2]]                F[0]        F[Lambda[1]]]
            [       F[Lambda[1]]        F[Lambda[1]] F[0] + F[Lambda[2]]]

        TESTS:

        We check that the product is commutative::

            sage: F = FusionAlgebra(QQ, ['A',2,1], 3)
            sage: all(x*y == y*x for x in F.basis() for y in F.basis())
            True
            sage: F = FusionAlgebra(QQ, ['C',2,1], 2)
            sage: all(x*y == y*x for x in F.basis() for y in F.basis())
            True
            sage: F = FusionAlgebra(QQ, ['D',4,1], 1)
            sage: all(x*y == y*x for x in F.basis() for y in F.basis())
            True
        """
        # We'll need to do something slightly different for twisted types
        classical = self._cartan_type.classical()

        # Compute the set of shifted weights of B(la)
        P = classical.root_system().weight_lattice()
        # We use Nakajima monomials because they are currently the fastest
        #   to iterate over
        B = CrystalOfNakajimaMonomials(classical, la)
        rho = P.rho()
        index_set = classical.index_set()
        alpha = P.simple_roots()
        def next_elt(x):
            ret = []
            wt, elt = x
            for i in index_set:
                next = elt.f(i)
                if next is not None:
                    ret.append((wt - alpha[i], next))
            return ret
        weights = TransitiveIdeal(next_elt, [(rho + mu + la, B.highest_weight_vector())])

        # Convert that into multiplicities
        mults = {}
        for wt,elt in weights:
            mults[wt] = mults.get(wt, 0) + 1

        ret = []
        Q = classical.root_system().root_lattice()
        theta = P.highest_root()
        thetacheck = Q.highest_root().associated_coroot()
        k = self._level
        n = classical.rank()
        wall = k + classical.dual_coxeter_number()
        for wt,mult in mults.items():
            sgn = 1
            wt, red = wt.to_dominant_chamber(reduced_word=True)
            sgn *= (-1)**len(red)
            while wt.scalar(thetacheck) > wall:
                wt = wt.reflection(thetacheck, use_coroot=True) + wall*theta
                sgn = -sgn
                wt, red = wt.to_dominant_chamber(reduced_word=True)
                sgn *= (-1)**len(red)
            if len(wt) == n and wt.scalar(thetacheck) != wall:
                ret.append((wt - rho, sgn * mult))

        return self.sum_of_terms(ret)

    @cached_method
    def conjugation_matrices(self):
        r"""
        Return the conjugation matrices `C_{\Omega}` of ``self`` indexed by
        their distinguished indices `\Omega`.

        EXAMPLES::

            sage: F = FusionAlgebra(QQ, ['D',4,1], 1)
            sage: M = F.conjugation_matrices()
            sage: sorted(M)
            [
            [0 0 0 1]  [0 0 1 0]  [0 1 0 0]  [1 0 0 0]
            [0 0 1 0]  [0 0 0 1]  [1 0 0 0]  [0 1 0 0]
            [0 1 0 0]  [1 0 0 0]  [0 0 0 1]  [0 0 1 0]
            [1 0 0 0], [0 1 0 0], [0 0 1 0], [0 0 0 1]
            ]
        """
        from sage.matrix.constructor import matrix
        from sage.matrix.constructor import identity_matrix
        from sage.sets.family import Family
        K = self.get_order()
        B = self.basis()
        B = [B[k] for k in K] # Order accordingly
        n = len(B)
        R = self.base_ring()
        I = identity_matrix(R, n)
        array = [[x*y for y in B] for x in B]
        ret = {}
        for i in K:
            C = matrix(R, [[val.coefficient(i) for val in row] for row in array])
            if C*C == I:
                ret[i] = C
        return Family(ret)

    def distinguished_indices(self):
        r"""
        Return the distinguished indices `\Omega` such that the associated
        conjugation matrix `C_{\Omega} = (N_{a,b}^{\Omega})_{a,b in I}`.

        EXAMPLES::

            sage: F = FusionAlgebra(QQ, ['D',4,1], 1)
            sage: sorted(F.distinguished_indices())
            [0, Lambda[1], Lambda[3], Lambda[4]]
        """
        return self.conjugation_matrices().keys()

    def conjugation(self, omega):
        r"""
        Return the conjugation morphism `\mathcal{C}_{\Omega}` of ``self``.

        INPUT:

        - ``omega`` -- the distinguished index `\Omega`

        EXAMPLES::

            sage: F = FusionAlgebra(QQ, ['D',4,1], 1)
            sage: omega = F.distinguished_indices()[0]; omega
            Lambda[4]
            sage: conj = F.conjugation(omega); conj
            Generic endomorphism of The fusion algebra ['D', 4, 1] of level 1 over Rational Field
            sage: F.an_element()
            2*F[0] + 3*F[Lambda[3]] + 2*F[Lambda[4]]
            sage: conj(F.an_element())
            2*F[0] + 3*F[Lambda[1]] + 2*F[Lambda[4]]
            sage: conj(conj(F.an_element()))
            2*F[0] + 3*F[Lambda[3]] + 2*F[Lambda[4]]
        """
        C = self.conjugation_matrices()[omega]
        K = self.get_order()
        conj = lambda x: self.sum_of_terms([(k, C[K.index(x)][i]) for i,k in enumerate(K)], distinct=True)
        return self.module_morphism(conj, codomain=self)

