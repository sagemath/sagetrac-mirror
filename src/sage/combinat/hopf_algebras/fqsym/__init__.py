# -*- coding: utf-8 -*-
r"""
The combinatorial Hopf algebra of Free Quasi-Symmetric functions

This module implements methods related to the Hopf algebra of permutations
called Free Quasi-Symmetric functions Hopf algebra or the Malvenuto-Reutenauer
Hopf algebra (see [MalReut]_, [NCSF-VI]_, [NCSF-VII]_ and [AguSot]_).

References
----------

.. [MalReut] Duality between quasi-symmetrical functions and the solomon
    descent algebra,
    Claudia Malvenuto and
    Christophe Reutenauer

.. [NCSF-VI] Noncommutative Symmetric Function VI: Free Quasi-Symmetric
    Functions and Related Algebras,
    Gérard Duchamp,
    Florent Hivert and
    Jean-Yves Thibon

.. [NCSF-VII] Noncommutative symmetric functions VII : free quasi-symmetric
    functions revisited,
    Gérard Duchamp,
    Florent Hivert,
    Jean-Christophe Novelli and
    Jean-Yves Thibon

.. [AguSot] Structure of the Malvenuto-Reutenauer Hopf algebra of permutations,
    Marcelo Aguiar and
    Frank Sottile

.. [AvaVien] The product of trees in the Loday-Ronco algebra through Catalan
    alternative tableaux,
    Jean-Christophe Aval and
    Xavier Viennot,

.. [AvNoThi] The # product in combinatorial Hopf algebras,
    Jean-Christophe Aval,
    Jean-Christophe Novelli and
    Jean-Yves Thibon


Description
-----------

The *Free Quasi-Symmetric function* or *Malvenuto-Reutenauer Hopf algebra* is
the permutations Hopf algebra.

AUTHOR:
-------

    - Jean-Baptiste Priez
"""
#*****************************************************************************
#       Copyright (C) 2013 Jean-Baptiste Priez <jbp@kerios.fr>.
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.categories.category import Category
from sage.categories.hopf_algebras_with_basis import HopfAlgebrasWithBasis
from sage.categories.realizations import Category_realization_of_parent
from sage.combinat.hopf_algebras import GenericGradedConnexeHopfAlgebras, GenericBasisOfGCHopfAlgebra, \
                                        words_like_getitem
from sage.combinat.hopf_algebras.categories.diese_product import DieseProductAlgebras
from sage.combinat.hopf_algebras.categories.polynomial_realization import PolynomialRealizationAlgebras
from sage.combinat.hopf_algebras.categories.scalar_product import ScalarProductAlgebras
from sage.combinat.ncsf_qsym.generic_basis_code import GradedModulesWithInternalProduct
from sage.combinat.permutation import Permutations
from sage.categories.hopf_algebras import HopfAlgebras
from sage.categories.bidendriform_bialgebras import BidendriformBialgebras


class FreeQuasiSymmetricFunctions(GenericGradedConnexeHopfAlgebras):
    r"""
    The combinatorial Hopf algebra of permutations is called the Free Quasi-
    Symmetric functions Hopf algebra or the Malvenuto-Reutenauer Hopf algebra
    (see [MalReut]_, [NCSF-VI]_, [NCSF-VII]_ and [AguSot]_).


    We construct this abstract Hopf algebra of free quasi-symmetric functions
    over the rational numbers::

        sage: FQS = FreeQuasiSymmetricFunctions(QQ); FQS # FQSym(QQ)
        The combinatorial Hopf algebra of Free Quasi-Symmetric Functions over the Rational Field

    We purpose several realizations of FQSym::

        sage: F = FQS.Fundamental()       # FQS.F()
        sage: G = FQS.FundamentalDual()   # FQS.G()
        sage: E = FQS.Elementary()        # FQS.E()
        sage: H = FQS.Homogene()          # FQS.H()
        sage: n = FQS.ElementaryDual()    # FQS.n()
        sage: m = FQS.HomogeneDual()      # FQS.m()

    And on these bases, we could use several operators:

        - product, coproduct::

            sage: F[3,1,2] * F[1]
            F[3, 1, 2, 4] + F[3, 1, 4, 2] + F[3, 4, 1, 2] + F[4, 3, 1, 2]
            sage: F[1,2,3].coproduct()
            F[] # F[1, 2, 3] + F[1] # F[1, 2] + F[1, 2] # F[1] + F[1, 2, 3] # F[]
            sage: E[3,1,2] * E[1,3,4,2]
            E[3, 1, 2, 4, 6, 7, 5]

        - antipode::

            sage: G[3,1,2].antipode()
            G[1, 3, 2] - G[2, 1, 3] + G[2, 3, 1] - 2*G[3, 1, 2]

        - internal product (composition of permutations)::

            sage: F[3,1,2].internal_product(F[3,1,2])
            F[2, 3, 1]

        - scalar product::

            sage: G[4,1,2,5,3].scalar_product(F[1,2,3]*F[1,2])
            1

        - #-product [AvaVien]_ and [AvNoThi]_::

            sage: G[3,1,2].diese_product(G[3,1,2])
            G[5, 1, 4, 2, 3] + G[5, 2, 4, 1, 3] + G[5, 3, 4, 1, 2]

        - ((bi)-dendriform operations `\prec`, `\succ`, `\Delta_{\prec}` and
    `\Delta_{\succ}`.) (C3 dependants)::

            sage: F[1,2]>>F[2,1]
            F[1, 2, 4, 3] + F[1, 4, 2, 3] + F[4, 1, 2, 3]
            sage: F[3,1,4,2].left_coproduct()
            F[2, 1, 3] # F[1]

        - expand the polynomial realization::

            sage: G[1,2].expand_to_polynomial(3)
            a1^2 + a1*a2 + a1*a3 + a2^2 + a2*a3 + a3^2

    FQSym is the free graded connected Hopf algebra, with
    generators are irreductible permutations.

    We use the Sage standard renaming idiom to get shorter outputs::

        sage: FQS.rename("FQSym")
        sage: FQS
        FQSym

    FQSym has many representations as a concrete Hopf algebra. Each of them
    has a distinguished basis, and its elements are expanded in this
    basis. Here is the fundamental representation::

        sage: F = FQS.F(); F
        FQSym on the Fundamental basis

    Elements of F are linear combinations of such permutations::

        sage: F.an_element()
        2*F[] + 2*F[1] + 3*F[1, 2]

    To construct an element one can therefore do::

        sage: F(Permutation([3,1,2]))
        F[3, 1, 2]

    As this is rather cumbersome, the following abuses of notation are
    allowed::

        sage: F([3,1,2])
        F[3, 1, 2]
        sage: F[[3,1,2]]
        F[3, 1, 2]
        sage: F[3,1,2]
        F[3, 1, 2]

    Unfortunately, due to a limitation in Python syntax, one cannot use::

        sage: F[]       # not implemented

    Instead, you can use::

        sage: F[[]]
        F[]

    Now, we can construct linear combinations of basis elements::

        sage: F[2,1,3] + 2 * (F[1,3,4,2] + F[2,1])
        2*F[1, 3, 4, 2] + 2*F[2, 1] + F[2, 1, 3]

    .. rubric:: Algebra structure

    To start with, F is a graded algebra, the grading being induced by the
    size of permutations. Due to this, the one is the basis element indexed by
    the empty permutation::

        sage: F.one()
        F[]
        sage: G.one()
        G[]
        sage: E.one()
        E[]

    As we have seen above, the ``F`` basis has a product::

        sage: F[3,1,2] * F[2,1]
        F[3, 1, 2, 5, 4] + F[3, 1, 5, 2, 4] + F[3, 1, 5, 4, 2] + F[3, 5, 1, 2, 4] + F[3, 5, 1, 4, 2] + F[3, 5, 4, 1, 2] + F[5, 3, 1, 2, 4] + F[5, 3, 1, 4, 2] + F[5, 3, 4, 1, 2] + F[5, 4, 3, 1, 2]

    .. rubric:: Hopf algebra structure

    ``F`` is further endowed with a coalgebra algebra structure. The
    coproduct is an algebra morphism, and therefore determined by its values on
    the generators; those are primitive::

        sage: F[3,1,2].coproduct()
        F[] # F[3, 1, 2] + F[1] # F[1, 2] + F[2, 1] # F[1] + F[3, 1, 2] # F[]
        sage: F[2,1].coproduct() * F[3,1,2].coproduct() == \
        ....: (F[2,1] * F[3,1,2]).coproduct()
        True

    The antipode is an anti-algebra morphism; in the ``F`` basis, it sends
    the generators to their opposite::

        sage: F[1,2,3].antipode()
        -F[3, 2, 1]
        sage: F[1,3,2].antipode()
        -F[2, 1, 3] - F[2, 3, 1] + F[3, 1, 2]

    The counit is defined by sending all elements of positive degree to zero::

        sage: (3 + F[3,1,2] + 2*F[1]).counit()
        3

    .. rubric:: Other concrete representations

    FQSym admits many other concrete realizations::

        sage: G = FQS.FundamentalDual()   # FQS.G()
        sage: E = FQS.Elementary()        # FQS.E()
        sage: H = FQS.Homogene()          # FQS.H()
        sage: n = FQS.ElementaryDual()    # FQS.n()
        sage: m = FQS.HomogeneDual()      # FQS.m()

    To change from one basis to another, one simply does::

        sage: G(F[3,1,2] + 2*G[4,1,3,2])
        G[2, 3, 1] + 2*G[4, 1, 3, 2]

    In general, one can mix up different basis in computations::

        sage: G[1,2]*E[1,2]
        G[1, 2, 3, 4] + G[1, 2, 4, 3] + G[1, 3, 2, 4] + G[1, 3, 4, 2] + G[1, 4, 2, 3] + G[1, 4, 3, 2] + G[2, 3, 1, 4] + G[2, 3, 4, 1] + G[2, 4, 1, 3] + G[2, 4, 3, 1] + G[3, 4, 1, 2] + G[3, 4, 2, 1]
        sage: FQS.reset_name()

    TESTS::

        sage: FQS = FQSym(QQ)
        sage: TestSuite(FQS).run()


    Here some examples of scalar product:

    EXAMPLES::

        sage: F = FQSym(QQ).F()
        sage: F[3,1,2].scalar_product(F[2,3,1])
        1
        sage: F[3,1,2].scalar_product(F[3,1,2])
        0
        sage: matr = lambda X,Y, n: matrix([[X(sigma).\
                            scalar_product(Y(mu))\
                for mu in Permutations(n)]\
                for sigma in Permutations(n)])
        sage: matr(F,F,3)
        [1 0 0 0 0 0]
        [0 1 0 0 0 0]
        [0 0 1 0 0 0]
        [0 0 0 0 1 0]
        [0 0 0 1 0 0]
        [0 0 0 0 0 1]
        sage: G = FQSym(QQ).G()
        sage: matr(F,G,3)
        [1 0 0 0 0 0]
        [0 1 0 0 0 0]
        [0 0 1 0 0 0]
        [0 0 0 1 0 0]
        [0 0 0 0 1 0]
        [0 0 0 0 0 1]
        sage: E = FQSym(QQ).E()
        sage: matr(F,E,3)
        [1 0 0 0 0 0]
        [1 1 0 0 0 0]
        [1 0 1 0 0 0]
        [1 1 0 0 1 0]
        [1 0 1 1 0 0]
        [1 1 1 1 1 1]
        sage: H = FQSym(QQ).H()
        sage: matr(E,H,3)
        [1 2 2 3 3 6]
        [0 1 0 1 1 3]
        [0 0 1 1 1 3]
        [0 0 0 0 1 2]
        [0 0 0 1 0 2]
        [0 0 0 0 0 1]
        sage: m = FQSym(QQ).m()
        sage: matr(m,H,3)
        [1 0 0 0 0 0]
        [0 1 0 0 0 0]
        [0 0 1 0 0 0]
        [0 0 0 0 1 0]
        [0 0 0 1 0 0]
        [0 0 0 0 0 1]
        sage: matr(m,E,3)
        [ 0  0  0  1  1  1]
        [ 0  0  0  0 -1  0]
        [ 0  0  0 -1  0  0]
        [ 0  0 -1 -1  0 -1]
        [ 0 -1  0  0 -1 -1]
        [ 1  1  1  1  1  1]
        sage: n = FQSym(QQ).n()
        sage: matr(n,E,3)
        [1 0 0 0 0 0]
        [0 1 0 0 0 0]
        [0 0 1 0 0 0]
        [0 0 0 0 1 0]
        [0 0 0 1 0 0]
        [0 0 0 0 0 1]
        sage: matr(n,H,3)
        [ 1  1  1  1  1  1]
        [-1  0 -1 -1  0  0]
        [-1 -1  0  0 -1  0]
        [ 0 -1  0  0  0  0]
        [ 0  0 -1  0  0  0]
        [ 1  1  1  0  0  0]
        sage: matr(n,G,3)
        [ 1  0  0  0  0  0]
        [-1  1  0  0  0  0]
        [-1  0  1  0  0  0]
        [ 0 -1  0  1  0  0]
        [ 0  0 -1  0  1  0]
        [ 1  0  0 -1 -1  1]
    """

    the_category = lambda self, R: Category.join((
         HopfAlgebras(R).Graded().Connected().WithRealizations(),
         BidendriformBialgebras(R).WithRealizations()
    ))

    def a_realization(self):
        return self.F()

    def dual(self):
        return self

    def _repr_(self):
        return "The combinatorial Hopf algebra of Free Quasi-Symmetric " + \
            "Functions over the %s" % self.base_ring()

    class Bases(Category_realization_of_parent):

        def super_categories(self):
            R = self.base().base_ring()
            return [self.base().Realizations(),
                HopfAlgebrasWithBasis(R).Graded().Connected().Realizations(),
                BidendriformBialgebras(R).WithBasis().Realizations(),
                DieseProductAlgebras(R).WithBasis().Realizations(),
                GradedModulesWithInternalProduct(R).WithBasis().Realizations(),
                PolynomialRealizationAlgebras(R).WithBasis().Realizations(),
                ScalarProductAlgebras(R).WithBasis().Realizations()]

        class Base(GenericBasisOfGCHopfAlgebra):
            _basis_indices = Permutations()
            _prefix = "** TO DEFINE **"

            __getitem__ = words_like_getitem