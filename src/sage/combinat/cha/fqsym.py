# -*- coding: utf-8 -*-
r"""
The combinatorial Hopf algebra of Free Quasi-Symmetric functions

This module implements method related to the Hopf algebra of permutations
called Free Quasi-Symmetric functions Hopf algebra or the Malvenuto-Reutenauer
Hopf algebra (see [MalReut]_, [NCSF-VI]_, [NCSF-VII]_ and [AguSot]_).

AUTHOR:

    - Jean-Baptiste Priez
    - Rémi Maurice

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
"""
#*****************************************************************************
#       Copyright (C) 2013 Jean-Baptiste Priez <jbp@kerios.fr>,
#                          Rémi Maurice <maurice@univ-mlv.fr>.
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.misc.lazy_import import LazyImport
from sage.categories.realizations import Category_realization_of_parent
from sage.combinat.permutation import Permutations
from sage.categories.graded_hopf_algebras_with_basis import \
    GradedHopfAlgebrasWithBasis
from sage.categories.graded_hopf_algebras import GradedHopfAlgebras
from sage.misc.abstract_method import abstract_method
from sage.combinat.ncsf_qsym.generic_basis_code import \
    GradedModulesWithInternalProduct
from sage.categories.category import Category

from categories.bidendriform_bialgebras import BidendriformBialgebras
from categories.diese_product_category import GradedAlgebrasWithDieseProduct
from tools.generic_basis import GenericBasis
import _fqsym


class FreeQuasiSymmetricFunctions(UniqueRepresentation, Parent):
    r'''
    The combinatorial Hopf algebra of permutations is called the Free Quasi-
    Symmetric functions Hopf algebra or the Malvenuto-Reutenauer Hopf algebra
    (see [MalReut]_, [NCSF-VI]_, [NCSF-VII]_ and [AguSot]_).

        The abstract algebra of non-commutative symmetric functions

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
            F[] # F[1, 2, 3] + F[1] # F[1, 2] + F[1, 2] # F[1] + F[1, 2, 3] # \
            F[]
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
            F[4, 1, 2, 3] + F[4, 1, 3, 2] + F[4, 3, 1, 2]
            sage: F.left_coproduct_on_basis([3,1,4,2])
            F[2, 1, 3] # F[1]

        - expand the polynomial realization::

            sage: G[1,2].expand(3)
            x_1^2 + x_1*x_2 + x_1*x_3 + x_2^2 + x_2*x_3 + x_3^2

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
        F[] + 2*F[1] + 3*F[1, 2]

    To construct an element one can therefore do::

        sage: F([Permutation([3,1,2])])
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
        F[3, 1, 2, 5, 4] + F[3, 1, 5, 2, 4] + F[3, 1, 5, 4, 2] + \
        F[3, 5, 1, 2, 4] + F[3, 5, 1, 4, 2] + F[3, 5, 4, 1, 2] + \
        F[5, 3, 1, 2, 4] + F[5, 3, 1, 4, 2] + F[5, 3, 4, 1, 2] + \
        F[5, 4, 3, 1, 2]

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
        G[1, 2, 3, 4] + G[1, 2, 4, 3] + G[1, 3, 2, 4] + G[1, 3, 4, 2] +\
        G[1, 4, 2, 3] + G[1, 4, 3, 2] + G[2, 3, 1, 4] + G[2, 3, 4, 1] +\
        G[2, 4, 1, 3] + G[2, 4, 3, 1] + G[3, 4, 1, 2] + G[3, 4, 2, 1]

    TESTS::

        sage: FQS = FQSym(QQ)
        sage: TestSuite(FQS).run()

    '''
    def __init__(self, R):
        '''
        TESTS::

            sage: FreeQuasiSymmetricFunctions(QQ)
            The combinatorial Hopf algebra of Free Quasi-Symmetric Functions\
            over the Rational Field
        '''
        from sage.categories.rings import Rings
        assert(R in Rings()), '%s must be a ring' % R
        Parent.__init__(
            self, base=R,
            category=Category.join((
                     GradedHopfAlgebras(R).WithRealizations(),
                     BidendriformBialgebras(R).WithRealizations()
            ))
        )

    F = Fundamental = LazyImport(
        "sage.combinat.cha._fqsym.fundamental_basis", "Fundamental")
    G = FundamentalDual = LazyImport(
        "sage.combinat.cha._fqsym.fundamental_dual_basis", "FundamentalDual")

    E = Elementary = LazyImport(
        "sage.combinat.cha._fqsym.elementary_basis", "Elementary")
    n = ElementaryDual = LazyImport(
        "sage.combinat.cha._fqsym.elementary_dual_basis", "ElementaryDual")

    H = Homogene = LazyImport(
        "sage.combinat.cha._fqsym.homogene_basis", "Homogene")
    m = HomogeneDual = LazyImport(
        "sage.combinat.cha._fqsym.homogene_dual_basis", "HomogeneDual")

    _shorthands = ["F", "G", "E", "n", "H", "m"]

    def _repr_(self):
        return "The combinatorial Hopf algebra of Free Quasi-Symmetric " + \
            "Functions over the %s" % self.base_ring()

    def dual(self):
        return self

    def a_realization(self):
        return self.Fundamental()

    class Bases(Category_realization_of_parent):

        def super_categories(self):
            R = self.base().base_ring()
            return [self.base().Realizations(),
                    GradedHopfAlgebrasWithBasis(R).Realizations(),
                    GradedModulesWithInternalProduct(R).Realizations(),
                    GradedAlgebrasWithDieseProduct(R).Realizations(),
                    BidendriformBialgebras.WithBasis(R).Realizations()]

        class ParentMethods:

            def build_morphisms(self):
                '''
                Define morphisms associated to the current basis
                '''

            def __init_extra__(self):
                self.build_morphisms()

            @abstract_method(optional=True)
            def expand_on_basis(self, sigma, i):
                '''
                This method return an iterator of all word associates to
                ``sigma`` in the polynomial realization

                MATH::

                    r(\mathbb{G}_{\sigma}) := \sum_{std(w) = \sigma} w

                TESTS::

                    sage: G = FQSym(QQ).G()
                    sage: sage: list(G.expand_on_basis([1,2],3))
                    [[1, 1], [1, 2], [1, 3], [2, 2], [2, 3], [3, 3]]
                '''

            def expand(self, elt, i):
                '''
                TESTS::

                    sage: G = FQSym(QQ).G()
                    sage: G.expand(G[1,2],3)
                    x_1^2 + x_1*x_2 + x_1*x_3 + x_2^2 + x_2*x_3 + x_3^2
                    sage: F = FQSym(QQ).F()
                    sage: F.expand(F[3,1,2],3)
                    x_2^2*x_1 + x_2*x_3*x_1 + x_3^2*x_1 + x_3^2*x_2
                '''
                from sage.algebras.free_algebra import FreeAlgebra
                FA = FreeAlgebra(self.base_ring(), i, ["x_%d" % i
                                                    for i in range(1, i + 1)])
                G = self.realization_of().G()
                return G.module_morphism(
                    lambda sigma: FA.sum(
                        [reduce(
                            lambda x, y: x * y,
                            [FA.gens_dict()["x_%d" % l] for l in w],
                            1
                        ) for w in G.expand_on_basis(sigma, i)]
                    ),
                    codomain=FA)(G(elt))

            @abstract_method
            def scalar_product_on_basis(self, sigma, mu):
                '''
                TESTS::

                    sage: F = FQSym(QQ).F()
                    sage: sigma = Permutation([1,4,2,3])
                    sage: F.scalar_product_on_basis(sigma, sigma.inverse())
                    1
                    sage: F.scalar_product_on_basis(sigma, sigma)
                    0
                '''

            def scalar_product(self, elt1, elt2):
                '''
                TESTS::

                    sage: FQS = FQSym(QQ)
                    sage: F = FQS.F(); G = FQS.G(); E = FQS.E()
                    sage: F[3,1,2].scalar_product(G[3,1,2])
                    1
                    sage: a = F[3, 1, 2] + F[3, 2, 1]
                    sage: b = F[2,3,1] + 2*F[3,2,1]
                    sage: a.scalar_product(b)
                    3
                    sage: E[3,1,2].scalar_product(E[2,1,3])
                    2
                '''
                F = self.realization_of().a_realization()
                return F.fundamental_scalar_product(F(elt1), F(elt2))

        class ElementMethods:

            def expand(self, i):
                '''
                the polynomial realization of ``self`` over `i` variables.
                TESTS::

                    sage: G = FQSym(QQ).G()
                    sage: G[1,2].expand(3)
                    x_1^2 + x_1*x_2 + x_1*x_3 + x_2^2 + x_2*x_3 + x_3^2
                '''
                return self.parent().expand(self, i)

            def scalar_product(self, Nelt):
                '''
                `\langle \mathbb{F}_{\sigma}, \mathbb{G}_{\mu}\rangle =
                 \delta_{\sigma, \mu}`
                EXAMPLES::

                    sage: F = FQSym(QQ).F()
                    sage: F[3,1,2].scalar_product(F[2,3,1])
                    1
                    sage: F[3,1,2].scalar_product(F[3,1,2])
                    0
                    sage: matr = lambda X,Y, n: matrix([[X(sigma).\
                    ....:                    scalar_product(Y(mu))\
                    ....:    for mu in Permutations(n)]\
                    ....:    for sigma in Permutations(n)])
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
                '''
                return self.parent().scalar_product(self, Nelt)

        class Base(GenericBasis):
            _basis_indices = Permutations()
            _prefix = "** TO DEFINE **"

            def _latex_term(self, m):
                """
                """
                if len(m) == 0:
                    return "1"
                prefix = self._print_options.get('latex_prefix')
                s = list(m)
                if max(s) < 10:
                    s = str(s).replace("[", "").\
                               replace("]", "").\
                               replace(",", "")
                else:
                    s = str(s).replace("[", "(").\
                               replace("]", ")")
                return "%s_{%s}" % (prefix, s)
