r"""
Examples of graded connected Hopf algebras with basis
"""
#*****************************************************************************
#  Copyright (C) 2015 Jean-Baptiste Priez <jbp@kerios.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.graded_hopf_algebras_with_basis import GradedHopfAlgebrasWithBasis
from sage.combinat.free_module import CombinatorialFreeModule
from sage.functions.other import binomial
from sage.misc.cachefunc import cached_method
from sage.sets.non_negative_integers import NonNegativeIntegers


class GradedConnectedCombinatorialHopfAlgebraWithPrimitiveGenerator(CombinatorialFreeModule):
    r"""
    This class illustrates an implementation of a graded Hopf algebra
    with basis that has one primitive generator of degree 1 and basis
    elements indexed by non-negative integers.

    This Hopf algebra example differs from what topologists refer to as
    a graded Hopf algebra because the twist operation in the tensor rule
    satisfies

    .. MATH::

        (\mu \otimes \mu) \circ (id \otimes \tau \otimes id) \circ
        (\Delta \otimes \Delta) = \Delta \circ \mu

    where `\tau(x\otimes y) = y\otimes x`.

    """
    def __init__(self, base_ring):
        """
        EXAMPLES::

            sage: from sage.categories.examples.graded_connected_hopf_algebras_with_basis import Example
            sage: H = Example(QQ)
            sage: TestSuite(H).run()

        """
        CombinatorialFreeModule.__init__(self, base_ring, NonNegativeIntegers(),
                                         category=GradedHopfAlgebrasWithBasis(base_ring).Connected())


    @cached_method
    def one_basis(self):
        """
        Returns 0, which index the unit of the Hopf algebra.

        OUTPUT:

        - the non-negative integer 0

        EXAMPLES::

            sage: from sage.categories.examples.graded_connected_hopf_algebras_with_basis import Example
            sage: H = Example(QQ)
            sage: H.one_basis()
            0
            sage: H.one()
            P0

        """
        return self.basis().keys()(0)

    def degree_on_basis(self, i):
        """
        The degree of a non-negative integer is itself

        INPUT:

        - ``i`` -- a non-negative integer

        OUTPUT:

        - a non-negative integer

        TESTS::

            sage: from sage.categories.examples.graded_connected_hopf_algebras_with_basis import Example
            sage: H = Example(QQ)
            sage: H.degree_on_basis(45)
            45

        """
        return i

    def _repr_(self):
        """
        Representation of the graded connected Hopf algebra

        EXAMPLES::

            sage: from sage.categories.examples.graded_connected_hopf_algebras_with_basis import Example
            sage: Example(QQ)
            An example of a graded connected Hopf algebra with basis over Rational Field

        """
        return "An example of a graded connected Hopf algebra with basis over %s" % self.base_ring()

    def _repr_term(self, i):
        """
        Representation for the basis element indexed by the integer ``i``.

        EXAMPLES::

            sage: from sage.categories.examples.graded_connected_hopf_algebras_with_basis import Example
            sage: H = Example(QQ)
            sage: H._repr_term(45)
            'P45'
        """
        return 'P' + repr(i)

    def product_on_basis(self, i, j):
        """
        The product of two basis elements.

        The product of elements of degree ``i`` and ``j`` is an element
        of degree ``i+j``.

        INPUT:

        - ``i``, ``j`` -- non-negative integers

        OUPUT:

        - a basis element indexed by ``i+j``

        TESTS::

            sage: from sage.categories.examples.graded_connected_hopf_algebras_with_basis import Example
            sage: H = Example(QQ)
            sage: H.monomial(4) * H.monomial(5)
            P9

        """
        return self.monomial(i+j)

    def coproduct_on_basis(self, i):
        """
        The coproduct of a basis element.

        .. MATH::

            \Delta(P_i) = \sum_{j=0}^i P_{i-j} \otimes P_j

        INPUT:

        - ``i`` -- a non-negative integer

        OUTPUT:

        - an element of the tensor square of ``self``

        TESTS::

            sage: from sage.categories.examples.graded_connected_hopf_algebras_with_basis import Example
            sage: H = Example(QQ)
            sage: H.monomial(3).coproduct()
            P0 # P3 + 3*P1 # P2 + 3*P2 # P1 + P3 # P0

        """
        return self.sum_of_terms(
            ((i-j, j), binomial(i, j))
            for j in range(i+1)
        )

Example = GradedConnectedCombinatorialHopfAlgebraWithPrimitiveGenerator
