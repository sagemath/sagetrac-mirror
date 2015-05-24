r"""
Non Negative Integer Semiring
"""
#*****************************************************************************
#  Copyright (C) 2010  Nicolas Borie <nicolas.borie at math.u-psud.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.sets.non_negative_integers import NonNegativeIntegers
from sage.categories.semirings import Semirings
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.sets.family import Family

class NonNegativeIntegerSemiring(NonNegativeIntegers):
    r"""
    A class for the semiring of the non negative integers

    This parent inherits from the infinite enumerated set of non
    negative integers and endows it with its natural semiring
    structure.

    EXAMPLES::

        sage: NonNegativeIntegerSemiring()
        Non negative integer semiring

    For convenience, ``NN`` is a shortcut for
    ``NonNegativeIntegerSemiring()``::

        sage: NN == NonNegativeIntegerSemiring()
        True

        sage: NN.category()
        Join of Category of semirings and Category of commutative monoids and Category of infinite enumerated sets and Category of facade sets

    Here is a piece of the Cayley graph for the multiplicative structure::

        sage: G = NN.cayley_graph(elements=range(9), generators=[0,1,2,3,5,7])
        sage: G
        Looped multi-digraph on 9 vertices
        sage: G.plot()
        Graphics object consisting of 48 graphics primitives

    This is the Hasse diagram of the divisibility order on ``NN``.

        sage: Poset(NN.cayley_graph(elements=[1..12], generators=[2,3,5,7,11])).show()

    Note: as for :class:`NonNegativeIntegers
    <sage.sets.non_negative_integers.NonNegativeIntegers>`, ``NN`` is
    currently just a "facade" parent; namely its elements are plain
    Sage ``Integers`` with ``Integer Ring`` as parent::

        sage: x = NN(15); type(x)
        <type 'sage.rings.integer.Integer'>
        sage: x.parent()
        Integer Ring
        sage: x+3
        18

    """
    def __init__(self):
        r"""
        TESTS::

            sage: NN = NonNegativeIntegerSemiring(); NN
            Non negative integer semiring
            sage: NN.category()
            Join of Category of semirings and Category of commutative monoids and Category of infinite enumerated sets and Category of facade sets
            sage: TestSuite(NN).run()
        """
        NonNegativeIntegers.__init__(self, category=(Semirings().Commutative(), InfiniteEnumeratedSets()) )

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: NonNegativeIntegerSemiring()
            Non negative integer semiring
        """
        return "Non negative integer semiring"

    def additive_semigroup_generators(self):
        r"""
        Returns the additive semigroup generators of ``self``.

        EXAMPLES::

            sage: NN.additive_semigroup_generators()
            Family (0, 1)
        """
        return Family([NN(0), NN(1)])

    def _latex_(self):
        r"""
        TESTS::

            sage: NN._latex_()
            '\\Bold{N}'
        """
        return '\\Bold{N}'

NN = NonNegativeIntegerSemiring()
