r"""
Signed Compositions
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
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

from composition import Compositions_n, Composition
import cartesian_product
from sage.rings.all import Integer
from sage.rings.arith import binomial
import __builtin__

class SignedCompositions(Compositions_n):
    """
    The class of signed compositions of `n`.

    EXAMPLES::

        sage: SC3 = SignedCompositions(3); SC3
        Signed compositions of 3
        sage: SC3.cardinality()
        18
        sage: len(SC3.list())
        18
        sage: SC3.first()
        [1, 1, 1]
        sage: SC3.last()
        [-3]
        sage: SC3.random_element()
        [1, -1, 1]
        sage: SC3.list()
        [[1, 1, 1],
         [1, 1, -1],
         [1, -1, 1],
         [1, -1, -1],
         [-1, 1, 1],
         [-1, 1, -1],
         [-1, -1, 1],
         [-1, -1, -1],
         [1, 2],
         [1, -2],
         [-1, 2],
         [-1, -2],
         [2, 1],
         [2, -1],
         [-2, 1],
         [-2, -1],
         [3],
         [-3]]

    TESTS::

        sage: SC = SignedCompositions(3)
        sage: TestSuite(SC).run()
    """
    def __repr__(self):
        """
        TESTS::

            sage: SignedCompositions(3)
            Signed compositions of 3
        """
        return "Signed compositions of %s"%self.n

    def __contains__(self, x):
        """
        TESTS::

            sage: [] in SignedCompositions(0)
            True
            sage: [0] in SignedCompositions(0)
            False
            sage: [2,1,3] in SignedCompositions(6)
            True
            sage: [-2, 1, -3] in SignedCompositions(6)
            True
        """
        if isinstance(x, __builtin__.list):
            for i in range(len(x)):
                if (not isinstance(x[i], (int, Integer))) and x[i] not in ZZ:
                    return False
                if x[i] == 0:
                    return False
        elif not isinstance(x, Composition):
            return False
        return sum([abs(i) for i in x]) == self.n

    def cardinality(self):
        r"""
        Return the number of elements in ``self``.

        The number of signed compositions of `n` is equal to

        .. MATH::

            \sum_{i=1}^{n+1} \binom{n-1}{i-1} 2^i

        TESTS::

            sage: SC4 = SignedCompositions(4)
            sage: SC4.cardinality() == len(SC4.list())
            True
            sage: SignedCompositions(3).cardinality()
            18
        """
        return sum([binomial(self.n-1, i-1)*2**(i) for i in range(1, self.n+1)])

    def __iter__(self):
        """
        TESTS::

            sage: SignedCompositions(0).list()   #indirect doctest
            [[]]
            sage: SignedCompositions(1).list()   #indirect doctest
            [[1], [-1]]
            sage: SignedCompositions(2).list()   #indirect doctest
            [[1, 1], [1, -1], [-1, 1], [-1, -1], [2], [-2]]
        """
        for comp in Compositions_n.__iter__(self):
            l = len(comp)
            a = [[1,-1] for i in range(l)]
            for sign in cartesian_product.CartesianProduct(*a):
                yield [ sign[i]*comp[i] for i in range(l)]

from sage.structure.sage_object import register_unpickle_override
register_unpickle_override('sage.combinat.composition_signed', 'SignedCompositions_n', SignedCompositions)

