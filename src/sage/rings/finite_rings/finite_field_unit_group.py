r"""
Unit Groups of Finite Fields

EXAMPLES::

    sage: F.<a> = GF(81)
    sage: U = F.unit_group(); U
    Unit group with structure C80 of Finite Field in a of size 3^4

The generator is a primitive root of unity in the field::

    sage: U.gen()
    a
    sage: U.gen().multiplicative_order()
    80

Units in the field can be converted into elements of the unit group represented
as elements of an abstract multiplicative group::

    sage: U(1)
    1
    sage: U(-1)
    u^40
    sage: U(2 + a + a^3)
    u^11

Subgroups of the unit group, consisting of the group of n-th roots of unity can also be constructed::

    sage: U5 = F.unit_group(5)
    sage: list(U5)
    [1, 2*a^2 + a + 2, 2*a^3 + a + 2, 2*a^2 + 2*a + 1, a^3 + 2*a^2 + 2*a]

The log function gives the exponent of a unit as a power of the chosen generator or of a given base::

    sage: U.log(1 + a)
    28
    sage: U.gen()^28
    a + 1
    sage: U.log(a^2, a + a^3)
    38
    sage: (a + a^3)^38
    a^2

AUTHOR:

- Francis Clarke
"""
# ****************************************************************************
#       Copyright (C) 2009 William Stein, Francis Clarke
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
# ****************************************************************************

from sage.groups.abelian_gps.abelian_group import AbelianGroup_class
from sage.groups.abelian_gps.abelian_group_element import AbelianGroupElement
from sage.arith.misc import gcd, power_mod
from sage.structure.sequence import Sequence


class FiniteFieldUnitGroup(AbelianGroup_class):
    """
    The unit group of a finite field.  Based on the implementation
    of the unit groups of number fields by John Cremona.

    Author: Francis Clarke
    """

    def __init__(self, field, n=None):
        """
        Create the unit group of a finite field,
        or the subgroup of n-th roots of unity.

        INPUT:

        - ``field`` - a finite field

        - ``n`` - a positive integer dividing the order of the field,
           or (default) None


        EXAMPLES::

            sage: GF(5).unit_group()
            Unit group with structure C4 of Finite Field of size 5
            sage: GF(25, 'a').unit_group(8)
            Group (with structure C8) of 8th roots of unity in Finite Field in a of size 5^2
        """
        # compute a generator
        z = field.zeta(n)
        if n is None:
            n = field.order() - 1

        # Store the field and the generator:
        self.__field = field
        self.__gens = Sequence([z], immutable=True, universe=self, check=False)

        # Construct the abtract group:
        AbelianGroup_class.__init__(self, (n,), 'u')

    def __call__(self, u):
        """
        Return the abstract group element corresponding to the unit u.

        INPUT:

        - ``u`` -- Any object from which an element of the unit group's
          field `K` may be constructed; an error is raised if
          an element of `K` cannot be constructed from `u`, or if the
          element constructed is not a  unit.

        EXAMPLES::

            sage: F25.<a> = GF(25)
            sage: mu8 = F25.unit_group(8)
            sage: mu8(1 + 3*a)
            u^3
        """
        K = self.__field
        q = K.order()
        n = self.order()

        try:
            u = K(u)
        except TypeError:
            raise ValueError("%s is not an element of %s" % (u, K))
        if u == 0:
            raise ValueError("%s is not a unit" % u)
        if n != q - 1 and u**n != 1:
            raise ValueError("%s is not a %s root of unity" % (u, n.ordinal_str()))
        m = u.log(self.gen())  # Do a discrete logarithm.
        return AbelianGroupElement(self, [m])

    def __contains__(self, u):
        """
        EXAMPLES::

            sage: 3 in GF(49, 'z').unit_group()
            True
            sage: 8 in GF(73).unit_group(12)
            True
            sage: 14 in GF(7).unit_group()
            False
        """
        K = self.__field
        q = K.order()
        n = self.order()

        try:
            u = K(u)
        except TypeError:
            return False
        if u == 0:
            return False
        return n == q - 1 or u**n == 1

    def _coerce_impl(self, x):
        """
        Canonical coercion of ``x`` into this unit group.
        """
        return self(x)

    def gens(self):
        """
        Return a generator for the unit group, as a list of length one.

        EXAMPLES::

            sage: GF(41).unit_group().gens()
            [6]
            sage: GF(49, 'a').unit_group(4).gens()
            [6*a + 4]
        """
        return self.__gens

    def ngens(self):
        """
        Return the number of generators of the unit group.

        EXAMPLES::

            sage: GF(11).unit_group().ngens()
            1
        """
        return 1

    def gen(self, i=0):
        """
        Return the ``i``-th generator for this unit group,
        where ``i`` had better be 0.

        EXAMPLES::

            sage: GF(41).unit_group().gen()
            6
        """
        if i < 0 or i >= len(self.__gens):
            raise IndexError
        return self.__gens[i]

    def _repr_(self):
        """
        Return the string representation of this unit group.

        EXAMPLES::

            sage: GF(31).unit_group()
            Unit group with structure C30 of Finite Field of size 31
            sage: GF(31).unit_group(10)
            Group (with structure C10) of 10th roots of unity in Finite Field of size 31
        """
        K = self.__field
        q = K.order()
        n = self.order()
        if n == q - 1:
            return 'Unit group with structure %s of %s' % (
                self._group_notation(self.invariants()), K)
        elif n == 1:
            return 'Trivial multiplicative group {1} in %s' % K
        else:
            if n == 2:
                nth = 'square'
            elif n == 3:
                nth = 'cube'
            else:
                nth = n.ordinal_str()
            return 'Group (with structure %s) of %s roots of unity in %s' % (
                self._group_notation(self.invariants()), nth, K)

    def field(self):
        """
        Return the finite field of which this is the group of units.

        EXAMPLES::

            sage: GF(13).unit_group().field()
            Finite Field of size 13
        """
        return self.__field

    def list(self):
        """
        List all the element of the unit group in the order given by taking powers
        of the generator.

        EXAMPLES::

            sage: GF(13).unit_group().list()
            [1, 2, 4, 8, 3, 6, 12, 11, 9, 5, 10, 7]
            sage: GF(23).unit_group(11).list()
            [1, 2, 4, 8, 16, 9, 18, 13, 3, 6, 12]
        """
        return [self.gen()**i for i in range(self.order())]

    def __iter__(self):
        """
        Return an iterator for the set of elements of the unit group

        EXAMPLES::

            sage: [x for x in GF(13).unit_group()]
            [1, 2, 4, 8, 3, 6, 12, 11, 9, 5, 10, 7]
            sage: [x^6 for x in GF(25, 'a').unit_group()]
            [1, 2, 4, 3, 1, 2, 4, 3, 1, 2, 4, 3, 1, 2, 4, 3, 1, 2, 4, 3, 1, 2, 4, 3]
        """
        for i in range(self.order()):
            yield self.gen()**i

    def log(self, u, base=None):
        """
        Return the exponents of the unit ``u`` with respect to the
        group generator, or to the given base.

        INPUT:

        - ``u`` -- Any object from which an element of the unit group's
          field `K` may be constructed; an error is raised if an element of `K`
          cannot be constructed from u, or if the element constructed is not a
          unit.

        - ``base`` -- a generator of the group, or any object from which such
          an element may be constructed.  For the default value of None the
          group's generator is used.

        OUTPUT:

        the exponent of ``u`` with respect to the unit group's generator
        or the given base.

        EXAMPLES::

            sage: U = GF(23).unit_group()
            sage: U.gen()
            5
            sage: U.log(7)
            19
            sage: power_mod(5, 19, 23) == 7
            True
            sage: F49.<g> = GF(49)
            sage: mu16 = F49.unit_group(16)
            sage: mu16.log(g + 3, 4*g + 3)
            12
            sage: (4*g + 3)^12
            g + 3
        """
        if base is None:
            return self(u).list()[0]
        else:
            n = self.order()
            k = self.log(base)
            if gcd(n, k) != 1:
                raise ValueError("%s is not a generator" % base)
            return (self.log(u) * power_mod(k, -1, n)).mod(n)
