r"""
Knots

AUTHORS:

- Miguel Angel Marco Buzunariz
- Amit Jamadagni
"""

##############################################################################
#       Copyright (C) 2014   Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
##############################################################################

from sage.knots.link import Link
from sage.rings.finite_rings.integer_mod import Mod

class Knot(Link):
    """
    A knot.

    A knot is defined as embedding of the circle `\mathbb{S}^1` in the
    3-dimensional sphere `\mathbb{S}^3`, considered up to ambient isotopy.
    They represent the physical idea of a knotted rope, but with the
    particularity that the rope is closed. That is, the ends of the rope
    are joined.

    .. SEEALSO::

        :class:`Link`

    REFERENCES:

    - :wikipedia:`Knot_(mathematics)`

    .. TODO::

        - Implement the connect sum of two knots.
        - Make a class Knots for the monoid of all knots and have this be an
          element in that monoid.
    """
    def __init__(self, data, check=True):
        """
        Initialize ``self``.

        TESTS::

            sage: B = BraidGroup(8)
            sage: L = Link(B([-1, -1, -1, -2,1, -2, 3, -2]))
            sage: TestSuite(L).run()
            sage: L = Link(B([1, 2, 1]))
            sage: TestSuite(L).run()
            sage: L = Link([[1, 1, 2, 2]])
            sage: TestSuite(L).run()
        """
        Link.__init__(self, data)
        if check:
            if self.number_components() != 1:
                raise ValueError("the input has more than 1 connected component")

    def dt_code(self):
        """
        Return the DT code of ``self``.

        The DT code is generated by the following way:

        Start moving along the knot, as we encounter the crossings we
        start numbering them, so every crossing has two numbers assigned to
        it once we have traced the entire knot. Now we take the even number
        associated with every crossing.

        The following sign convention is to be followed:

        Take the even number with a negative sign if it is an overcrossing
        that we are encountering.

        OUTPUT: DT code representation of the knot

        EXAMPLES::

            sage: K = Knot([[1,5,2,4],[5,3,6,2],[3,1,4,6]])
            sage: K.dt_code()
            [4, 6, 2]
            sage: B = BraidGroup(4)
            sage: K = Knot(B([1, 2, 1, 2]))
            sage: K.dt_code()
            [4, -6, 8, -2]
            sage: K = Knot([[[1, -2, 3, -4, 5, -1, 2, -3, 4, -5]], [1, 1, 1, 1, 1]])
            sage: K.dt_code()
            [6, 8, 10, 2, 4]
        """
        b = self.braid().Tietze()
        N = len(b)
        label = [0 for i in range(2 * N)]
        string = 1
        next_label = 1
        type1 = 0
        crossing = 0
        while next_label <= 2 * N:
            string_found = False
            for i in range(crossing, N):
                if abs(b[i]) == string or abs(b[i]) == string - 1:
                    string_found = True
                    crossing = i
                    break
            if not string_found:
                for i in range(0, crossing):
                    if abs(b[i]) == string or abs(b[i]) == string - 1:
                        string_found = True
                        crossing = i
                        break
            assert label[2 * crossing + next_label % 2] != 1, "invalid knot"

            label[2 * crossing + next_label % 2] = next_label
            next_label = next_label + 1
            if type1 == 0:
                if b[crossing] < 0:
                    type1 = 1
                else:
                    type1 = -1
            else:
                type1 = -1 * type1
                if ((abs(b[crossing]) == string and b[crossing] * type1 > 0)
                    or (abs(b[crossing]) != string and b[crossing] * type1 < 0)):
                    if next_label % 2 == 1:
                        label[2 * crossing] = label[2 * crossing] * -1
            if abs(b[crossing]) == string:
                string = string + 1
            else:
                string = string - 1
            crossing = crossing + 1
        code = [0 for i in range(N)]
        for i in range(N):
            for j in range(N):
                if label[2 * j + 1] == 2 * i + 1:
                    code[i] = label[2 * j]
                    break
        return code

    def arf_invariant(self):
        """
        Return the Arf invariant.

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: K = Knot(B([-1, 2, 1, 2]))
            sage: K.arf_invariant()
            0
            sage: B = BraidGroup(8)
            sage: K = Knot(B([-2, 3, 1, 2, 1, 4]))
            sage: K.arf_invariant()
            0
            sage: K = Knot(B([1, 2, 1, 2]))
            sage: K.arf_invariant()
            1
        """
        a = self.alexander_polynomial()
        if ((Mod(a(-1), 8) == 1) or (Mod(a(-1), 8) == 7)):
            return 0

        return 1

