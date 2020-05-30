r"""
Oscillating Tableaux

This is an implementation of the abstract base class
:class:`sage.combinat.path_tableau.path_tableaux`.

In this implementation we have oscillating tableaux which are sequences
of partitions such that at each step either one box is added or one box
is removed.

The purpose is to define promotion, evacuation and the action of
the cactus group on oscillating tableaux.

For further information see [RSW2013]_ and [PRW2018]_ .

AUTHORS:
- Bruce Westbury (2018): initial version
"""

#*****************************************************************************
#       Copyright (C) 2018 Bruce Westbury <bruce.westbury@gmail.com>,
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from six import add_metaclass
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.combinat.perfect_matching import PerfectMatching
from sage.combinat.path_tableaux.path_tableau import PathTableau, PathTableaux
from sage.combinat.tableau import Tableau
from sage.combinat.partition import Partition
from sage.rings.integer import Integer
from itertools import zip_longest
r"""
Here we illustrate one of the main theorems of [PRW2018]_ that
promotion for oscillating tableaux corresponds to inverse rotation of perfect matchings.

    TESTS::

        sage: all(OscillatingTableau(pm) == OscillatingTableau(pm.rotate()).promotion() for pm in PerfectMatchings(6))
        True


REFERENCES:

.. [RSW2013] Martin Rubey, Bruce Sagan, Bruce W. Westbury.
   *Descent sets for symplectic groups *,
   :arxiv:`1303.5850`

.. [PRW2018] Stephan Pfannerer, Martin Rubey, Bruce W. Westbury.
   *Promotion on oscillating and alternating tableaux and rotation of matchings and permutations *,
   :arxiv:`1804.06736`

EXAMPLES::

    sage: ot = OscillatingTableau([1,2,-2,-1])
    sage: CylindricalDiagram(ot)
     [[], [1], [1, 1], [1], []]
     ['', [], [1], [1, 1], [1], []]
     ['', '', [], [1], [1, 1], [1], []]
     ['', '', '', [], [1], [1, 1], [1], []]
     ['', '', '', '', [], [1], [1, 1], [1], []]

    sage: TestSuite(ot).run()

    sage: t = OscillatingTableau([1,2,1,-1,-2,-1])
    sage: len(t.orbit())
    9
    sage: TestSuite(t).run()
"""

@add_metaclass(InheritComparisonClasscallMetaclass)
class OscillatingTableau(PathTableau):
    r"""
    An oscillating tableau is a sequence of partitions such that at
    each step either a single box is added or a single box is removed.
    These are also known as up-down tableaux.

    Usually the sequence will start and end with the empty partition
    but this is not a requirement.

    These arise in the representation theory of the symplectic groups.
    Let 'C' be the crystal of the vector representation of 'Sp(2n)'.
    This has '2n' vertices which we label '1, 2,\ldots ,n, -n,\ldots ,-1'.
    Consider the 'r'-th tensor power of 'C'. This has vertices words of
    length 'r' in this alphabet. The analogue of the Robinson-Schensted
    correspondence is a bijection between these words and pairs
    '(P,Q)' where 'P', the insertion tableau, is a symplectic tableau
    and 'Q', the recording tableau, is an oscillating tableau.
    The oscillating tableau has initial shape the empty partition
    and the final shape is the shape of 'P'.

    There is a bijection between oscillating tableaux whose initial
    and final shape is the empty partition and perfect matchings.
    This is called the Sundaram bijection.

    INPUT:

    - a sequence of partitions
    - a sequence of lists
        each list is converted to a partition
    - a sequence of nonzero integers
        this is converted to a sequence of partitions by starting with
        the empty partition and at each step adding a box in row i
        or removing a box in row -i
    - a perfect matching
        this is converted to a sequence of partitions by the
        Sundaram bijection


    EXAMPLES::

        sage: OscillatingTableau([Partition([]), Partition([1]), Partition([])])
        [[], [1], []]

        sage: OscillatingTableau([[], [1], []])
        [[], [1], []]

        sage: pm = PerfectMatching([[1,2]])
        sage: OscillatingTableau(pm)
        [[], [1], []]

        sage: OscillatingTableau([1,-1])
        [[], [1], []]

    """

    @staticmethod
    def __classcall_private__(self, ot):
        """
        This is the preprocessing for creating paths.

        INPUT:

        - a sequence of partitions
        - a sequence of lists
        - a sequence of nonzero integers
        - a perfect matching

        EXAMPLES::

            sage: pm = PerfectMatching([[1,2]])
            sage: OscillatingTableau(pm)
            [[], [1], []]

            sage: pm = PerfectMatching([[1,3],[2,4]])
            sage: OscillatingTableau(pm)
            [[], [1], [1, 1], [1], []]

            sage: OscillatingTableau([1,-1])
            [[], [1], []]

        TESTS::

            sage: OscillatingTableau([1,0,-1])
            Traceback (most recent call last):
            ...
            ValueError: list may not contain zero

            sage: t = OscillatingTableau([1,1,2,2,-1])
            Traceback (most recent call last):
            ...
            ValueError: [1, 2] is not an element of Partitions

            sage: OscillatingTableau([Permutation(1,2),Permutation(2,1)])
            Traceback (most recent call last):
            ...
            TypeError: 'sage.rings.integer.Integer' object is not iterable
        """

        if isinstance(ot, OscillatingTableau):
            return ot

        w = None

        if isinstance(ot,(tuple,list)):
            try:
                w = tuple([ Partition(a) for a in ot ])
            except TypeError:
                pass

            try:
                ot = tuple([Integer(a) for a in ot])
                if any([a==0 for a in ot]):
                    raise ValueError("list may not contain zero")

                w = [Partition([])]*(len(ot)+1)

                for i,a in enumerate(ot):
                    if a > 0:
                        w[i+1] = w[i].add_cell(a-1)
                    else:
                        pt = list(w[i])
                        pt[abs(a)-1] -= 1
                        w[i+1] = Partition(pt)
            except TypeError:
                pass

        if isinstance(ot,PerfectMatching):
            tb = Tableau([])
            n = ot.size()
            w = [Partition([])]*(n+1)

            for i in range(n,0,-1):
                c = [ k for k, r in enumerate(tb) if i in r ]
                if len(c) == 1:
                    k = c[0]
                    tc = [ list(x) for x in tb ]
                    tc[k].remove(i)
                    if tc[k] == []:
                        tc.remove(tc[k])
                    tb = Tableau(tc)
                    w[i-1] = tb.shape()
                elif len(c) == 0:
                    x = ot.partner(i)
                    tb = tb.bump(x)
                    w[i-1] = tb.shape()
                else:
                    raise RuntimeError("tableau must be standard")

        if w is None:
            raise ValueError("invalid input {!s}".format(ot))

        return OscillatingTableaux()(w)

    def check(self):
        r"""
        Check that ``self`` is a valid oscillating tableau.

        EXAMPLES::

            sage: ot = OscillatingTableau([[], [1], [2], [1], [], [1], []])
            sage: ot.check()

        TESTS::

            sage: ot = OscillatingTableau([[], [1], [3], [1], [], [1], []])
            Traceback (most recent call last):
            ...
            ValueError: adjacent partitions differ by more than one box
            sage: ot = OscillatingTableau([[], [1], [2], [2,1], [4]])
            Traceback (most recent call last):
            ...
            ValueError: next partition is not obtained by adding a cell
            sage: ot = OscillatingTableau([[], [1], [2], [3], [1,1]])
            Traceback (most recent call last):
            ...
            ValueError: next partition is not obtained by removing a cell
        """
        n = len(self)
        for i in range(n-1):
            h = self[i]
            t = self[i+1]

            if abs(t.size()-h.size()) != 1:
                    raise ValueError("adjacent partitions differ by more than one box")
            if t.size() == h.size()+1:
                if not t.contains(h):
                    raise ValueError("next partition is not obtained by adding a cell")
            if h.size() == t.size()+1:
                if not h.contains(t):
                    raise ValueError("next partition is not obtained by removing a cell")

    def local_rule(self,i):
        r"""
        This has input a list of objects. This method first takes
        the list of objects of length three consisting of the `(i-1)`-st,
        `i`-th and `(i+1)`-term and applies the rule. It then replaces
        the `i`-th object  by the object returned by the rule.

        EXAMPLES::

            sage: t = OscillatingTableau([[],[1],[2],[2,1],[1,1]])
            sage: t.local_rule(3)
            [[], [1], [2], [1], [1, 1]]

        TESTS::

            sage: t = OscillatingTableau([[],[1],[2],[1],[]])
            sage: t.local_rule(0)
            Traceback (most recent call last):
            ...
            ValueError: 0 is not a valid integer
            sage: t._local_rule(4)
            Traceback (most recent call last):
            ...
            ValueError: 4 is not a valid integer
        """
        def _rule(x):
            """
            This is the rule on a sequence of three partitions.
            """
            result = [ i-j+k for i,j,k in zip_longest(*x,fillvalue=0)]
            result = [ abs(_) for _ in result ]
            result.sort(reverse=True)
            return Partition(result)

        if not (i > 0 and i < len(self)-1 ):
            raise ValueError("%d is not a valid integer" % i)

        with self.clone() as result:
            result[i] = _rule(self[i-1:i+2])

        return result

    def is_skew(self):
        r"""
        Return ``True`` if Tableau is skew and ``False`` if not.

        EXAMPLES::

            sage: T = OscillatingTableau([[],[1],[2],[1],[]])
            sage: T.is_skew()
            False

            sage: T = OscillatingTableau([[1],[2],[1],[]])
            sage: T.is_skew()
            True
        """
        return self[0] != Partition([])

    def conjugate(self):
        r"""
        Return the oscillating tableau obtained by conjugating all partitions in ''self''.

        EXAMPLES::

            sage: OscillatingTableau([[],[1],[2],[1],[]]).conjugate()
            [[], [1], [1, 1], [1], []]
        """
        return OscillatingTableau([p.conjugate() for p in self])

    def crossing_number(self):
        """
        Return the crossing number.

        EXAMPLES::

            sage: T = OscillatingTableau([[],[1],[2],[1],[]])
            sage: T.crossing_number()
            1

        TESTS::

            sage: [ OscillatingTableau(pm).crossing_number() for pm in PerfectMatchings(6) ]
            [1, 2, 1, 2, 1, 2, 2, 3, 2, 2, 1, 2, 2, 2, 1]
        """
        return max( a.length() for a in self )

    def nesting_number(self):
        """
        Return the nesting number.

        EXAMPLES::

            sage: T = OscillatingTableau([[],[1],[2],[1],[]])
            sage: T.nesting_number()
            2

        TESTS::

            sage: [ OscillatingTableau(pm).nesting_number() for pm in PerfectMatchings(6) ]
            [1, 1, 2, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2, 2, 3]
        """

        v = (a for a in self if a.length() > 0)
        if v == None:
            return 0
        else:
            return max( a[0] for a in v )

    def to_perfect_matching(self):
        """
        Construct the perfect matching.

        EXAMPLES::

            sage: OscillatingTableau([1,2,-2,-1]).to_perfect_matching()
            [(1, 3), (2, 4)]

        TESTS::

            sage: all(pm == OscillatingTableau(pm).to_perfect_matching()  for pm in PerfectMatchings(6) )
            True

            sage: OscillatingTableau([[], [1], [1, 1], [2, 1], [1, 1]]).to_perfect_matching()
            Traceback (most recent call last):
            ...
            ValueError: the final shape of [[], [1], [1, 1], [2, 1], [1, 1]] must be the empty shape
        """
        if self[-1]:
            raise ValueError("the final shape of {!s} must be the empty shape".format(self))

        return self.sundaram()[1]

    def to_word(self):
        """
        Converts an oscillating tableau to a word in the alphabet
        ...,-2,-1,1,2,...

        EXAMPLES::

            sage: T = OscillatingTableau([[],[1],[2],[1],[]])
            sage: T.to_word()
            [1, 1, -1, -1]
            sage: T = OscillatingTableau([[],[1],[1,1],[1],[]])
            sage: T.to_word()
            [1, 2, -2, -1]
            sage: OscillatingTableau([[2,1],[2,2]]).to_word()
            [2]

        TESTS::

            sage: ots = [ OscillatingTableau(pm) for pm in PerfectMatchings(6) ]
            sage: all( OscillatingTableau(ot.to_word())==ot for ot in ots )
            True
        """
        n = len(self)
        result = [0]*(n-1)
        l = [ len(_) for _ in self ]

        for i in range(n-1):
            if l[i] > l[i+1]:
                result[i] = -l[i]
            elif l[i] < l[i+1]:
                result[i] = l[i+1]
            else:
                for j in range(l[i]):
                    d = self[i+1][j]-self[i][j]
                    if d:
                        result[i] = (j+1)*d
                        break

        return result

    def descents(self):
        r"""
        Return the descent set. This is defined in [RSW2013]_

        It is defined in terms of the word 'w' by 'i' is a descent if
        'w[i-1] > 0' and 'w[i] < 0' or '0 < w[i-1] < w[i]' or 'w[i-1] < w[i] < 0'.

        Equivalently, it is defined in terms of the sequence of partitions by:
        'i' is a descent if
        * step 'i' is an expansion and step 'i+1' is a a contraction, or
        * steps 'i' and 'i+1' are both expansions and the row for the box added
        at step 'i' is strictly above the row for the box added
        at step 'i+1', or
        * steps 'i' and 'i+1' are both contractions and the row for the box added
        at step 'i' is strictly below the row for the box added
        at step 'i+1'

        EXAMPLES::

            sage: T = OscillatingTableau([[],[1],[2],[1],[]])
            sage: T.descents()
            {2}

        TESTS::

            all(OscillatingTableau(pm).conjugate().descents()
                == set(pm.to_permutation().inverse().descents())
                    for pm in PerfectMatchings(6))
            True
        """
        result = set()
        w = self.to_word()

        for i, x in enumerate(zip(w,w[1:])):
            u, v = x
            if (u > 0 and v < 0) or (0 < u < v) or (u < v < 0):
                result.add(i+1)

        return result

    def sundaram(self):
        r"""
        This implements the bijection due to S. Sundaram between
        oscillating tableaux with empty initial shape and pairs
        '(S,M)' where 'S' is a partial standard tableaux whose shape
        is the final shape of the oscillating tableau and a
        perfect matching on the complement of the set of entries
        of 'S'.

        In particular if the final shape of the oscillating tableau
        is empty this gives the empty tableau and a perfect matching.
        This is an inverse bijection to the map which constructs an
        oscillating tableau from a perfect matching.

        INPUT: An oscillating tableau which starts with the empty partition.

        OUTPUT: A pair '(S,M)'; 'S' is a Tableau and 'M' is a PerfectMatching

        EXAMPLES::

            sage: t = OscillatingTableau([[],[1],[1,1],[1],[]])
            sage: t.sundaram()
            ([], [(1, 3), (2, 4)])

            sage: t = OscillatingTableau([[],[1],[1,1],[2,1],[2],[1],[2],[2,1],[2,1,1],[2,1]])
            sage: t.sundaram()
            ([[2, 7], [8]], [(1, 4), (3, 5), (6, 9)])

            sage: s = OscillatingTableau([[],[1],[2],[2,1],[1,1],[1],[1,1],[2,1],[3,1],[2,1]])
            sage: s.sundaram()
            ([[3, 7], [6]], [(1, 5), (2, 4), (8, 9)])

        TESTS::

            sage: T = OscillatingTableau([[1],[2],[1],[]])
            sage: T.sundaram()
            Traceback (most recent call last):
            ...
            ValueError: this has only been implemented for straight oscillating tableaux

            sage: T = OscillatingTableau([[],[1],[2],[1]])
            sage: T.sundaram()
            ([[1]], [(2, 3)])

            sage: all(OscillatingTableau(pm).sundaram()[1]==pm for pm in PerfectMatchings(6))
            True
        """
        if self.is_skew():
            raise ValueError("this has only been implemented for straight oscillating tableaux")
        tb = Tableau([])
        pm = set([])

        for i, u in enumerate(zip(self,self[1:])):
            mu, lb = u
            if lb.contains(mu):
                cell = [ c for c in lb.corners() if c not in mu.corners() ][0]
                tb = tb.add_entry(cell,i+1)
            else:
                cell = [ c for c in mu.corners() if c not in lb.corners() ][0]
                tb, x = tb.reverse_bump(cell)
                pm.add((i+1,x))

        return (Tableau(tb),PerfectMatching(pm))

    def _test_descents(self, **options):
        r"""
        Check the property of descent sets proved in [RSW2013]_
        Namely that the Sundaram map preserves descent sets.

        TESTS::

            sage: t = OscillatingTableau([1,1,2,3,-1,-3,-2,-1])
            sage: t._test_descents()
        """
        tester = self._tester(**options)

        lhs = self.descents()
        tb, pm = self.conjugate().sundaram()
        rhs = set(tb.standardization().standard_descents())
        rhs = rhs.union(pm.to_permutation().inverse().descents())
        wt = tb.weight()
        for i, x in enumerate(zip(wt,wt[1:])):
            if x[0]>x[1]:
                rhs.add(i+1)

        tester.assertTrue(lhs == rhs)

###############################################################################

class OscillatingTableaux(PathTableaux):
    """
    This is the parent class for oscillating tableaux.

    An oscillating tableau is a sequence of partitions such that at
    each step either a single box is added or a single box is removed.
    These are also know as up-down tableaux.

    This class can be called to construct an oscillating tableau directly.
    The input is a sequence of partitions.

    EXAMPLES:

        sage: OscillatingTableaux()([Partition([]), Partition([1]), Partition([])])
        [[], [1], []]

    """
    Element = OscillatingTableau

