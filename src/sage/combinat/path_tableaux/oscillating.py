r"""
Oscillating Tableaux

This is an implementation of the abstract base class
:class:`sage.combinat.path_tableau.path_tableaux`.

In this implementation we have oscillating tableaux which are sequences
of partitions such that at each step either one box is added or one box
is removed.

The purpose is to define promotion, evacuation and the action of
the cactus group on oscillating tableaux.

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
from sage.structure.list_clone import ClonableArray
from sage.structure.parent import Parent
from sage.combinat.perfect_matching import PerfectMatching
from sage.combinat.path_tableaux.path_tableau import PathTableau, PathTableaux
from sage.combinat.tableau import Tableau
from sage.combinat.partition import Partition
from sage.modules.free_module_element import vector
from sage.rings.integer import Integer

"""
Here we illustrate one of the main theorems of [PRW2018]_ that
promotion for oscillating tableaux corresponds to inverse rotation of perfect matchings.

    TEST::

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

    sage: t = OscillatingTableau([1,2,-2,-1])
    sage: SkewTableau(t.cylindrical_diagram()).pp()
     [][1][1, 1][1] []
      . [][1][1, 1][1] []
      .  . [][1][1, 1][1] []
      .  .  . [][1][1, 1][1] []
      .  .  .  . [][1][1, 1][1] []

    sage: TestSuite(t).run()

    sage: t = OscillatingTableau([1,2,1,-1,-2,-1])
    sage: len(t.orbit())
    9
    sage: TestSuite(t).run()
"""

@add_metaclass(InheritComparisonClasscallMetaclass)
class OscillatingTableau(ClonableArray,PathTableau):
    """
    An oscillating tableau is a sequence of partitions such that at
    each step either a single box is added or a single box is removed.
    These are also known as up-down tableaux.

    Usually the sequence will start and end with the empty partition
    but this is not a requirement.

    These arise in the representation theory of the symplectic groups.
    Let C be the crystal of the vector representation of Sp(2n).
    This has 2n vertices which we label 1, 2, ... ,n, -n, ... ,-1.
    Consider the r-th tensor power of C. This has vertices words of
    length r in this alphabet. The analogue of the Robinson-Schensted
    correspondence is a bijection between these words and pairs
    (P,Q) where P, the insertion tableau, is a symplectic tableau
    and Q, the recording tableau, is an oscillating tableau.
    The oscillating tableau has initial shape the empty partition
    and the final shape is the shape of P.

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
                    raise ValueError("List may not contain zero.")

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
                    tc = map(list,tb)
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
            raise ValueError("invalid input %s" % ot)

        return OscillatingTableaux()(w)

    def _hash_(self):
        return hash(tuple(map(tuple, self)))

    def check(self):
        n = len(self)
        for i in range(n-1):
            h = self[i]
            t = self[i+1]

            if abs(t.size()-h.size()) != 1:
                    raise ValueError("Adjacent partitions differ by more than one box.")
            if t.size() == h.size()+1:
                if not t.contains(h):
                    raise ValueError("Next partition is not obtained by adding a cell.")
            if h.size() == t.size()+1:
                if not h.contains(t):
                    raise ValueError("Next partition is not obtained by removing a cell.")

    def _local_rule(self,i):
        """
        This has input a list of objects. This method first takes
        the list of objects of length three consisting of the `(i-1)`-st,
        `i`-th and `(i+1)`-term and applies the rule. It then replaces
        the `i`-th object  by the object returned by the rule.

        EXAMPLES::

            sage: t = CatalanTableau([0,1,2,3,2,1,0])
            sage: t._local_rule(3)
            [0, 1, 2, 1, 2, 1, 0]
        """

        def _rule(x):
            """
            This is the rule on a sequence of three partitions.
            """
            y = map(list,x)
            m = max([ len(u) for u in y ])
            z = map( lambda u: vector(u + [0]*(m-len(u)) ), y )
            result = z[0]-z[1]+z[2]
            result = map(abs,result)
            result.sort(reverse=True)
            return Partition(result)

        if not (i > 0 and i < len(self) ):
            raise ValueError("%d is not a valid integer" % i)

        with self.clone() as result:
            result[i] = _rule(self[i-1:i+2])

        return result

    def is_skew(self):
        """
        Returns True if Tableau is skew and False if not.

        EXAMPLE:
        sage: T = OscillatingTableau([[],[1],[2],[1],[]])
        sage: T.is_skew()
        False
        """
        return self[0] != Partition([])

    def crossing_number(self):
        """
        Returns the crossing number.
        sage: T = OscillatingTableau([[],[1],[2],[1],[]])
        sage: T.crossing_number()
        1
        """
        return max( a.length() for a in self )

    def nesting_number(self):
        """
        Returns the nesting number.

        EXAMPLE:
        sage: T = OscillatingTableau([[],[1],[2],[1],[]])
        sage: T.nesting_number()
        2
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

        """
        return self.sundaram()[1]

    def to_word(self):
        """
        Converts an oscillating tableau to a word in the alphabet
        ...,-2,-1,1,2,...

        EXAMPLE::

            sage: T = OscillatingTableau([[],[1],[2],[1],[]])
            sage: T.to_word()
            [1, 1, -1, -1]
            sage: T = OscillatingTableau([[],[1],[1,1],[1],[]])
            sage: T.to_word()
            [1, 2, -2, -1]
            sage: OscillatingTableau([[2,1],[2,2]]).to_word()
            [2]

        """
        n = len(self)
        result = [0]*(n-1)
        l = map(len, self)

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
        """
        Returns the descent set. This is defined in [RSW2013]_

        EXAMPLE:
        sage: T = OscillatingTableau([[],[1],[2],[1],[]])
        sage: T.descents()
        {2}

        """
        result = set()
        w = self.to_word()

        for i in range(1,len(w)):
            if w[i-1] > 0 and w[i] < 0:
                result.add(i)
            if 0 < w[i-1] < w[i]:
                result.add(i)
            if w[i-1] < w[i] < 0:
                result.add(i)

        return result

    def sundaram(self):
        """
        This implements the bijection due to S. Sundaram between
        oscillating tableaux with empty initial shape and pairs
        (S,M) where S is a partial standard tableaux whose shape
        is the final shape of the oscillating tableau and a
        perfect matching on the complement of the set of entries
        of S.

        INPUT: A straight oscillating tableau.

        OUTPUT: A pair (S,M); S is a Tableau and M is a PerfectMatching

        sage: t = OscillatingTableau([[],[1],[1,1],[1],[]])
        sage: t.sundaram()
        ([], [(1, 3), (2, 4)])

        sage: t = OscillatingTableau([[],[1],[1,1],[2,1],[2],[1],[2],[2,1],[2,1,1],[2,1]])
        sage: t.sundaram()
        ([[2, 7], [8]], [(1, 4), (3, 5), (6, 9)])

        sage: s = OscillatingTableau([[],[1],[2],[2,1],[1,1],[1],[1,1],[2,1],[3,1],[2,1]])
        sage: s.sundaram()
        ([[3, 7], [6]], [(1, 5), (2, 4), (8, 9)])

        """
        if self.is_skew():
            raise ValueError("This has only been implemented for straight oscillating tableaux.")
        tb = Tableau([])
        pm = set([])

        for i in range(1,len(self)):
            lb = self[i]
            mu = self[i-1]
            if lb.contains(mu):
                cell = [ c for c in lb.corners() if c not in mu.corners() ][0]
                tb = tb.add_entry(cell,i)
            else:
                cell = [ c for c in mu.corners() if c not in lb.corners() ][0]
                tb, x = tb.reverse_bump(cell)
                pm.add((i,x))

        return (tb,PerfectMatching(pm))

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

