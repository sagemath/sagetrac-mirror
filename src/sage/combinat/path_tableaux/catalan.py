
#*****************************************************************************
#       Copyright (C) 2018 Bruce Westbury <bruce.westbury@gmail.com>,
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.list_clone import ClonableArray
from sage.combinat.path_tableaux.path_tableau import PathTableau, PathTableaux, CylindricalDiagram
from sage.combinat.combinatorial_map import combinatorial_map
from sage.combinat.dyck_word import DyckWord
from sage.combinat.perfect_matching import PerfectMatching
from sage.combinat.skew_tableau import SkewTableau
from sage.combinat.tableau import Tableau, StandardTableau
from sage.rings.integer import Integer

###############################################################################
r"""
Catalan Tableaux

This is an implementation of the abstract base class
:class:`sage.combinat.pathtableau.pathtableaux`.
This is the simplest implementation of PathTableaux and is included to
provide a convenient test case and for pedagogical purposes.

In this implementation we have sequences of nonnegative integers. These
are required to be the heights Dyck words (except that we do not require
the sequence to start or end at height zero). These are in bijection
with skew standard tableaux with at most two rows. Sequences which start
and end at height zero are in bijection with noncrossing perfect matchings.

AUTHORS:

- Bruce Westbury (2018): initial version

Here we illustrate the slogan that promotion = rotation.

EXAMPLES::

    sage: t = CatalanTableau([0,1,2,3,2,1,0])
    sage: t.to_perfect_matching()
    [(0, 5), (1, 4), (2, 3)]

    sage: t = t.promotion()
    sage: t.to_perfect_matching()
    [(0, 3), (1, 2), (4, 5)]

    sage: t = t.promotion()
    sage: t.to_perfect_matching()
    [(0, 1), (2, 5), (3, 4)]

    sage: t = t.promotion()
    sage: t.to_perfect_matching()
    [(0, 5), (1, 4), (2, 3)]

EXAMPLES::

    sage: t = CatalanTableau([0,1,2,3,2,1,0])
    sage: print(CylindricalDiagram(t))
    The cylindrical growth diagram:
    [0, 1, 2, 3, 2, 1, 0]
    ['', 0, 1, 2, 1, 0, 1, 0]
    ['', '', 0, 1, 0, 1, 2, 1, 0]
    ['', '', '', 0, 1, 2, 3, 2, 1, 0]
    ['', '', '', '', 0, 1, 2, 1, 0, 1, 0]
    ['', '', '', '', '', 0, 1, 0, 1, 2, 1, 0]
    ['', '', '', '', '', '', 0, 1, 2, 3, 2, 1, 0]
    sage: TestSuite(t).run()
"""

class CatalanTableau(ClonableArray,PathTableau):
    """
    An instance is the sequence of nonnegative
    integers given by the heights of a Dyck word.



    INPUT:

        - a sequence of nonnegative integers
        - a two row standard skew tableau
        - a Dyck word
        - a noncrossing perfect matching

    EXAMPLES::

        sage: CatalanTableau([0,1,2,1,0])
        [0, 1, 2, 1, 0]

        sage: w = DyckWord([1,1,0,0])
        sage: CatalanTableau(w)
        [0, 1, 2, 1, 0]

        sage: p = PerfectMatching([(1,2),(3,4)])
        sage: CatalanTableau(p)
        [0, 1, 0, 1, 0]

        sage: t = Tableau([[1,2],[3,4]])
        sage: CatalanTableau(t)
        [0, 1, 2, 1, 0]
    """

    @staticmethod
    def __classcall_private__(cls, ot):
        """This is the preprocessing for creating paths.

        """
        w = None

        if isinstance(ot, DyckWord):
            w = ot.heights()

        if isinstance(ot, PerfectMatching):
            if ot.is_noncrossing():
                u = [1]*ot.size()
                for a in ot.arcs():
                    u[a[1]-1] = 0
                w = DyckWord(u).heights()
            else:
                raise ValueError("the perfect matching must be non crossing")

        if isinstance(ot, Tableau):
            if len(ot) <= 2:
                if ot.is_standard():
                    u = [1] * ot.size()
                    for i in ot[1]:
                        u[i-1] = 0
                    w = DyckWord(u).heights()
                else:
                    raise ValueError("the tableau must be standard")
            else:
                raise ValueError("the tableau must have at most two rows")

        if isinstance(ot, SkewTableau):
            if len(ot) <= 2:
                # The check that ot is standard is not implemented
                u = [1] * ot.size()
                for i in ot[1]:
                    if i is not None:
                        u[i-1] = 0
                w = DyckWord(u).heights()
            else:
                raise ValueError("the skew tableau must have at most two rows")

        if isinstance(ot, (list,tuple)):
            try:
                w = tuple([Integer(a) for a in ot])
            except TypeError:
                raise ValueError("%s is not a sequence of integers" % ot)

        if w is None:
            raise ValueError("invalid input %s" % ot)

        return CatalanTableaux()(w)

    def check(self):

        n = len(self)
        if any(a < 0 for a in self):
           raise ValueError( "%s has a negative entry" % (str(self)) )
        for i in range(n-1):
            if abs(self[i+1]-self[i]) != 1:
                raise ValueError( "%s is not a Dyck path" % str(self) )

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
            This is the rule on a sequence of three letters.
            """
            return abs(x[0]-x[1]+x[2])

        if not (i > 0 and i < len(self) ):
            raise ValueError("%d is not a valid integer" % i)

        with self.clone() as result:
            result[i] = _rule(self[i-1:i+2])

        return result

    def is_skew(self):
        """
        Return ``True`` if ``self`` is skew and ``False`` if not.

        EXAMPLES::

            sage: CatalanTableau([0,1,2,1]).is_skew()
            False

            sage: CatalanTableau([1,0,1,2,1]).is_skew()
            True
        """
        return self[0] != 0

    @combinatorial_map(name='to Dyck word')
    def to_DyckWord(self):
        """
        Converts ``self`` to a Dyck word.

        EXAMPLES::

            sage: c = CatalanTableau([0,1,2,1,0])
            sage: c.to_DyckWord()
            [1, 1, 0, 0]
        """
        return DyckWord(heights_sequence = list(self))

    def descents(self):
        """
        Return the descent set of ``self``.

        EXAMPLES::

            sage: CatalanTableau([0,1,2,1,2,1,0,1,0]).descents()
            {3, 6}
        """
        result = set()

        for i in range(1,len(self)-1):
            if self[i] < self[i-1] and self[i] < self[i+1]:
                result.add(i)

        return result

    def to_word(self):
        """
        Return the word in the alphabet `\{0,1\}` associated to ``self``.

        EXAMPLES::

            sage: CatalanTableau([1,0,1,2,1]).to_word()
            [0, 1, 1, 0]
        """
        return [ (self[i+1]-self[i]+1)/2 for i in range(self.size()-1) ]

    def to_perfect_matching(self):
        """
        Return the perfect matching associated to ``self``.

        EXAMPLES::

            sage: CatalanTableau([0,1,2,1,2,1,0,1,0]).to_perfect_matching()
            [(0, 5), (1, 2), (3, 4), (6, 7)]
        """
        if self.is_skew():
            raise ValueError( "%s does not start at 0" % (str(self)) )
        w = self.to_word()
        y = DyckWord(w)
        pairs = set()
        for i, a in enumerate(y):
            c = y.associated_parenthesis(i)
            if i < c:
                pairs.add((i,c))
        return PerfectMatching(pairs)

    def to_tableau(self):
        """
        Return the skew tableau associated to ``self``.

        EXAMPLES::

            sage: T = CatalanTableau([0,1,2,3,2,3])
            sage: T.to_tableau()
            [[1, 2, 3, 5], [4]]

            sage: U = CatalanTableau([2,3,2,3])
            sage: U.to_tableau()
            [[None, None, 1, 3], [2]]
        """
        w = self.to_word()
        top = [ i+1 for i, a in enumerate(w) if a == 1 ]
        bot = [ i+1 for i, a in enumerate(w) if a == 0 ]
        if self.is_skew():
            return SkewTableau([[None]*self[0]+top,bot])
        else:
            return StandardTableau([top,bot])

class CatalanTableaux(PathTableaux):
    """
    The parent class for CatalanTableau
    """
    Element = CatalanTableau

