#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This is an implementation of the Category PathTableaux.
This is the simplest implementation of PathTableaux and is included to
provide a convenient test case and for pedagogical purposes.

In this implementation we have sequences of nonnegative integers. These
are required to be the heights Dyck words (except that we do not require
the sequence to start or end at height zero). These are in bijection
with skew standard tableaux with at most two rows. Sequences which start
and end at height zero are in bijection with noncrossing perfect matchings.

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
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.list_clone import ClonableArray
from sage.structure.parent import Parent

from sage.categories.pathtableaux import PathTableaux
#from sage.categories.sets_cat import Sets
#from sage.combinat.catalan import CatalanTableau
#from sage.combinat.catalan import CatalanTableaux


###############################################################################

@add_metaclass(InheritComparisonClasscallMetaclass)
class CatalanTableau(ClonableArray):
    """
    An instance is the sequence of nonnegative
    integers given by the heights of a Dyck word. The acceptable inputs
    are:
        - a sequence of nonnegative integers
        - a two row standard skew tableau
        - a Dyck word
        - a noncrossing perfect matching

    EXAMPLES:

        sage: CatalanTableau([0,1,2,1,0])
        [0, 1, 2, 1, 0]

        sage: w = DyckWord([1,1,0,0])
        sage: CatalanTableau(w)
        [0, 1, 2, 1, 0]

        sage: p = PerfectMatching([(1,2),(3,4)])
        sage: CatalanTableau(p)
        /home/bruce/sage-8.1/src/bin/sage-ipython:76: DeprecationWarning: is_non_crossing is deprecated. Please use is_noncrossing instead.
        See http://trac.sagemath.org/23982 for details.
        [1, 0, 1, 0]

        sage: t = SkewTableau([[1,2],[3,4]])
        sage: CatalanTableau(t)
        [1, 1, 0, 0]


    """
    @staticmethod
    def __classcall_private__(self, ot):

        w = None

        if isinstance(ot,DyckWord):
            w = ot.heights()

        if isinstance(ot,PerfectMatching):
            if ot.is_non_crossing():
                w = [1]*ot.size()
                for a in ot.arcs():
                    w[a[1]-1] = 0
            else:
                raise ValueError("The perfect matching must be non crossing.")

        if isinstance(ot,SkewTableau):
            if len(ot) == 2:
                if ot.is_standard():
                    w = [1]*ot.size()
                    for i in ot[1]:
                        w[i-1] = 0
                else:
                    raise ValueError("The tableau must be standard.")
            else:
                raise ValueError("The tableau must have two rows.")

        if isinstance(ot,(list,tuple)):
            try:
                w = tuple([ Integer(a) for a in ot ])
            except TypeError:
                raise ValueError("%s is not a sequence of integers." % str(ot) )

        if w == None:
            raise ValueError( "Sorry, not sorry; I don't know what to do with %s." % str(ot) )

        return CatalanTableaux()(w)

    def check(self):
        """
        This overwrites the abstract method.

        This checks that heights are nonnegative and that succesive heights
        differ by +1 or -1.

        EXAMPLES:
        sage: CatalanTableau([0,1,2,3,2,3])
        [0, 1, 2, 3, 2, 3]

        sage: CatalanTableau([0,1,0,-1,0])
        ValueError: [0, 1, 0, -1, 0] has a negative entry.

        sage: CatalanTableau([0,1,3,3,2,3])
        ValueError: [0, 1, 3, 3, 2, 3] is not a Dyck path.

        """
        n = len(self)
        if any([ a < 0 for a in self]):
           raise ValueError( "%s has a negative entry." % (str(self)) )
        for i in range(n-1):
            if abs(self[i+1]-self[i]) > 1:
                raise ValueError( "%s is not a Dyck path." % (str(self)) )

    @staticmethod
    def _rule(x):
        """
        This overwrites the abstract method.
        """
        return abs(x[0]-x[1]+x[2])

    def is_skew(self):
        """
        Returns True if Tableau is skew and False if not.

        EXAMPLES:
        sage: CatalanTableau([0,1,2,1]).is_skew()
        False
        sage: CatalanTableau([1,0,1,2,1]).is_skew()
        True

        """
        return self[0] != 0

    def descents(self):
        """
        Returns the descent set.

        EXAMPLE:

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
        Converts to a word in the alphabet 0,1

        EXAMPLE:

        sage: CatalanTableau([1,0,1,2,1]).to_word()
        [0, 1, 1, 0]

        """
        return [ (self[i+1]-self[i]+1)/2 for i in range(self.size()-1) ]

    def to_dyck_word(self):
        """
        This converts to a Dyck path.

        EXAMPLE:

        sage: CatalanTableau([1,0,1,2,1]).to_dyck_word()
        [0, 1, 1, 0]

        """
        return DyckWord(self.to_word())

    def to_perfect_matching(self):
        """
        This converts to a perfect matching.

        EXAMPLE:

        sage: CatalanTableau([0,1,2,1,2,1,0,1,0]).to_perfect_matching()
        [(0, 5), (1, 2), (3, 4), (6, 7)]

        """
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
        Converts to a skew tableau.
        """
        top = [ i for i, a in enumerate(self) if a == 1 ]
        bot = [ i for i, a in enumerate(self) if a == 0 ]
        return SkewTableau([[None]*self[0]+top,bot])

    def draw(self):
        """
        This draws the Dyck path.
        """
        return line([ (i,a) for i, a in enumerate(self)])

"""

Here we illustrate the slogan that promotion = rotation.

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

Here we test the Category methods.

sage: t = CatalanTableau([0,1,2,3,2,1,0])
sage: SkewTableau(t.cylindrical_diagram()).pp()
  0  1  2  3  2  1  0
  .  0  1  2  1  0  1  0
  .  .  0  1  0  1  2  1  0
  .  .  .  0  1  2  3  2  1  0
  .  .  .  .  0  1  2  1  0  1  0
  .  .  .  .  .  0  1  0  1  2  1  0
  .  .  .  .  .  .  0  1  2  3  2  1  0



"""
###############################################################################

class CatalanTableaux(UniqueRepresentation,Parent):
    """
    This constructs the Parent class.
    """
    @staticmethod
    def __classcall_private__(cls):
        return super(CatalanTableaux, cls).__classcall__(cls)

    def __init__(self):

        Parent.__init__(self, category=PathTableaux())

    def __contains__(self, ot):

        return isinstance(ot, (list, tuple, CatalanTableau))

    def _element_constructor_(self, ot, check=True):

        if isinstance(ot, CatalanTableaux) and ot.parent() == self:
            return ot

        return self.element_class(self, list(ot))

    Element = CatalanTableau

