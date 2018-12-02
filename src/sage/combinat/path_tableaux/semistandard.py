r"""
This is an impementation of the Category PathTableaux.

In this implementation we have sequences of partitions. These are in
bijection with dual semistandard tableaux. This gives an effective
version of operations on tableaux constructed using jeu-de-taquin.
In the standard constructions of these operations one usually assumes
the tableau is standard.

For rectification and evacuation the operations here
agree with the standard construction. For promotion the construction
here agrees with the standard construction on rectangular standard
tableaux, but, in general, they are different.

The operations here also give the Bender-Knuth involutions and
dual equivalence graphs.

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
from sage.combinat.path_tableaux.path_tableau import PathTableau, PathTableaux
from sage.combinat.partition import Partition
from sage.combinat.tableau import SemistandardTableau, from_chain
from sage.combinat.skew_tableau import SkewTableau, SemistandardSkewTableaux
from sage.combinat.gelfand_tsetlin_patterns import GelfandTsetlinPattern
from sage.modules.free_module_element import vector
from sage.combinat.combinatorial_map import combinatorial_map

"""
Here we show that we have implemented the jeu-de-taquin operations.
In the Tableau and SkewTableau classes promotion and evacuation have only been implemented
for standard tableaux.

We extend promotion, evacuation and dual equivalence graphs to semistandard tableaux.
We give an efficient implementation of the Bender-Knuth moves and use these for
rectification of a skew semistandard tableaux.

TESTS::

sage: ta = StandardTableau([[1,2,3,4,5],[6,7,8]])
sage: tp = DualSemistandardTableau(ta)
sage: len(tp.orbit())
14
sage: TestSuite(tp).run()


sage: ST = StandardTableaux([4,3,1])
sage: all( DualSemistandardTableau(t.evacuation()) == DualSemistandardTableau(t).evacuation() for t in ST )
True

sage: ST = SemistandardSkewTableaux([[6,5,3],[4,3,1]],max_entry=2)
sage: all(DualSemistandardTableau(t)._test_bender_knuth() for t in ST)
False

"""
@add_metaclass(InheritComparisonClasscallMetaclass)
class DualSemistandardTableau(PathTableau):
    """
       An instance is the sequence of partitions correspond to the
       chain of partitions of a dual semistandard skew tableau.

    The acceptable inputs are:
        - a sequence such that each term defines a partition
        - a semistandard tableau
        - a semistandard skew tableau
        - a Gelfand-Tsetlin pattern

    EXAMPLES:

    sage: DualSemistandardTableau([[],[1],[2],[2,1]])
    [[], [1], [2], [2, 1]]

    sage: t = SkewTableau([[None,None,None,4,4,5,6,7],[None,2,4,6,7,7,7],[None,4,5,8,8,9],[None,6,7,10],[None,8,8,11],[None],[4]])
    sage: DualSemistandardTableau(t)
    [[6, 1, 1], [6, 1, 1], [6, 2, 1], [6, 2, 1], [7, 3, 2, 1, 1], [7, 3, 3, 1, 1, 1], [7, 4, 3, 2, 1, 1, 1], [7, 4, 4, 2, 2, 2, 2, 1], [7, 5, 5, 3, 3, 2, 2, 1], [7, 5, 5, 3, 3, 3, 2, 1], [7, 5, 5, 4, 3, 3, 2, 1], [7, 5, 5, 5, 3, 3, 2, 1]]

    """
    @staticmethod
    def __classcall_private__(self, ot):
        """
            The acceptable inputs are:
        - a sequence such that each term defines a partition
        - a semistandard tableau
        - a semistandard skew tableau
        - a Gelfand-Tsetlin pattern

        EXAMPLES::

        sage: DualSemistandardTableau([[],[1],[2],[2,1]])
        [[], [1], [2], [2, 1]]

        sage: t = SemistandardTableau([[1,1,3,3,3],[2]])
        sage: DualSemistandardTableau(t)
        [[], [1, 1], [2, 1], [2, 1, 1, 1, 1]]
        sage: t = SkewTableau([[None,None,None,4,4,5,6,7],[None,2,4,6,7,7,7],[None,4,5,8,8,9],[None,6,7,10],[None,8,8,11],[None],[4]])
        sage: DualSemistandardTableau(t)
        [[6, 1, 1], [6, 1, 1], [6, 2, 1], [6, 2, 1], [7, 3, 2, 1, 1], [7, 3, 3, 1, 1, 1], [7, 4, 3, 2, 1, 1, 1], [7, 4, 4, 2, 2, 2, 2, 1], [7, 5, 5, 3, 3, 2, 2, 1], [7, 5, 5, 3, 3, 3, 2, 1], [7, 5, 5, 4, 3, 3, 2, 1], [7, 5, 5, 5, 3, 3, 2, 1]]
        sage: gt = GelfandTsetlinPattern([[5,3,1],[4,2],[3]])
        sage: DualSemistandardTableau(gt)
        [[1, 1, 1], [2, 2, 1, 1], [3, 2, 2, 1, 1]]

        """

        if isinstance(ot, DualSemistandardTableau):
            return ot

        w = None

        if isinstance(ot,(list,tuple)):
            try:
                w = tuple([ Partition(a) for a in ot ])
            except TypeError:
                raise ValueError("%s is not a sequence of partitions." % str(ot) )

        if isinstance(ot,(SkewTableau,SemistandardTableau)):
            w = ot.conjugate().to_chain()

        if isinstance(ot,GelfandTsetlinPattern):
            u = list(ot)
            u.reverse()
            v = map(Partition,u)
            w = [t.conjugate() for t in v]

        if w == None:
            raise ValueError( "%s is not a valid input" % str(ot) )

        return DualSemistandardTableaux()(w)

    def check(self):
        n = len(self)
        for i in range(n-1):
            h = self[i]
            t = self[i+1]
            if not t.contains(h):
                raise ValueError( "%s must contain %s" % (str(t),str(h)) )
            for r, s in zip(h,t):
                if s > r+1:
                    raise ValueError( "%s / %s is not a vertical strip" % (str(t),str(h)) )
            for a in t[len(h):]:
                if a > 1:
                    raise ValueError( "%s / %s is not a vertical strip" % (str(t),str(h)) )

    def _local_rule(self,i):
        """
        This has input a list of objects. This method first takes
        the list of objects of length three consisting of the `(i-1)`-st,
        `i`-th and `(i+1)`-term and applies the rule. It then replaces
        the `i`-th object  by the object returned by the rule.

        EXAMPLES::

            sage: t = DualSemistandardTableau([[],[1,1],[2,2,1],[3,2,1]])
            sage: t._local_rule(2)
            [[], [1, 1], [2, 1], [3, 2, 1]]

        """

        def _rule(x):
            y = map(list,x)
            m = max([ len(u) for u in y ])
            z = map( lambda u: vector(u + [0]*(m-len(u)) ), y )
            result = list(z[0]-z[1]+z[2])
            result.sort(reverse=True)
            return Partition(result)

        if not (i > 0 and i < len(self) ):
            raise ValueError("%d is not a valid integer" % i)

        with self.clone() as result:
            result[i] = _rule(self[i-1:i+2])

        return result

    def evaluation(self):
        """
        The evaluation, or content, of a tableau.

        EXAMPLES::

            sage: t = DualSemistandardTableau([[],[1],[2],[2,1]])
            sage: t.evaluation()
            [1, 1, 1]

        """
        z = [ p.size() for p in self ]
        return [ z[i+1] - z[i] for i in range(len(self)-1) ]

    @combinatorial_map(name='to Tableau')
    def to_tableau(self):
        """
        Returns the conjugate skew tableau. This will be semistandard.

        EXAMPLES::

            sage: t = DualSemistandardTableau([[],[1],[2],[2,1]])
            sage: t.to_tableau()
            [[1, 3], [2]]


        """
        ch = [ p.conjugate() for p in self]
        if self.is_skew():
            return SkewTableau(chain=ch)
        else:
            return from_chain(ch)

    @combinatorial_map(name='to Gelfand-Tsetlin pattern')
    def to_gelfand_tsetlin_pattern(self):
        """
        Returns the Gelfand-Tsetlin pattern. This is not implemented
        for skew tableaux as skew Gelfand-Tsetlin patterns are not
        implemented.

        EXAMPLES::

            sage: t = DualSemistandardTableau([[],[1],[2],[2,1]])
            sage: t.to_gelfand_tsetlin_pattern()
            [[2, 1, 0], [1, 1], [1], []]

        """
        if self.is_skew():
            raise ValueError("not implemented for skew tableaux")

        u = [a.conjugate() for a in self]
        s = max(len(u),len(u[0]))
        for i in range(s):
            u[i] = list(u[i]) + [0]*(i-len(u[i]))
        u.reverse()

        return GelfandTsetlinPattern(u)

    def is_skew(self):
        """
        Returns True if Tableau is skew and False if not.

        EXAMPLE::

            sage: DualSemistandardTableau([[],[1],[2],[3,1]]).is_skew()
            False
            sage: DualSemistandardTableau([[1],[2],[3,1]]).is_skew()
            True

        """
        return self[0] != Partition([])

    def rectify(self):
        """
        The rectification of a skew tableau.

        EXAMPLE::

            sage: DualSemistandardTableau([[1],[2],[3,1]]).rectify()
            [[], [1], [2, 1]]

        """
        if not self.is_skew():
            return self

        la = self[0].conjugate()
        mu = [ la[:i] for i in range(len(self)+1) ]
        mu = [ Partition(a).conjugate() for a in mu ]

        return DualSemistandardTableau(mu).commutor(self)[0]

    def _test_bender_knuth(self, **options):
        tester = self._tester(**options)
        for i in range(1,len(self)-1):
            lhs = self._local_rule(i).to_tableau()
            rhs = self.to_tableau().bender_knuth_involution(i)
            tester.assertTrue( lhs == rhs )

    def _test_rectify(self, **options):
        tester = self._tester(**options)

        lhs = self.rectify().to_tableau()
        t = self.to_tableau()
        if isinstance(t,SemistandardTableau):
            rhs == t
        else:
            rhs = SkewTableau(t).rectify()
        tester.assertTrue( lhs == rhs )

    def _test_evacuation(self, **options):
        tester = self._tester(**options)
        if self.is_skew():
            tester.assertTrue(True)
        else:
            lhs = self.evacuation().to_tableau()
            rhs = self.to_tableau().evacuation()
            tester.assertTrue( lhs == rhs )

########################################################################

class DualSemistandardTableaux(PathTableaux):

    Element = DualSemistandardTableau
