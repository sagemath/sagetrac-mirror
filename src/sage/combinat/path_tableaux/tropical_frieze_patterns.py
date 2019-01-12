r"""
Tropical Frieze Patterns

This is an implementation of the abstract base class
:class:`sage.combinat.pathtableau.pathtableaux`.

In this implementation we have sequences of nonnegative integers.
This follows [Pec2012]_

REFERENCES:

.. [Pec2012] Oliver Pechenik.
   *Cyclic sieving of increasing tableaux and small Schroder paths*,
   :arxiv:`1209.1355`
   
AUTHORS:

- Bruce Westbury (2019): initial version
"""
#*****************************************************************************
#       Copyright (C) 2019 Bruce Westbury <bruce.westbury@gmail.com>,
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
from sage.combinat.path_tableaux.path_tableau import PathTableau, PathTableaux
from sage.combinat.combinatorial_map import combinatorial_map
from sage.combinat.tableau import Tableau, StandardTableau
from sage.rings.integer import Integer

###############################################################################

"""

EXAMPLES::

    sage: t = TropicalFriezePattern([0,1,2,3,2,1,0])
    sage: t.orbit()
    {[0, -1, 0, -1, 0, 1, 0],
     [0, -1, 0, 1, 0, 1, 0],
     [0, -1, 0, 1, 2, 1, 0],
     [0, 1, 0, -1, 0, 1, 0],
     [0, 1, 0, 1, 0, 1, 0],
     [0, 1, 0, 1, 2, 1, 0],
     [0, 1, 2, 1, 0, 1, 0],
     [0, 1, 2, 1, 2, 1, 0],
     [0, 1, 2, 3, 2, 1, 0]}

    sage: SkewTableau(t.cylindrical_diagram()).pp()
      0  1  2  3  2  1  0
      .  0  1  2  1  0 -1  0
      .  .  0  1  0 -1  0  1  0
      .  .  .  0 -1  0  1  2  1  0
      .  .  .  .  0  1  2  3  2  1  0
      .  .  .  .  .  0  1  2  1  0 -1  0
      .  .  .  .  .  .  0  1  0 -1  0  1  0

    sage: TestSuite(t).run()

    sage: t = TropicalFriezePattern([0,1,2,1,1,0])
    sage: t.orbit()
    {[0, -1, 0, 0, 1, 0],
     [0, -1, 0, 1, 1, 0],
     [0, 0, -1, 0, 1, 0],
     [0, 0, 1, 0, 1, 0],
     [0, 0, 1, 2, 1, 0],
     [0, 1, 0, 0, 1, 0],
     [0, 1, 0, 1, 1, 0],
     [0, 1, 1, 0, 1, 0],
     [0, 1, 1, 2, 1, 0],
     [0, 1, 2, 1, 1, 0],
     [0, 1, 2, 2, 1, 0]}

    sage: SkewTableau(t.cylindrical_diagram()).pp()
      0  1  2  1  1  0
      .  0  1  0  0 -1  0
      .  .  0 -1  0  0  1  0
      .  .  .  0  1  1  2  1  0
      .  .  .  .  0  0  1  0 -1  0
      .  .  .  .  .  0  1  0  0  1  0

    sage: TestSuite(t).run()
    Failure in _test_promotion:
    Traceback (most recent call last):
        ...
    AssertionError: False is not true
    ------------------------------------------------------------
    The following tests failed: _test_promotion

"""


@add_metaclass(InheritComparisonClasscallMetaclass)
class TropicalFriezePattern(ClonableArray,PathTableau):
    """
    An instance is the sequence of nonnegative
    integers.
    """

    @staticmethod
    def __classcall_private__(cls, fp):
        """This is the preprocessing for creating paths.

        INPUT:

            - a sequence of nonnegative integers

        EXAMPLES::

            sage: TropicalFriezePattern([0,1,2,1,1,0])
            [0, 1, 2, 1, 1, 0]

        """
        w = None

        if isinstance(fp, (list,tuple)):
            try:
                w = tuple([Integer(a) for a in fp])
            except TypeError:
                raise ValueError("%s is not a sequence of integers" % fp)

        if w is None:
            raise ValueError("invalid input %s" % fp)

        return TropicalFriezePatterns()(w)

    def check(self):

        n = len(self)
        if any(a < -1 for a in self):
           raise ValueError( "%s has an entry below -1" % (str(self)) )
        for i in range(n-1):
            if abs(self[i+1]-self[i]) > 1:
                raise ValueError( "successive terms differ by more than 1" )
            if self[i] == -1 and self[i+1] == -1:
                raise ValueError( "there are successive -1s" )

    def _local_rule(self,i):
        """
        This has input a list of objects. This method first takes
        the list of objects of length three consisting of the `(i-1)`-st,
        `i`-th and `(i+1)`-term and applies the rule. It then replaces
        the `i`-th object  by the object returned by the rule.

        EXAMPLES::

            sage: t = TropicalFriezePattern([0,1,2,1,1,0])
            sage: t._local_rule(3)
            [0, 1, 2, 2, 1, 0]
        """

        def _rule(x):
            """
            This is the rule on a sequence of three letters.
            """
            return max(x[0]+x[2],0)-x[1]

        if not (i > 0 and i < len(self) ):
            raise ValueError("%d is not a valid integer" % i)

        with self.clone() as result:
            result[i] = _rule(self[i-1:i+2])

        return result

    def is_skew(self):
        """
        Return ``True`` if ``self`` is skew and ``False`` if not.

        EXAMPLES::

            sage: TropicalFriezePattern([0,1,2,1]).is_skew()
            False

            sage: TropicalFriezePattern([1,0,1,2,1]).is_skew()
            True
        """
        return self[0] != 0

    @combinatorial_map(name='to Schroder path')
    def to_height_word(self):
        """
        Converts ``self`` to the height sequence of a small Schroder path

        EXAMPLES::

            sage: TropicalFriezePattern([0,1,2,1,1,0]).to_height_word()
            [0, 1, 2, 3, 2, 2, 1, 0]
        """
        return [0]+[ a+1 for a in self ]+[0]

    @combinatorial_map(name='to increasing tableau')
    def to_increasing_tableau(self):
        """
        Return the increasing tableau associated to ``self``.

        EXAMPLES::

           sage: TropicalFriezePattern([0,1,2,1,1,0]).to_increasing_tableau()
           [[1, 2, 3, 5], [4, 5, 6, 7]]

        """
        if self.is_skew():
            raise ValueError("only implemented for straight shapes")

        h = self.to_height_word()
        top = []
        bot = []
        for i in range(1,len(h)):
            if h[i] == h[i-1]+1:
                top = top + [i]
            elif h[i] == h[i-1]:
                top = top + [i]
                bot = bot + [i]
            elif h[i] == h[i-1]-1:
                bot = bot + [i]
            else:
                raise RuntimeError("this can't happen")

        return Tableau([top,bot])

    @combinatorial_map(name='to flag tableau')
    def to_flag_tableau(self):
        """
        Return the increasing tableau associated to ``self``.

        EXAMPLES::

           sage: TropicalFriezePattern([0,1,2,1,1,0]).to_flag_tableau()
           [[1, 2, 3], [4, 6, 7], [5]]
        """
        if self.is_skew():
            raise ValueError("only implemented for straight shapes")

        h = self.to_height_word()
        result = [[],[]]
        for i in range(1,len(h)):
            if h[i] == h[i-1]+1:
                result[0] = result[0] + [i]
            elif h[i] == h[i-1]:
                result = result + [[i]]
            elif h[i] == h[i-1]-1:
                result[1] = result[1] + [i]
            else:
                raise RuntimeError("this can't happen")

        return StandardTableau(result)

class TropicalFriezePatterns(PathTableaux):

    Element = TropicalFriezePattern

