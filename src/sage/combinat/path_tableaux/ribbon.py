r"""
Ribbon Tableaux

This is an implementation of the abstract base class
:class:`sage.combinat.pathtableau.pathtableaux`.

AUTHORS:

- Bruce Westbury (2020): initial version

EXAMPLES::

    """

#*****************************************************************************
#       Copyright (C) 2018 Bruce Westbury <bruce.westbury@gmail.com>,
#
# This program is free software: you can redistribute it and/or modifyde
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.combinat.path_tableaux.path_tableau import PathTableau, PathTableaux
from sage.categories.sets_cat import Sets
from sage.structure.parent import Parent
from sage.structure.list_clone import ClonableArray
from sage.combinat.partition import Partition
from sage.combinat.tableau import Tableau
from sage.combinat.skew_tableau import SkewTableau, StandardSkewTableaux, SemistandardSkewTableaux
from sage.combinat.ribbon_tableau import RibbonTableau
from sage.combinat.combinatorial_map import combinatorial_map
from contextlib import suppress

###############################################################################

class RibbonPathTableau(PathTableau):
    r"""
    An instance is a sequence of partitions representing a ribbon tableau.

    INPUT:

    EXAMPLES::

        sage: RibbonPathTableau([],3)
        []
        sage: RibbonPathTableau([[],[3]],3)
        [[], [3]]

    TESTS::


        sage: RibbonPathTableau([[],[3]],2)
        Traceback (most recent call last):
        ...
        ValueError: partition sizes must increase by 2 at each step

        sage: RibbonPathTableau([],2)
        []

        sage: RibbonPathTableau([[1]],2)
        [[1]]

        sage: RibbonPathTableau([[1]])
        Traceback (most recent call last):
        ...
        ValueError: the ribbon size k must be specified for [[1]]

        sage: RibbonPathTableau([],0)
        Traceback (most recent call last):
        ...
        ValueError: 0 is not a valid ribbon size

        sage: RibbonPathTableau(Tableau([[1, 2, 3], [4]]),2)
        Traceback (most recent call last):
        ...
        ValueError: partition sizes must increase by 2 at each step

        sage: RibbonPathTableau([1/2],2)
        Traceback (most recent call last):
        ...
        ValueError: all parts of [1/2] should be nonnegative integers

        sage: RibbonPathTableau([],1/2)
        Traceback (most recent call last):
        ...
        ValueError: 1/2 is not a valid ribbon size

        sage: [ RibbonPathTableau(T,1) for T in StandardTableaux(4) ]
        [[[], [1], [2], [3], [4]],
        ...
        [[], [1], [1, 1], [1, 1, 1], [1, 1, 1, 1]]]

        sage: [ RibbonPathTableau(T,1) for T in StandardSkewTableaux([[2,2,1],[1]]) ]
        [[[1], [1, 1], [2, 1], [2, 1, 1], [2, 2, 1]],
         [[1], [1, 1], [1, 1, 1], [2, 1, 1], [2, 2, 1]],
         [[1], [2], [2, 1], [2, 1, 1], [2, 2, 1]],
         [[1], [2], [2, 1], [2, 2], [2, 2, 1]],
         [[1], [1, 1], [2, 1], [2, 2], [2, 2, 1]]]

        sage: RibbonPathTableau(Partition([2]),2)
        Traceback (most recent call last):
        ...
        ValueError: invalid input [2]

    """
    @staticmethod
    def __classcall_private__(cls, rt, k=None):
        r"""
        """
        if k == None:
            if len(rt) == 0:
                raise ValueError(f"the ribbon size k must be specified for {rt}")

            if isinstance(rt, (list,tuple)):
                if len(rt) == 1:
                    raise ValueError(f"the ribbon size k must be specified for {rt}")
                k = Partition(rt[1]).size() - Partition(rt[0]).size

            if isinstance(rt, (SkewTableau,Tableau)):
                k = rt.weight()[0]

        if k < 1:
            raise ValueError(f"{k} is not a valid ribbon size")

        return RibbonPathTableaux(k)(rt)

    def __init__(self, parent, rt, check=True):
        r"""
        Construct an element of ``self`` from ``rt``.
        """
        if isinstance(rt, RibbonPathTableau) and rt.parent() == self:
            return rt

        w = None

        if isinstance(rt, (list,tuple)):
            w = tuple([ Partition(a) for a in rt ])

        if isinstance(rt, (SkewTableau,Tableau)):
            w = tuple([ Partition(a) for a in rt.to_chain() ])

        if w is None:
            raise ValueError(f"invalid input {rt}")

        ClonableArray.__init__(self, parent, w, check=check)

    def check(self):
        r""" Checks that ``self`` is a valid¨ ribbon tableau.

        TESTS::

            sage: RibbonPathTableau([],0) # indirect test
            Traceback (most recent call last):
            ...
            ValueError: 0 is not a valid ribbon size

        """
        k = self.parent().k
        if k and any(v.size() != u.size()+k for u,v in zip(self,self[1:])):
            raise ValueError(f"partition sizes must increase by {k} at each step")

        for u, v in zip(self,self[1:]):
            if not v.contains(u):
                raise ValueError(f"partition {v} does not contain {u}")
            cu = set(u.cells())
            cv = set(v.cells())
            cb = list({a[0]-a[1] for a in cv if not a in cu})
            cb.sort()
            if not cb == list(range(cb[0],cb[-1]+1)):
                raise ValueError(f"the skew shape {v}\{u} is not a ribbon shape")

    def _local_rule(self,i):
        """
        This has input a list of objects. This method first takes
        the list of objects of length three consisting of the `(i-1)`-st,
        `i`-th and `(i+1)`-term and applies the rule. It then replaces
        the `i`-th object  by the object returned by the rule.

        EXAMPLES::

        TESTS::

        Test this agrees with Darij's Bender-Knuth involution.

        """

        def _rule(x):
            """
            This is the rule on a sequence of three partitions.
            """
            k = self.parent().k
            stx = list(SemistandardSkewTableaux([x[2],x[0]],[k,k]))
            if len(stx) == 1:
                return x[1]
            rtx = []
            for T in stx:
                with suppress(ValueError):
                    RibbonPathTableau(T)
                rtx.append(T)
            if len(rtx) == 1:
                return x[1]
            else:
                rtx.remove(x)
                return rtx[0][1]

        if not (i > 0 and i < len(self)-1):
            raise ValueError(f"{i} is not a valid integer")

        with self.clone() as result:
            result[i] = _rule(self[i-1:i+2])

        return result

    @combinatorial_map(name='to ribbon tableau')
    def to_ribbon(self):
        r"""
        Return ''self'' as a ribbon tableau.

        EXAMPLES::

            sage: T = Tableau([[1, 1, 2, 2], [3, 3, 5, 5], [4, 4, 6, 6]])
            sage: RibbonPathTableau(T).to_ribbon()
            [[1, 1, 2, 2], [3, 3, 5, 5], [4, 4, 6, 6]]

        """
        from sage.combinat.tableau import from_chain
        return RibbonTableau(from_chain(self))

class RibbonPathTableaux(PathTableaux):
    """
    The parent class for RibbonPathTableau.
    """

    @staticmethod
    def __classcall_private__(cls, k=None):
        if k<1:
            raise ValueError(f"{k} is not a valid ribbon size")

        return super(RibbonPathTableaux, cls).__classcall__(cls, k)

    def __init__(self, k):

        self.k = k
        Parent.__init__(self, category=Sets())

    def __contains__(self,rt):
        r"""
        Return ''True'' if ''rt'' is in ''self''.
        """

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: RibbonPathTableaux(4)
            4-Ribbon Tableaux

        """
        return f"{self.k}-Ribbon Tableaux"

    Element = RibbonPathTableau
