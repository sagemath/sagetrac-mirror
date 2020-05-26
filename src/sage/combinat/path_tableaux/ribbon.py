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
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.list_clone import ClonableArray
from sage.combinat.path_tableaux.path_tableau import PathTableau, PathTableaux
from sage.categories.sets_cat import Sets
from sage.structure.parent import Parent
from sage.combinat.partition import Partition
from sage.combinat.tableau import Tableau
from sage.combinat.skew_tableau import SkewTableau
from sage.combinat.ribbon_tableau import RibbonTableau
from sage.combinat.combinatorial_map import combinatorial_map

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

sage: 

        sage: RibbonPathTableau([],0)
        Traceback (most recent call last):
        ...
        ValueError: 0 is not a valid ribbon size

    """
    @staticmethod
    def __classcall_private__(cls, rt, k=None):
        r"""
        """
        return RibbonPathTableaux()(rt)

    def __init__(self, parent, rt, check=True):
        r"""
        Initialize a ribbon tableau.

        INPUT:

        Can be any of:

        * a list or tuple of partitions
        * a tableau or skew tableau
        * a ribbon tableau

        """
        w = None

        if isinstance(rt, (list,tuple)):
            w = tuple([ Partition(a) for a in rt ])

        if isinstance(rt, (SkewTableau,Tableau)):
            w = tuple([ Partition(a) for a in rt.to_chain() ])

        if w is None:
            raise ValueError(f"invalid input {rt}")

#        if not k:
#            if not w:
#                raise ValueError("k must be defined for the empty sequence")
#            else:
#                k = w[0].size()

        ClonableArray.__init__(self, parent, w, check=check)

    def check(self):
        r""" Checks that ``self`` is a valid path.

        TESTS::

        """
        k = self.parent().k

        for u, v in zip(self,self[1:]):
            if not v.contains(u):
                return False
            if not v.size() == u.size() + k:
                return False
            cu = set(u.cells())
            cv = set(v.cells())
            cb = list({a[0]-a[1] for a in cv if not a in cu})
            cb.sort()
            if not cb == list(range(cb[0],cb[0]+k)):
                return False

        return True

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
            rtx = set(RibbonTableaux([x[2],x[0]],k)) # Do we need weight as well?
            rtx.remove(x[1])
            if rtx:
                return rtx.pop()
            else:
                return x[1]


        if not (i > 0 and i < len(self)-1):
            raise ValueError(f"{i} is not a valid integer")

        with self.clone() as result:
            result[i] = _rule(self[i-1:i+2])

        return result

    @combinatorial_map(name='to ribbon tableau')
    def to_ribbon(self):
        r"""
        Return ''self'' as a ribbon tableau.
        """
        return RibbonTableau(list(self))

class RibbonPathTableaux(PathTableaux):
    """
    The parent class for RibbonPathTableau.
    """

    @staticmethod
    def __classcall_private__(cls, k):
        return super(RibbonPathTableaux, cls).__classcall__(cls, k)

    def __init__(self, k):
        if k<1:
            raise ValueError(f"{k} is not a valid ribbon size")
            
        if not k:
            if not self:
                raise ValueError("k must be defined for the empty sequence")
            elif len(self) == 1:
                k = self[0].size()
            else:
                k = self[1].size() - self[0].size()
        self.k = k
        Parent.__init__(self, category=Sets())

    Element = RibbonPathTableau
