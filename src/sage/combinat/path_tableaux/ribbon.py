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
#from sage.combinat.combinatorial_map import combinatorial_map

###############################################################################

class RibbonTableau(PathTableau):
    r"""
    An instance is a sequence of partitions representing a ribbon tableau.

    INPUT:
        
    """
    @staticmethod
    def __classcall_private__(cls, rt):
        r"""
        """
        return RibbonTableaux()(rt)

    def __init__(self, parent, rt, check=True):
        r"""
        Initialize a Catalan tableau.

        INPUT:

        Can be any of:
        """
        
        ClonableArray.__init__(self, parent, w, check=check)

    def check(self):
        r""" Checks that ``self`` is a valid path.

        TESTS::
            
        """
        
    def _local_rule(self,i):
        """
        This has input a list of objects. This method first takes
        the list of objects of length three consisting of the `(i-1)`-st,
        `i`-th and `(i+1)`-term and applies the rule. It then replaces
        the `i`-th object  by the object returned by the rule.

        EXAMPLES::

        TESTS::

        """

        def _rule(x):
            """
            This is the rule on a sequence of three letters.
            """
            return abs(x[0]-x[1]+x[2])

        if not (i > 0 and i < len(self)-1):
            raise ValueError(f"{i} is not a valid integer")

        with self.clone() as result:
            result[i] = _rule(self[i-1:i+2])

        return result

class RibbonTableaux(PathTableaux):
    """
    The parent class for RibbonTableau.
    """
    Element = RibbonTableau