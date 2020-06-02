r"""
Vacillating Tableaux
====================

AUTHORS: - Bruce Westbury (2020): initial version

This implements the combinatorics of vacillating tableaux.
This module defines two classes:
    - the Parent class is VacillatingTableaux
    - the Element class is VacillatingTableau

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
from sage.categories.sets_cat import Sets
from sage.combinat.partition import Partition
from sage.combinat.set_partition import SetPartition

###############################################################################

@add_metaclass(InheritComparisonClasscallMetaclass)
class VacillatingTableau(ClonableArray):
    r"""
    A vacillating tableau is a sequence of partitions. A sequence of partitions
    is a valid vacillating tableau if:

    * it has odd length,
    * it starts with the empty partition,
    * every even step either does nothing or adds a box
    * every odd step either does nothing or removes a box

    EXAMPLES::

        sage: VacillatingTableau([[],[],[1],[1],[1],[1],[2],[2],[2,1],[1,1],[1,1]])
        [[], [], [1], [1], [1], [1], [2], [2], [2, 1], [1, 1], [1, 1]]

    TESTS::

        sage: VacillatingTableau([[],[],[1],[1],[1],[1],[2],[2],[2.1],[1,1],[1,1]])
        Traceback (most recent call last):
        ...
        ValueError: all parts of [2.10000000000000] should be nonnegative integers
    """
    @staticmethod
    def __classcall_private__(self, vt):
        """
        This constructs the sequence of partitions and then delegates the
        construction of the vacillating tableau.
        """
        try:
            vt = tuple(Partition(a) for a in vt)
        except TypeError:
            raise ValueError(f"{vt} is not a sequence of partitions.")
        return VacillatingTableaux()(vt)

        raise ValueError(f"{vt} is not valid input")

    def check(self):
        """
        This checks that the sequence of partitions in `self` is a valid
        vacillating tableau.

        TESTS::

            sage: VacillatingTableau([[1]])
            Traceback (most recent call last):
            ...
            ValueError: the first partition [1] must be the empty partition
            sage: VacillatingTableau([[],[1]])
            Traceback (most recent call last):
            ...
            ValueError: the length of [[], [1]] must be odd
            sage: VacillatingTableau([[],[1],[2]])
            Traceback (most recent call last):
            ...
            ValueError: the partition [] must contain the partition [1]
            sage: VacillatingTableau([[],[2],[2]])
            Traceback (most recent call last):
            ...
            ValueError: the partition [] must contain the partition [2]
            sage: VacillatingTableau([[],[1],[]])
            Traceback (most recent call last):
            ...
            ValueError: the partition [] must contain the partition [1]
            sage: VacillatingTableau([[],[],[2]])
            Traceback (most recent call last):
            ...
            ValueError: the even partition [2] is too large
            sage: VacillatingTableau([[],[],[1],[1],[2],[],[]])
            Traceback (most recent call last):
            ...
            ValueError: the odd partition [] is too small
        """
        if self[0] != Partition([]):
            raise ValueError(f"the first partition {self[0]} must be the empty partition")
        if len(self) % 2 == 0:
            raise ValueError(f"the length of {self} must be odd")
        k = len(self) // 2
        # Check the odd steps
        for i in range(k):
            if not self[2*i].contains(self[2*i+1]):
                raise ValueError(f"the partition {self[2*i]} must contain the partition {self[2*i+1]}")
            if self[2*i+1].size() + 1 < self[2*i].size():
                raise ValueError(f"the odd partition {self[2*i+1]} is too small")
        # Check the even steps
        for i in range(k):
            if not self[2*i+2].contains(self[2*i+1]):
                raise ValueError(f"the partition {self[2*i+2]} must contain the partition {self[2*i+1]}")
            if self[2*i+2].size() > self[2*i+1].size() +1:
                raise ValueError(f"the even partition {self[2*i+2]} is too large")

###############################################################################

class VacillatingTableaux(UniqueRepresentation,Parent):
    """
    The parent class for VacillatingTableau
    """
    def __init__(self):
        """
        Initializes the class of all VacillatingTableaux

        TESTS::

            sage: VacillatingTableau([[]]).parent() # indirect test
            <sage.combinat.vacillating.VacillatingTableaux_with_category object at ...>
        """
        Parent.__init__(self, category=Sets())

    def from_setpartition(self,S):
        """
        Construct a vacillating tableau from a set partition.
        """

    Element = VacillatingTableau
