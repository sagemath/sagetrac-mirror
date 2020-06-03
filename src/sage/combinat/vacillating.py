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
from sage.combinat.partition import Partition
from sage.combinat.tableau import SemistandardTableau
from sage.combinat.set_partition import SetPartitions
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.sets.recursively_enumerated_set import RecursivelyEnumeratedSet
from sage.typeset.ascii_art import AsciiArt
from sage.typeset.unicode_art import UnicodeArt

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

    def size(self):
        """
        Return the size of ``self``. This is ``k`` when ``self`` has length ``2k+1``.

        EXAMPLES::

            sage: VacillatingTableau([[],[],[1],[1],[1],[1],[2],[2],[2,1]]).size()
            4

        """
        return len(self) // 2

    def final_shape(self):
        """
        Return the final shape of ``self``.

        EXAMPLES::

            sage: VacillatingTableau([[],[],[1],[1],[1],[1],[2],[2],[2,1]]).final_shape()
            [2, 1]
        """
        return self[-1]

    def crossing_number(self):
        """
        Return the crossing number of ``self``.

        EXAMPLES::

            sage: VacillatingTableau([[],[],[1],[1],[1],[1],[2],[2],[2,1]]).crossing_number()
            2
        """
        if self[-1] == []:
            return 0
        else:
            return max(a[0] for a in self if a)

    def is_noncrossing(self):
        """
        Return True if ``self`` is noncrossing.
        """
        return self.crossing_number() < 2

    def nesting_number(self):
        """
        Return the nesting number of ``self``.

        EXAMPLES::

            sage: VacillatingTableau([[],[],[1],[1],[1],[1],[2],[2],[2,1]]).nesting_number()
            2
        """
        return max(len(a) for a in self)

    def is_nonnesting(self):
        """
        Return True if ``self`` is nonnesting.
        """
        return self.nesting_number() < 2

    def _latex_(self):
        r"""
        Return a LaTeX representation of ``self``.

        EXAMPLES::

            sage: V = VacillatingTableau([[],[],[1]])
            sage: latex(V) # indirect test
        """
        result = "\\left( \\begin{array}{"+'l'*len(self)+"}\n"+latex(self[0])+"\n"
        for p in self[1:]:
            result += "& "+latex(p)+"\n"
        result += "\\end{array} \\right)"
        return result

    def _ascii_art_(self):
        r"""
        Return an ascii art representation of ``self``.

        EXAMPLES::

            sage: V = VacillatingTableau([[],[],[1],[1],[1],[1],[2],[2],[2,1]])
            sage: AsciiArt(V) # indirect test
            [                           ** ]
            [ -, -, *, *, *, *, **, **, *  ]
        """
        return AsciiArt([AsciiArt(p) for p in self])

    def _unicode_art_(self):
        r"""
        Return an ascii art representation of ``self``.

        EXAMPLES::

            sage: V = VacillatingTableau([[],[],[1],[1],[1],[1],[2],[2],[2,1]])
            sage: UnicodeArt(V) # indirect test
            ⎡                                 ┌┬┐ ⎤
            ⎢       ┌┐  ┌┐  ┌┐  ┌┐  ┌┬┐  ┌┬┐  ├┼┘ ⎥
            ⎣ ∅, ∅, └┘, └┘, └┘, └┘, └┴┘, └┴┘, └┘  ⎦
        """
        return UnicodeArt([UnicodeArt(p) for p in self])

    def conjugate(self):
        """
        Return the conjugate of ``self``. This is given by simply conjugating each
        partition in the sequence. This interchanges the crossing and nesting numbers.

        EXAMPLES::

            sage: VacillatingTableau([[],[],[1],[1],[1],[1],[2],[2],[2,1]]).conjugate()
            [[], [], [1], [1], [1], [1], [1, 1], [1, 1], [2, 1]]

        TESTS::


        """
        return VacillatingTableau([a.conjugate() for a in self])

    def stabilise(self,n=None):
        """
        Convert ``self`` to a list of partitions all of size `n`.
        If the parameter `n` is not specified it defaults to the minimal valid
        number.

        EXAMPLES::

            sage: V = VacillatingTableau([[],[],[1],[1],[1],[1],[2],[2],[2,1]])
            sage: V.stabilise()
            [[5], [5], [4, 1], [4, 1], [4, 1], [4, 1], [3, 2], [3, 2], [2, 2, 1]]

        TESTS::

            sage: VacillatingTableau([[],[],[1],[1],[1],[1],[2],[2],[2,1]]).stabilise(4)
            Traceback (most recent call last):
            ...
            ValueError: 4 is too small
        """
        if n == None:
            if all(p == [] for p in self):
                return list(self)
            else:
                n = max(p.size()+p[0] for p in self if p != [])

        if any(n-p.size() < p[0] for p in self if p != []):
            raise ValueError(f"{n} is too small")

        return [ Partition([n-p.size()]+list(p)) for p in self ]

    def to_set_partition(self):
        """
        Return a pair (P,T) where P is a set partition of [n] and T
        is a partial tableau (a semistandard tableau with no repeated entries)
        with content(T) a subset of max(P). The shape of T
        is the final shape of ``self``. In particular, if the final shape of
        ``self`` is the empty partition then this is just a set partition of [n].
        This map between vacillating tableau whose final shape is the empty partition
        and set partitions is a bijection.
        
        EXAMPLES::
            
            sage: V = VacillatingTableau([[],[],[1],[1],[2],[2],[2],[2],[2,1],[2,1],[2,1,1],[2,1],[2,1],[1,1],[2,1]])
            sage: V.to_set_partition()
            ({{1}, {2, 6}, {3}, {4, 7}, {5}}, [[1, 7], [5]])            
        """
        P = set() # empty set
        T = SemistandardTableau([]) # empty tableau
        for i, (lam, mu) in enumerate(zip(self,self[1:])):
            if mu.size() == lam.size()+1:
                assert(mu.contains(lam))
                cell = [c for c in mu.cells() if not c in lam.cells()][0]
                T = T.add_entry(cell,(i+1)//2)
            elif lam.size() == mu.size()+1:
                assert(lam.contains(mu))
                cell = [c for c in lam.cells() if not c in mu.cells()][0]
                T, j = T.reverse_bump(cell)
                P.add((j,i//2+1))
            else:
                assert(lam==mu)

        return SetPartitions().from_arcs(P, self.size()), T

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
        Parent.__init__(self, category=InfiniteEnumeratedSets())

    def __iter__(self):
        """
        An iterator function.

        TESTS::

            sage: VacillatingTableaux()[0]
            {[[]]}

            sage: VacillatingTableaux()[1]
            {[[], [], []], [[], [], [1]]}

            sage: [ len(VacillatingTableaux()[i]) for i in range(5) ]
            [1, 2, 7, 31, 164]
        """
        def succ(tableau):
            """
            Return a list of the successors to the VacillatingTableau ``tableau``.
            """
            p = tableau[-1]
            first_step = [p]+[p.remove_cell(*c) for c in p.removable_cells()]
            second_step = sum([ [[r,r]]+[[r,r.add_cell(*c)] for c in r.addable_cells()] for r in first_step ],[])
            return [VacillatingTableau(list(tableau)+x) for x in second_step]
        R = RecursivelyEnumeratedSet([VacillatingTableau([[]])], succ, structure='graded')
        return R.graded_component_iterator()

    Element = VacillatingTableau
