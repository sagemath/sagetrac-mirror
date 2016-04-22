# -*- coding: utf-8 -*-
r"""
Tabloid

AUTHORS:

- Jackson Criswell (2016): initial version



This file consists of the following major classes:

Element classes:

* :class:`Tableau`
* :class:`SemistandardTableau`
* :class:`StandardTableau`

Factory classes:

* :class:`Tableaux`
* :class:`SemistandardTableaux`
* :class:`StandardTableaux`

Parent classes:

* :class:`Tableaux_all`
* :class:`Tableaux_size`
* :class:`SemistandardTableaux_all` (facade class)
* :class:`SemistandardTableaux_size`
* :class:`SemistandardTableaux_size_inf`
* :class:`SemistandardTableaux_size_weight`
* :class:`SemistandardTableaux_shape`
* :class:`SemistandardTableaux_shape_inf`
* :class:`SemistandardTableaux_shape_weight`
* :class:`StandardTableaux_all` (facade class)
* :class:`StandardTableaux_size`
* :class:`StandardTableaux_shape`

For display options, see :meth:`Tableaux.global_options`.

.. TODO:

    
"""

#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
#                     2011 Jason Bandlow <jbandlow@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.sets.family import Family
from sage.sets.non_negative_integers import NonNegativeIntegers
from sage.structure.global_options import GlobalOptions
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.list_clone import ClonableList
from sage.structure.parent import Parent
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.rings.infinity import PlusInfinity
from sage.arith.all import factorial, binomial
from sage.rings.integer import Integer
from sage.combinat.composition import Composition, Compositions
from integer_vector import IntegerVectors
import sage.libs.symmetrica.all as symmetrica
import sage.misc.prandom as random
import permutation
import itertools
from sage.groups.perm_gps.permgroup import PermutationGroup
from sage.misc.all import uniq, prod
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.sets_cat import Sets
from sage.combinat.combinatorial_map import combinatorial_map
##JACKSONSTART
from sage.sets.set import Set
from sage.combinat.permutation import Permutation
from sage.combinat.permutation import Permutations 
from sage.combinat.combinatorial_algebra import CombinatorialFreeModule
from sage.rings.rational_field import QQ

##JACKSONEND

from sage.combinat.tableau import *
##############
class Tabloid2(Tableau):

    def __init__(self,parent,t):
        tabloid=[]
        for row in t:
            tabloid+=[[Set(row)]]

        ClonableList.__init__(self, parent,t)


    def size(self):
        """
        Return the size of the tableaux contained in the tabloid.



        """""

        return sum(shape(self))

    def shape(self):
        """
        Return the shape of the tableaux contained in the tabloid.



        """
        from sage.combinat.partition import Partition
        shape=[]
        for row in self:
            shape+=len(row[0][0])
        return Partition(shape)

class Tabloid(Tableau):
    """
    Class to model a tabloid of specified shape, as described in Sagan (2000).

   

    INPUT:

    - ``t`` -- a tableau, a list of iterables, or an empty list

    - ``shape`` -- a partition, or list representation of a partition

    OUTPUT:

    - A Tabloid object constructed from ``t``.

    Tabloids can be consider as 1 column generalized tableau, whose entries
    are disjoint sets of integers in [n], in weakly decreasing order by size

    EXAMPLES::



    .. SEEALSO:

        - :class:`Tableaux`
        - :class:`Tableau`
        - :class:`SemistandardTableaux`
        - :class:`StandardTableaux`
        - :class:`StandardTableau`

    TESTS::


    """
    @staticmethod
    def __classcall_private__(self, t):
        r"""
        This ensures that a SemistandardTableau is only ever constructed as an
        element_class call of an appropriate parent.

        TESTS::

            sage: t = SemistandardTableau([[1,1],[2]])
            sage: TestSuite(t).run()

            sage: t.parent()
            Semistandard tableaux
            sage: t.category()
            Category of elements of Semistandard tableaux
            sage: type(t)
            <class 'sage.combinat.tableau.SemistandardTableaux_all_with_category.element_class'>
        """
        if isinstance(t, Tabloid):
            return t
        elif t in Tabloids():
            return Tabloids_all().element_class(Tabloids_all(), t)
        tabloid=[]
        for row in t:
            tabloid+=[[Set(row)]]
        return Tableau(tabloid)
        # t is not a semistandard tableau so we give an appropriate error message
        if t not in Tableaux():
            raise ValueError('%s is not a tableau' % t)

        if not all(isinstance(c,(int,Integer)) and c>0 for row in t for c in row):
            raise ValueError("entries must be positive integers"%t)


    def __init__(self, parent, t):
        r"""
        Initialize a tabloid.

        TESTS::

            sage: t = Tableaux()([[1,1],[2]])
            sage: s = SemistandardTableaux(3)([[1,1],[2]])
            sage: s==t
            True
            sage: s.parent()
            Semistandard tableaux of size 3 and maximum entry 3
            sage: r = SemistandardTableaux(3)(t); r.parent()
            Semistandard tableaux of size 3 and maximum entry 3
            sage: isinstance(r, Tableau)
            True
            sage: s2 = SemistandardTableaux(3)([(1,1),(2,)])
            sage: s2 == s
            True
            sage: s2.parent()
            Semistandard tableaux of size 3 and maximum entry 3
        """
        tabloid=[]
        for row in t:
            tabloid+=[[Set(row)]]
        super(Tabloid, self).__init__(parent, tabloid)

        # Tableau() has checked that t is tableau, so it remains to check that
        # the entries of t are positive integers which are weakly increasing
        # along rows
        from sage.sets.positive_integers import PositiveIntegers
        PI = PositiveIntegers()

        for row in t:
            if any(c not in PI for c in row):
                raise ValueError("the entries of a tabloid must be non-negative integers")
 


    def size(self):
        """
        Return the size of the tableaux contained in the tabloid.



        """""

        return sum(shape(self))

    def shape(self):
        """
        Return the shape of the tableaux contained in the tabloid.



        """
        from sage.combinat.partition import Partition
        shape=[]
        for row in self:
            shape+=len(row[0][0])
        return Partition(shape)

