# -*- coding: utf-8 -*-
r"""
Semistandard Set Valued Tableaux

AUTHORS:

- Jeremy Meza, Oliver Pechenik, Wencin Poh (2019): initial version
"""

#*****************************************************************************
#       Copyright (C) 
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import print_function, absolute_import
from six.moves import range, zip, map
from six import add_metaclass

from sage.combinat.partition import Partition, Partitions, _Partitions, Partitions_n
from sage.combinat.tableau import Tableau, Tableaux, SemistandardTableaux, SemistandardTableau
from sage.combinat.integer_vector import IntegerVectors
from sage.combinat.composition import Composition, Compositions
#from sage.combinat.integer_vector import IntegerVectors, integer_vectors_nk_fast_iter
from sage.combinat.combinatorial_map import combinatorial_map
#from sage.combinat.posets.posets import Poset
#from sage.combinat import permutation

from sage.categories.regular_crystals import RegularCrystals
from sage.categories.classical_crystals import ClassicalCrystals
from sage.combinat.root_system.cartan_type import CartanType

from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.sets.family import Family
#from sage.sets.non_negative_integers import NonNegativeIntegers

from sage.categories.enumerated_sets import EnumeratedSets
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.sets_cat import Sets

#from sage.structure.list_clone import ClonableArray
#from sage.structure.sage_object import SageObject
from sage.structure.element_wrapper import ElementWrapper
from sage.structure.global_options import GlobalOptions
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.list_clone import ClonableList
from sage.structure.parent import Parent
from sage.structure.richcmp import richcmp, richcmp_method
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass

#from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
from sage.rings.infinity import PlusInfinity
from sage.rings.integer import Integer

#from sage.arith.all import factorial, binomial
#import sage.libs.symmetrica.all as symmetrica
#from sage.groups.perm_gps.permgroup import PermutationGroup

#import sage.misc.prandom as random
#from sage.misc.all import prod
from sage.misc.misc import powerset
from sage.misc.lazy_attribute import lazy_attribute

#####################################
# Semistandard Set-Valued Tableaux  #
#####################################

class SemistandardSetValuedTableau(Tableau):
    """
    A class to model a semistandard set-valued tableau.
    """
    @staticmethod
    def __classcall_private__(self, t):
        r"""
        This ensures that a :class:`SemistandardSetValuedTableau` is only
        constructed as an ``element_class`` call of an appropriate parent.

        EXAMPLES::
        
            sage: t = SemistandardSetValuedTableau([[[1],[1],[8]],[[2],[6,2]],[[3,7,4]]])
            sage: t.shape()
            [3, 2, 1]

        TESTS::

            sage: SemistandardSetValuedTableau([[[1,2],[2,3]],[[4,6]]])
            [[[1, 2], [2, 3]], [[4, 6]]]

        """

        if isinstance(t, SemistandardSetValuedTableau):
            return t

        t_list = list(t)
        for i in range(len(t)):
            for j in range(len(t[i])):
                #t_list[i][j] =  tuple(sorted(t[i][j]))
                t_list[i][j] =  sorted(t[i][j])

        SSVT = SemistandardSetValuedTableaux(Tableau(t_list).shape())
        return SSVT.element_class(SSVT, t_list)

    def check(self):
        """
        Check that ``self`` is a valid semistandard set-valued tableau.
        """
        super(SemistandardSetValuedTableau, self).check()

        # Tableau() has checked that t is tableau, so it remains to check that
        # the entries of t are positive integers which are weakly increasing
        # along rows
        for row in self:
            try:
                row_t = [tuple(_) for _ in row]
            except TypeError:
                raise ValueError("the entries of a semistandard set-valued tableau must be iterables")
            if any(max(row[c]) > min(row[c+1]) for c in range(len(row)-1)):
                raise ValueError("the entries in each row of a semistandard set-valued tableau must be weakly increasing")

        # and strictly increasing down columns
        if self:
            for row, next in zip(self, self[1:]):
                if not all(max(row[c]) < min(next[c]) for c in range(len(next))):
                    raise ValueError("the entries of each column of a semistandard set-valued tableau must be strictly increasing")

    def excess(self):
        r"""
        Return the excess statistic for ``self``.

        The excess of a semistandard set-valued tableaux ``T`` is the total number 
        of integers in ``T`` minus the size of ``T``.

        EXAMPLES::

            sage: T = SemistandardSetValuedTableau([[[1,2],[2,3]],[[3,4,5]]])
            sage: T.excess()
            4
        """
        tot = sum([len(cells) for cells in self.entries()])
        return tot - self.size()

    def weight(self):
        r"""
        Return the weight of the set-valued tableau ``self``. Trailing zeroes are
        omitted when returning the weight.

        The weight of a tableau `T` is the sequence `(a_1, a_2, a_3, \ldots )`,
        where `a_k` is the number of entries of `T` equal to `k`. This
        sequence contains only finitely many nonzero entries.

        EXAMPLES::

            sage: SemistandardSetValuedTableau([[[1],[1],[8]],[[2],[6,2]],[[3,7,4]]]).weight()
            [2, 2, 1, 1, 0, 1, 1, 1]

            sage: SemistandardSetValuedTableau([]).weight()
            []
        """
        if len(self) == 0:
            return []
        vec = sum([sum(row,[]) for row in self],[])
        m = max(vec)
        wt = [0] * m
        for i in vec:
            if i > 0:
                wt[i - 1] += 1
        return wt

    @combinatorial_map(order=2, name='Bender-Knuth involution')
    def bender_knuth_involution(self, k, rows=None):
        r"""
        Return the image of ``self`` under the `k`-th Bender--Knuth
        involution, assuming ``self`` is a semistandard set-valued tableau.
        This function was introduced by Ikeda and Shimazaki [IS2014]_ in this context.

        Let `T` be a tableau and fix `k`. Then a free `k` in `T` means a cell of
        `T` that contains `k` and whose direct lower
        neighbor does not contain `k + 1` (in particular,
        this lower neighbor might not exist). A free `k + 1`
        in `T` is a cell of `T` that contains `k + 1`
        and whose direct upper neighbor does not contain `k`
        (in particular, this neighbor might not exist). Note that a cell
        that contains both `k` and `k + 1` is both a free `k` and a free `k + 1`. 
        It is clear that for any row `r` of `T`, the free `k`'s and
        free `k + 1`'s in `r` together form a contiguous interval of `r`.

        The *`k`-th Bender--Knuth involution at row `i`* changes the entries of
        the cells in this interval in such a way that if it used to have
        `a` entries of `k` and `b` entries of `k + 1`, it will now
        have `b` entries of `k` and `a` entries of `k + 1`. For fixed `k`, the
        `k`-th Bender--Knuth switches for different `i` commute. The
        composition of the `k`-th Bender--Knuth switches for all rows is
        called the *`k`-th Bender--Knuth involution*. This is used to show that
        the symmetric Grothendieck polynomials defined as generating functions for
        semistandard set-valued tableaux are in fact symmetric polynomials.

        INPUT:

        - ``k`` -- an integer

        - ``rows`` -- (Default ``None``) When set to ``None``, the method
          computes the `k`-th Bender--Knuth involution as defined above.
          When an iterable, this computes the composition of the `k`-th
          Bender--Knuth switches at row `i` over all `i` in ``rows``. When set
          to an integer `i`, the method computes the `k`-th Bender--Knuth
          switch at row `i`. Note the indexing of the rows starts with `1`.

        OUTPUT:

        The image of ``self`` under either the `k`-th Bender--Knuth
        involution, the `k`-th Bender--Knuth switch at a certain row, or
        the composition of such switches, as detailed in the INPUT section.

        EXAMPLES::

            sage: t = SemistandardSetValuedTableau([[[1],[1,2,3]],[[4,6]]])
            sage: t
            [[[1], [1, 2, 3], [[4, 6]]
            sage: t.bender_knuth_involution(1)
            [[[1, 2], [2, 3], [[4, 6]]
            sage: t.bender_knuth_involution(1).bender_knuth_involution(2)
            [[[1, 2, 3], [3], [[4, 6]]

        The Bender--Knuth involution is an involution::

            sage: t = SemistandardSetValuedTableau([[[1],[1],[8]],[[2],[6,2]],[[3,7,4]]])
            sage: all(t.bender_knuth_involution(k).bender_knuth_involution(k) == t for k in range(1,8))
            True
        """
        ell = len(self)    # ell is the number of rows of self.
        # Sanitizing the rows input so that it always becomes a list of
        # nonnegative integers. We also subtract 1 from these integers
        # because the i-th row of a tableau T is T[i - 1].
        def rem(tup,num):
            mylist = list(tup)
            mylist.remove(num)
            return tuple(mylist)

        def app(tup,num):
            mylist = list(tup)
            mylist.append(num)
            return tuple(mylist)

        if rows is None:
            rows = list(range(ell))
        elif rows in ZZ:
            rows = [rows - 1]
        else:
            rows = [i - 1 for i in rows]
        # Now, rows should be iterable.

        # result_tab is going to be the result tableau (as a list of lists);
        # we will build it up step by step, starting with a deep copy of self.
        result_tab = self.to_list()
        for i in rows:
            if i >= ell:
                continue
            # Setup the previous and next rows
            if i == 0:
                prev_row = [[None]] * len(result_tab[i])
            else:
                prev_row = result_tab[i-1]
            if i == ell - 1:
                next_row = [[None]] * len(result_tab[i])
            else:
                next_row = result_tab[i+1] + [[None]] * (len(result_tab[i]) - len(result_tab[i+1]))
            a = 0 #counter for free k
            b = 0 #counter for free k+1
            c=0 #counter for free box with both k and k+1
            sk = None # The column number of the first free k
            sk1 = None # The column number of the first free k+1
            skboth = None #The column number of the (unique) free box with both k and k+1
            for j, val in enumerate(result_tab[i]):
                if k in val and k+1 not in next_row[j]:
                    if k+1 in val:
                        c += 1
                        skboth = j
                    else:
                        if sk is None:
                            sk = j
                        a += 1
                elif k+1 in val and k not in prev_row[j]:
                    if sk1 is None:
                        sk1 = j
                    b += 1
            if skboth is None:
                if sk1 is not None:
                    if a > b:
                        for j in range(sk1-(a-b), sk1):
                            result_tab[i][j] = rem(result_tab[i][j],k)
                            result_tab[i][j] = app(result_tab[i][j],k+1)
                    elif a < b:
                        for j in range(sk1, sk1+b-a):
                            result_tab[i][j] = rem(result_tab[i][j],k+1)
                            result_tab[i][j] = app(result_tab[i][j],k)
                elif sk is not None:
                    for j in range(sk, sk+a):
                        result_tab[i][j] = rem(result_tab[i][j],k)
                        result_tab[i][j] = app(result_tab[i][j],k+1)
            else:
                if sk1 is not None:
                    if a > b:
                        result_tab[i][sk1-(a-b)-1] = app(result_tab[i][sk1-(a-b)-1],k+1)
                        for j in range(sk1-(a-b), sk1):
                            result_tab[i][j] = rem(result_tab[i][j],k)
                            result_tab[i][j] = app(result_tab[i][j],k+1)
                        result_tab[i][sk1] = rem(result_tab[i][sk1],k)
                    elif a < b:
                        result_tab[i][sk1-1] = rem(result_tab[i][sk1-1],k+1)
                        for j in range(sk1, sk1+b-a-1):
                            result_tab[i][j] = rem(result_tab[i][j],k+1)
                            result_tab[i][j] = app(result_tab[i][j],k)
                        result_tab[i][sk1+b-a-1] = app(result_tab[i][sk1+b-a-1],k)
                elif sk is not None:
                    result_tab[i][sk] = app(result_tab[i][sk],k+1)
                    for j in range(sk+1, sk+a):
                        result_tab[i][j] = rem(result_tab[i][j],k)
                        result_tab[i][j] = app(result_tab[i][j],k+1)
                    result_tab[i][sk+a] = rem(result_tab[i][sk+a],k)
        return SemistandardSetValuedTableau(result_tab)

    def pp(self):
        r"""
        Returns the pretty print of self.

        EXAMPLES::

            sage: T = [[[1,2,3,4,6],[6],[6,7],[8,9],[9,11,12],[12]], [[7],[7],[8,9,10],[10,11,13,14],[14]], [[8,9],[9,10],[11,13],[16,17,18]]]
            sage: SemistandardSetValuedTableau(T).pp()
            [ 1,2,3,4,6 ][  6   ][  6,7   ][     8,9     ][ 9,11,12 ][ 12 ]
            [     7     ][  7   ][ 8,9,10 ][ 10,11,13,14 ][   14    ]
            [    8,9    ][ 9,10 ][ 11,13  ][  16,17,18   ]

            sage: U = SemistandardSetValuedTableau([[[1],[1],[8]],[[2],[6,2]],[[3,7,4]]])
            sage: U.pp()
            [   1   ][  1  ][ 8 ]
            [   2   ][ 2,6 ]
            [ 3,4,7 ]

        """
        max_len = max(len(row) for row in self)
        str_len = [[sum([len(str(elt)) for elt in cell])+len(cell)-1 for cell in row] for row in self]
        col_max = [max(row[j] for row in str_len if j<len(row)) for j in range(max_len)]
        S = ""
        for row in self:
            for j in range(len(row)):
                s = ""
                for k in range(len(row[j])):
                    s += str(row[j][k])
                    if k < len(row[j])-1:
                        s+=","
                S+="[ "+'{st:{c}^{n}}'.format(st=s,c=" ",n=col_max[j])+" ]"
            S += "\n"
        print(S)

class CrystalElementSemistandardSetValuedTableau(SemistandardSetValuedTableau):
    r"""
    Class for elements of ``crystals.SemistandardSetValuedTableaux``

    INPUT:

    - ``shape`` -- the shape of set-valued tableaux

    - ``n`` -- maximum entry in cells of tableaux
    
    If no maximum entry is specified, then 

    EXAMPLES::
        
        sage: B = crystals.SemistandardSetValuedTableaux([2,1],3); B
        Crystal of set-valued tableaux of type A_2 and shape [2,1]
        sage: C = crystals.SemistandardSetValuedTableaux([2,1,1,1]); C
        Crystal of set-valued tableaux of type A_3 and shape [2,1,1,1]    
    """

    # @staticmethod
    # def __classcall_private__(cls, shape, n=None):
    #     r"""
    #     Classcall to mend the input.

    #     TESTS::

    #         sage: B = crystals.SemistandardSetValuedTableaux([2,1],3); B
    #         Crystal of set-valued tableaux of type A_2 of shape [2,1]
    #         sage: C = crystals.SemistandardSetValuedTableaux([2,1,1,1]); C
    #         Crystal of set-valued tableaux of type A_3 of shape [2,1,1,1]
    #     """ 
    #     if not isinstance(shape,Partition) and not isinstance(shape,list):
    #         raise ValueError("shape should be a partition")
    #     else:
    #         shape = Partition(shape)
    #     if n is None:
    #         n = shape.size()
    #     elif n<=0:
    #         raise ValueError("n should be a positive integer")
    #     return super(SetValuedTableaux, cls).__classcall__(cls, shape, n)

    # def __init__(self, shape, n):
    #     r"""
    #     Initialize crystal of semistandard set-valued tableaux of a fixed shape and given maximum entry. 
        
    #     EXAMPLES::
            
    #         sage: B = crystals.SemistandardSetValuedTableaux([2,1],3)
    #         sage: B.n
    #         3
    #         sage: B.shape
    #         [2,1]

    #     TESTS::

    #         sage: B = crystals.SemistandardSetValuedTableaux([2,1],3)
    #         sage: TestSuite(B).run()
    #     """
    #     Parent.__init__(self, category = ClassicalCrystals())
    #     self.n = n
    #     self.shape = shape
    #     cartan_type = CartanType(['A',n-1])
    #     self._cartan_type = cartan_type
    #     # (this enumerates all highest weight vectors!)
    #     self.module_generators = [] 

    # def _repr_(self):
    #     r"""
    #     EXAMPLES::

    #         sage: B = crystals.SemistandardSetValuedTableaux([2,1],3); B
    #         Crystal of set-valued tableaux of type A_2 of shape [2,1]
    #     """
    #     return "Crystal of set-valued tableaux of type A_{} of shape {}".format(self.n-1, self.shape)

    # # temporary workaround while an_element is overriden by Parent
    # _an_element_ = EnumeratedSets.ParentMethods._an_element_

    # class Element(ElementWrapper):
    
    #     def __init__(self,parent,tab):
    #         r"""
    #         Initialize self as a crystal element.

    #         EXAMPLES::

    #         sage: SVT = crystals.SemistandardSetValuedTableau([2,1],3)
    #         sage: T = SVT([ [[1,2],[2,3]],[[3]] ]); T
    #         [[[1, 2], [2,3]], [[3]] ]
    #         """
    #         self.n = parent.n
    #         self.value = SemistandardSetValuedTableau(tab)
    #         ElementWrapper.__init__(self, parent, self.value)

    def _get_signs(self, i):
        """
        Auxiliary function for `e_i` and `f_i` methods.

        Assign each column of self a +1, -1 or 0 according to
        +1 if there is an unmatched `i+1` aka left paren '('
        -1 if there is an unmatched `i` aka right paren ')'
        0 if all the `i`s and `i+1`s are matched
        
        Return list of +1, -1, 0 with length equal to number of columns of self.
        """
        #st = SemistandardSetValuedTableau(self)
        st = self.value
        signs = []
        for col in st.conjugate():
            word = sum(col, ())
            if i in word and i+1 in word:
                signs += [0]
            elif i in word: #i in word, i+1 not in word
                signs += [-1]
            elif i+1 in word: # i not in word, i+1 is
                signs += [+1]
            else: # neither i nor i+1 in word
                signs += [0]
        return signs

    def _bracket(self, i, right=True):
        """
        auxiliary function for `e_i` and `f_i` methods.

        Return index of column in self with rightmost `i` to be changed to `i+1` 
        or leftmost `i+1` to be changed to `i`.

        If right is True (default), then return index of rightmost `i`.
        If right is False, then return index of leftmost `i+1`.

        If no `i` can be changed to `i+1` or vice versa, return -1.
        """
        x = self._get_signs(i) # x is a list of +1, -1, 0
        if not right:
            x = [-j for j in x][::-1]
        count = 0
        index = -1
        for j in range(len(x)):
            if x[j] == -1:
                if count == 0:
                    index = j
                else:
                    count -= 1
            if x[j] == +1:
                count += 1
        if right:
            return index
        else:
            return -1 if index < 0 else len(x)-1-index

    def e(self, i):
        r"""
        Returns the action of `e_i` on ``self``.

        EXAMPLES::

            sage: 
        """
        if i not in self.index_set():
            raise ValueError("i must be in the index set")
        col = self._bracket(i,right=False)
        if col == -1:
            return None
        tab = [[list(entry) for entry in r] for r in self.value]
        t = self.value.conjugate()
        column = t[col]
        row = min([ j for j in range(len(column)) if i+1 in column[j] ])
        # checks that there is a cell to the left and that the cell contains i and i+1
        if col>0 and all(x in tab[row][col-1] for x in [i,i+1]):
            tab[row][col-1].remove(i+1)
        else:
            tab[row][col].remove(i+1)
        tab[row][col] = sorted(tab[row][col]+[i])
        return self.parent()(tab)
        
    def f(self, i):
        r"""
        Returns the action of `f_i` on ``self``.

        EXAMPLES::

            sage: 
        """
        if i not in self.index_set():
            raise ValueError("i must be in the index set")
        col = self._bracket(i,right=True)
        if col == -1:
            return None
        tab = [[list(entry) for entry in r] for r in self.value]
        t = self.value.conjugate()
        column = t[col]
        row = min([ j for j in range(len(column)) if i in column[j] ])
        # checks that there is a cell to the right and that the cell contains i and i+1
        if col<len(self.value[row])-1 and all(x in tab[row][col+1] for x in [i,i+1]):
            tab[row][col+1].remove(i)
        else:
            tab[row][col].remove(i)
        tab[row][col] = sorted(tab[row][col]+[i+1])
        return self.parent()(tab)

    def reading_word(self):
        pass

class SemistandardSetValuedTableaux(Tableaux):
    r"""
    Return the class of semistandard set-valued tableaux.

    .. WARNING::

        Giving no shape or size is currently not implemented. Will need to add a
        class ``SemistandardSetValuedTableaux_all``

    INPUT:
    
    Positional arguments:

    - ``p`` -- The first argument is either a non-negative integer or a partition
    
    Keyword arguments:

    - ``max_entry`` -- a positive integer that is the maximum allowed entry in the tableau

    OUTPUT:

    - With no shape or size, the class of all semistandard set-valued tableaux.

    - With no ``max_entry``, the class of all semistandard set-valued tableaux of 
      specified size or shape.

    - With a ``max_entry``, one of the following:
        - With a non-negative integer amount, ``p``, the class of all semistandard 
          set-valued tableaux of size ``p`` and maximum integer ``max_entry``.

        - With a partition argument, ``p``, the class of all semistandard 
          set-valued tableaux of shape ``p`` and maximum integer ``max_entry``.

    A semistandard set-valued tableau is a tableau whose cells are filled with 
    nonempty sets of integers such that the tableau is semistandard with respect 
    to the ordering on sets A, B given by A <= B if max(A) < min(B).

    EXAMPLES::

        sage: SSVT = SemistandardSetValuedTableaux(3); SSVT
        Semistandard set-valued tableaux of size 3
        sage: SSVT = SemistandardSetValuedTableaux(3, max_entry=2); SSVT
        Semistandard set-valued tableaux of size 3 and max entry 2
        sage: list(SSVT)
        [[[[1], [1], [1]]],
         [[[1], [1], [1, 2]]],
         [[[1], [1], [2]]],
         [[[1], [1, 2], [2]]],
         [[[1], [2], [2]]],
         [[[1, 2], [2], [2]]],
         [[[2], [2], [2]]],
         [[[1], [1]], [[2]]],
         [[[1], [1, 2]], [[2]]],
         [[[1], [2]], [[2]]]]
        sage: SSVT = SemistandardSetValuedTableaux([2,2], max_entry=3); SSVT
        Semistandard set-valued tableaux of shape [2, 2] and max entry 3
    """
    @staticmethod
    def __classcall_private__(cls, *args, **kwargs):
        r"""
        TESTS::

            sage: SemistandardSetValuedTableaux()
            Semistandard set-valued tableaux
            sage: SemistandardSetValuedTableaux(max_entry=4)
            Semistandard set-valued tableaux of max entry 4
            sage: SemistandardSetValuedTableaux(2)
            Semistandard set-valued tableaux of size 2
            sage: SemistandardSetValuedTableaux(2, max_entry=4)
            Semistandard set-valued tableaux of size 2 and max entry 4
            sage: SemistandardSetValuedTableaux([4,2,1])
            Semistandard set-valued tableaux of shape [4, 2, 1]
            sage: SemistandardSetValuedTableaux([4,2,1], max_entry=3)
            Semistandard set-valued tableaux of shape [4, 2, 1] and max entry 3
            sage: SemistandardSetValuedTableaux([1], max_entry=None)
            Semistandard set-valued tableaux of shape [1] and max entry 1
            sage: from sage.rings.infinity import PlusInfinity
            sage: SemistandardSetValuedTableaux([1], max_entry=PlusInfinity())
            Semistandard set-valued tableaux of shape [1]
            sage: SemistandardSetValuedTableaux([])
            Semistandard set-valued tableaux of shape []

            sage: SemistandardSetValuedTableaux(7, max_entry=0)
            Traceback (most recent call last):
            ...
            ValueError: the maximum entry must be a positive integer, None, or PlusInfinity
            sage: SemistandardSetValuedTableaux([2], max_entry=0)
            Traceback (most recent call last):
            ...
            ValueError: the maximum entry must be a positive integer, None, or PlusInfinity

            sage: SemistandardSetValuedTableaux({1,2,3})
            Traceback (most recent call last):
            ...
            ValueError: the argument must be a non-negative integer or a partition
        """
        max_entry = kwargs.get('max_entry')
        if args:
            p = args[0]

            # if p is size
            if isinstance(p,(int,Integer)) and p>=0:
                if max_entry is None or max_entry==PlusInfinity():
                    return SemistandardSetValuedTableaux_size_inf(Integer(p))
                elif isinstance(max_entry, (int,Integer))==False:
                    raise ValueError("the maximum entry must be a positive integer, None, or PlusInfinity")
                elif isinstance(max_entry, (int,Integer)) and max_entry <= 0:
                    raise ValueError("the maximum entry must be a positive integer, None, or PlusInfinity")
                return SemistandardSetValuedTableaux_size(Integer(p), max_entry)

            # if p is shape
            if p in _Partitions:
                if max_entry is None or max_entry==PlusInfinity():
                    return SemistandardSetValuedTableaux_shape_inf(_Partitions(p))
                elif isinstance(max_entry, (int,Integer)) and max_entry <= 0:
                    raise ValueError("the maximum entry must be a positive integer, None, or PlusInfinity")
                return SemistandardSetValuedTableaux_shape(_Partitions(p), max_entry)
            else:
                raise ValueError("the argument must be a non-negative integer or a partition")            
        else:
            if max_entry is not None and max_entry!=PlusInfinity() and isinstance(max_entry, (int,Integer))==False:
                raise ValueError("the maximum entry must be a positive integer, None, or PlusInfinity")
            if isinstance(max_entry, (int,Integer)) and max_entry <= 0:
                raise ValueError("the maximum entry must be a positive integer, None, or PlusInfinity")
            return SemistandardSetValuedTableaux_all(max_entry=max_entry)

    Element = SemistandardSetValuedTableau

    def __init__(self, *args, **kwargs):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: S = SemistandardSetValuedTableaux(2, max_entry=3)
            sage: TestSuite(S).run()
        """
        if 'max_entry' in kwargs:
            if kwargs['max_entry']==PlusInfinity():
                self.max_entry = None
            else:
                self.max_entry = kwargs['max_entry']
            kwargs.pop('max_entry')
        else:
            self.max_entry = None

        Tableaux.__init__(self, *args, **kwargs)

    def __contains__(self,t):
        r"""
        Determine if t is an element of self.

        TESTS::

            sage: T1 = [ [ [1],[1,2] ],[ [2] ] ]
            sage: T2 = Tableau([ [ [1],[1,3] ],[ [2] ] ])
            sage: T3 = Tableau([ [ [1,2],[2] ],[ [2] ] ])
            sage: T1 in SemistandardSetValuedTableaux(3)
            True
            sage: T2 in SemistandardSetValuedTableaux(3)
            True
            sage: T3 in SemistandardSetValuedTableaux(3)
            True
            sage: T2 in SemistandardSetValuedTableaux([2,1])
            True
            sage: T1 in SemistandardSetValuedTableaux([2,1,1])
            False
            sage: T1 in SemistandardSetValuedTableaux(3, max_entry=2)
            True
            sage: T2 in SemistandardSetValuedTableaux(3, max_entry=2)
            False

            sage: ssvts = []
            sage: for st in StandardTableaux([3,1,1]):
            ....:     ssvts.append([[[_] for _ in row] for row in st])
            ....:        
            sage: all(i in SemistandardSetValuedTableaux([3,1,1]) for i in ssvts)
            True
        """
        if isinstance(t,SemistandardSetValuedTableau):
            return True
        # t is assumed to be at least a list of lists with shape given by a partition
        elif Tableaux.__contains__(self,t):         
            for row in t:
                for cell in row:
                    # checks that cell contains a set
                    if not isinstance(cell,(list,tuple,set)) or len(cell)>len(set(cell)):
                        return False
                    # checks that the set consists of nonempty set of integers
                    if len(cell)==0 or not all(isinstance(elt,(int,Integer)) and elt>0 for elt in cell):
                        return False
                # checks cells are weakly increasing along rows
                for i in range(len(row)-1):
                    left,right = row[i],row[i+1]
                    if max(left) > min(right):
                        return False
            # checks cells are strictly increasing along columns
            for up, down in zip(t, t[1:]):
                if any(max(up[i])>=min(down[i]) for i in range(len(down))):
                    return False 
            # if self has max_entry, checks if that all elements in each cell are within max_entry
            return self.max_entry is None or max(max(max(cell) for cell in row) for row in t)<=self.max_entry

class SemistandardSetValuedTableaux_all(SemistandardSetValuedTableaux):
    """
    Class of all semistandard set-valued tableaux.
    """
    def __init__(self, max_entry=None):
        r"""
        Initializes the class of all semistandard set-valued tableaux.

        TESTS::

            sage: SSVT = SemistandardSetValuedTableaux()
            sage: [[[1,2],[2]],[[3]]] in SSVT
            True
            sage: [[[1,2],[2]],[[2]]] in SSVT
            False
            sage: [[[1,2],[1,3]],[[3]]] in SSVT
            False
            sage: [[[1,2],[]],[[3]]] in SSVT
            False
            sage: TestSuite(SSVT).run()  # long time
        """
        if max_entry==PlusInfinity() or max_entry is None:
            SemistandardSetValuedTableaux.__init__(self,category=InfiniteEnumeratedSets())
        else:
            SemistandardSetValuedTableaux.__init__(self,max_entry=max_entry,category=InfiniteEnumeratedSets())

    def _repr_(self):
        """
        TESTS::

            sage: SemistandardSetValuedTableaux()
            Semistandard set-valued tableaux
            sage: SemistandardSetValuedTableaux(max_entry=4)
            Semistandard set-valued tableaux of max entry 4
        """
        if self.max_entry is None:
            return "Semistandard set-valued tableaux"
        return "Semistandard set-valued tableaux of max entry {}".format(self.max_entry)

    def __contains__(self, x):
        """
        TESTS::     
            sage: SSVT1 = SemistandardSetValuedTableaux()
            sage: [[[1,2]],[[2]],[[4]]] in SSVT1
            False
            sage: [[[1,2],[1,3]],[[2]],[[4]]] in SSVT1
            False
            sage: [[[1,2]],[[3]],[[4]]] in SSVT1
            True
            sage: [[[1,2],[2]],[[3]],[[4]]] in SSVT1
            True
            sage: [[[1,2],[3],[40]]] in SSVT1
            True
            sage: [[[1]],[[3,40]],[[50]]] in SSVT1
            True

            sage: SSVT2 = SemistandardSetValuedTableaux(max_entry=4)
            sage: [[[1,2]],[[2]],[[4]]] in SSVT2
            False
            sage: [[[1,2],[1,3]],[[2]],[[4]]] in SSVT2
            False
            sage: [[[1,2]],[[3]],[[4]]] in SSVT2
            True
            sage: [[[1,2],[2]],[[3]],[[4]]] in SSVT2
            True
            sage: [[[1,2],[3],[4]]] in SSVT2
            True
            sage: [[[1]],[[3,4]],[[5]]] in SSVT2
            False
        """

        return SemistandardSetValuedTableaux.__contains__(self, x)

class SemistandardSetValuedTableaux_size(SemistandardSetValuedTableaux):
    """
    Class of all semistandard set-valued tableaux of fixed size ``n``.
    """
    def __init__(self, n, max_entry=None):
        r"""
        Initializes the class of all semistandard set-valued tableaux of size ``n``.

        .. WARNING::

            Input is not checked; please use :class:`SemistandardSetValuedTableaux` to
            ensure the options are properly parsed.

        TESTS::

            sage: TestSuite( SemistandardSetValuedTableaux(0,max_entry=4) ).run()
            sage: TestSuite( SemistandardSetValuedTableaux(3,max_entry=4) ).run()
        """
        if max_entry is None:
            self.max_entry = n
        SemistandardSetValuedTableaux.__init__(self,max_entry=max_entry,category=FiniteEnumeratedSets())
        self._size = Integer(n)

    def _repr_(self):
        """
        TESTS::

            sage: SemistandardSetValuedTableaux(3,max_entry=4)
            Semistandard set-valued tableaux of size 3 and max entry 4
        """
        return "Semistandard set-valued tableaux of size {} and max entry {}".format(self._size, self.max_entry)

    def __contains__(self, x):
        """
        TESTS::

            sage: SSVT = SemistandardSetValuedTableaux(3,max_entry=4)
            sage: [[[1,2]],[[2]],[[4]]] in SSVT
            False
            sage: [[[1,2],[1,3]],[[2]],[[4]]] in SSVT
            False
            sage: [[[1,2]],[[3]],[[4]]] in SSVT
            True
            sage: [[[1,2],[2]],[[3]],[[4]]] in SSVT
            False
            sage: [[[1,2],[3],[4]]] in SSVT
            True
            sage: [[[1]],[[3,4]],[[5]]] in SSVT
            False
        """
        return SemistandardSetValuedTableaux.__contains__(self, x) and \
        sum(map(len, x)) == self._size and \
        max([max([max(cell) for cell in row]) for row in x]) <= self.max_entry

    def __iter__(self):
        """
        EXAMPLES::

            sage: SSVT = SemistandardSetValuedTableaux(3,max_entry=2); SSVT
            Semistandard set-valued tableaux of size 3 and max entry 2
            sage: SSVT.list()
            [[[[1], [1], [1]]],
            [[[1], [1], [1, 2]]],
            [[[1], [1], [2]]],
            [[[1], [1, 2], [2]]],
            [[[1], [2], [2]]],
            [[[1, 2], [2], [2]]],
            [[[2], [2], [2]]],
            [[[1], [1]], [[2]]],
            [[[1], [1, 2]], [[2]]],
            [[[1], [2]], [[2]]]]
            sage: len(SSVT)
            10
        """
        if self._size == 0:
            yield self.element_class(self, [])

        for part in Partitions(self._size):
            for ssvt in SemistandardSetValuedTableaux_shape(part, self.max_entry):
                yield self.element_class(self, ssvt)

class SemistandardSetValuedTableaux_size_inf(SemistandardSetValuedTableaux, DisjointUnionEnumeratedSets):
    """
    Class of all semistandard set-valued tableaux of fixed size `n` and no maximum entry.
    """
    def __init__(self, n):
        r"""
        Initializes the class of all semistandard set-valued tableaux of size ``n`` and no maximum entry.

        .. WARNING::

            Input is not checked; please use :class:`SemistandardSetValuedTableaux` to
            ensure the options are properly parsed.

        TESTS::

            sage: TestSuite( SemistandardSetValuedTableaux(0) ).run()
            sage: TestSuite( SemistandardSetValuedTableaux(3) ).run()
        """
        SemistandardSetValuedTableaux.__init__(self,category=InfiniteEnumeratedSets())
        DisjointUnionEnumeratedSets.__init__(self,
                                             Family(Partitions_n(n),SemistandardSetValuedTableaux_shape_inf),
                                             facade=True, keepkey=False)
        self._size = Integer(n)

    def _repr_(self):
        """
        TESTS::

            sage: SemistandardSetValuedTableaux(3)
            Semistandard set-valued tableaux of size 3
        """
        return "Semistandard set-valued tableaux of size {}".format(self._size)

    def __contains__(self, x):
        """
        TESTS::

            sage: SSVT = SemistandardSetValuedTableaux(3)
            sage: [[[1,2]],[[2]],[[3]]] in SSVT 
            False
            sage: [[[1,2]],[[3]],[[4]]] in SSVT 
            True
            sage: [[[1,2],[2]],[[3]],[[4]]] in SSVT 
            False
            sage: [[[1,2],[2],[3]]] in SSVT 
            True
        """
        return SemistandardSetValuedTableaux.__contains__(self, x) and sum(map(len, x)) == self._size

class SemistandardSetValuedTableaux_shape(SemistandardSetValuedTableaux):
    """
    Class of all semistandard set-valued tableaux of a fixed shape `p`.
    """
    def __init__(self, p, max_entry):
        r"""
        Initialize the class of all semistandard set-valued tableaux of given shape ``p``.

        .. WARNING::

            Input is not checked; please use :class:`SemistandardSetValuedTableaux` to
            ensure the options are properly parsed.

        TESTS::

            sage: TestSuite( SemistandardSetValuedTableaux([2,1,1],max_entry=3) ).run()
        """
        if max_entry is None:
            self.max_entry = len(p)
        SemistandardSetValuedTableaux.__init__(self,max_entry=max_entry,category=ClassicalCrystals())
        self._shape = p
        self._cartan_type = CartanType(['A', max_entry-1])
        self.Element = CrystalElementSemistandardSetValuedTableau

    def __contains__(self, x):
        """
        EXAMPLES::
        
            sage: SSVT = SemistandardSetValuedTableaux([2,1,1],max_entry=6)
            sage: [[[1,2],[2,3]],[[2]],[[3]]] in SSVT
            False
            sage: [[[1,2],[2,3]],[[3,4]],[[6]]] in SSVT
            True
            sage: [[[1,2],[2,3]],[[3,6]],[[6,7]]] in SSVT
            False
            sage: [[[1,2],[1,3]],[[3,4]],[[6]]] in SSVT
            False
            sage: [[[1,2],[2,3]],[[3,4],[4,5,6]],[[6]]] in SSVT
            False
        """
        return SemistandardSetValuedTableaux.__contains__(self, x) and \
               [len(row) for row in x] == self._shape and \
               max([max([max(cell) for cell in row]) for row in x]) <= self.max_entry

    def _repr_(self):
        """
        TESTS::

            sage: SemistandardSetValuedTableaux([3,1],max_entry=2)
            Semistandard set-valued tableaux of shape [3, 1] and max entry 2 
        """
        return "Semistandard set-valued tableaux of shape {} and max entry {}".format(self._shape, self.max_entry)

    @lazy_attribute
    def module_generators(self):
        """
        Return the generators of ``self`` as a crystal.

        TESTS::

            sage: 
        """
        # if self._skew is not None:
        #     raise NotImplementedError("only for non-skew shapes")
        # list_dw = []
        # if self._max_entry is None:
        #     max_entry = sum(self._shape)
        # else:
        #     max_entry = self._max_entry
        # for weight in (Partition(self._shape).dominated_partitions(rows=max_entry)):
        #     list_dw.extend([self.element_class(self, T, check=False,
        #                                        preprocessed=True)
        #                     for T in ShiftedPrimedTableaux(weight=tuple(weight),
        #                                                    shape=self._shape)])
        # return tuple(list_dw)

    def shape(self):
        """
        Return the shape of the semistandard set valued tableau ``self``.

        TESTS::

            sage: SemistandardSetValuedTableaux([3,1],max_entry=2).shape()
            [3, 1]
        """
        return self._shape

    def __iter__(self):
        r"""
        An iterator for the semistandard set-valued tableaux associated to the
        shape ``p`` of ``self``.

        .. WARNING::

            if ``max_entry`` is None, will raise a ValueError.

        EXAMPLES::

            sage: SSVT = SemistandardSetValuedTableaux([2,2],max_entry=3); SSVT
            Semistandard set-valued tableaux of shape [2, 2] and max entry 3
            sage: SSVT.list()
            [[[[1], [1]], [[2], [2]]],
            [[[1], [1]], [[2], [2, 3]]],
            [[[1], [1]], [[2], [3]]],
            [[[1], [1]], [[2, 3], [3]]],
            [[[1], [1, 2]], [[2], [3]]],
            [[[1], [1, 2]], [[2, 3], [3]]],
            [[[1], [1]], [[3], [3]]],
            [[[1], [1, 2]], [[3], [3]]],
            [[[1], [2]], [[2], [3]]],
            [[[1], [2]], [[2, 3], [3]]],
            [[[1], [2]], [[3], [3]]],
            [[[1, 2], [2]], [[3], [3]]],
            [[[2], [2]], [[3], [3]]]]
            sage: len(SSVT)
            13

        """
        if self.max_entry is None:
            raise ValueError("max_entry must be specified to iterate")

        def jumps(row):
            """
            Return list of indices i where row[i] < row[i+1] and also
            the last index of row.

            INPUT:
            
            - ``row`` -- a list of list of integers
            """
            return [i for i in range(len(row)-1) if max(row[i])<min(row[i+1])] + [len(row)-1]
        
        def addable(cell,right,below,max_entry):
            """
            Return a list of numbers that can added to the cell ``cell`` of a
            semistandard set-valued tableau.

            INPUT:

            - ``cell`` -- a nonempty list of integers
            - ``right`` -- a list of integers or None
            - ``below`` -- a list of integers or None
            - ``max_entry`` -- a nonnegative integer

            """
            cell_val = max(cell)+1
            right_val = min(right)+1 if right is not None else max_entry+1
            below_val = min(below) if below is not None else max_entry+1

            if right_val < below_val:
                return list(range(cell_val, right_val))
            else:
                return list(range(cell_val, below_val))

        from sage.misc.mrange import cartesian_product_iterator as CPI
        import copy
        for t in SemistandardTableaux(shape=self._shape, max_entry=self.max_entry):
            tab = [[[entry] for entry in row] for row in t]
            addable_indices = []
            addable_sets = []
            for i in range(len(tab)):
                J = jumps(tab[i])
                for j in J:
                    cell = tab[i][j]
                    right, below = None, None
                    if i+1 < len(tab):
                        if j < len(tab[i+1]):
                            below = tab[i+1][j]
                    if j+1 < len(tab[i]):
                        right = tab[i][j+1]

                    # if a cell can be added with extra letters
                    addable_indices += [(i,j)]
                    elts_to_add = [tuple(_) for _ in powerset(addable(cell,right,below,self.max_entry))]
                    addable_sets.append(elts_to_add)
            
            for cp in CPI(addable_sets):
                tab_copy = copy.deepcopy(tab)
                for k in range(len(addable_sets)):
                    i,j = addable_indices[k]
                    tab_copy[i][j] += cp[k]
                yield self.element_class(self, tab_copy)

class SemistandardSetValuedTableaux_shape_inf(SemistandardSetValuedTableaux):
    """
    Class of all semistandard set-valued tableaux of a fixed shape `p` and no maximum entry.
    """
    def __init__(self, p):
        r"""
        Initializes the class of all semistandard set-valued tableaux of a given 
        shape and no maximum entry.

        .. WARNING::

            Input is not checked; please use :class:`SemistandardSetValuedTableaux` to
            ensure the options are properly parsed.

        TESTS::

            sage: TestSuite( SemistandardSetValuedTableaux([2,1,1]) ).run()
        """
        SemistandardSetValuedTableaux.__init__(self,category=InfiniteEnumeratedSets())
        self._shape = p

    def __contains__(self, x):
        """
        EXAMPLES::
        
            sage: SSVT = SemistandardSetValuedTableaux([2,1,1])
            sage: [[[1,2],[2,3]],[[2]],[[3]]] in SSVT
            False
            sage: [[[1,2],[2,3]],[[3,4]],[[6]]] in SSVT
            True
            sage: [[[1],[2,3,4]],[[3,10]],[[100]]] in SSVT
            True
            sage: [[[1,2],[1,3]],[[3,4]],[[6]]] in SSVT
            False
            sage: [[[1,2],[2,3]],[[3,4],[4,5,6]],[[6]]] in SSVT
            False
        """
        return SemistandardSetValuedTableaux.__contains__(self, x) and [len(row) for row in x] == self._shape

    def _repr_(self):
        """
        TESTS::

            sage: SemistandardSetValuedTableaux([3,1])
            Semistandard set-valued tableaux of shape [3, 1] 
        """
        return "Semistandard set-valued tableaux of shape {}".format(self._shape)

######################
#  Helper functions  #
######################

def reconstruct_tableau(seq):
    pass

def crowding_reverse_insertion(P,Q):
    pass

def crowding_map(P,Q):
    pass

def insertion_sequence(T):
    r"""
    Returns a sequence of words to insert in the uncrowding map.

    The algorithm assumes that entries in cells of T are sorted in increasing order.
    
        EXAMPLES:
            sage: T = SemistandardSetValuedTableau([ [ [2,3,1],[6,3,5],[7],[11,8,12,10] ], [ [4,5,6,7],[7],[9,12,8] ], [ [8],[9,8],[13] ], [ [9,10] ] ])
            sage: insertion_sequence(T)
            [[10, 9], [8, 9, 13, 8], [7, 7, 12, 9, 8, 6, 5, 4], [3, 6, 7, 12, 11, 10, 8, 5, 3, 2, 1]]
    """
    if T not in SemistandardSetValuedTableaux():
        raise ValueError("Semistandard set-valued tableau not given")

    if len(T)==0:
        return []
    
    seq = []
    for row in T[::-1]:
        S = []
        M = [cell[-1] for cell in row[:-1]]
        S += M
        H = row[-1][::-1]
        for cell in row[:-1][::-1]:
            H += cell[:-1][::-1]
        S += H
        seq += [S]
    return seq

def uncrowding_insertion(seq,P=None,Q=None):
    r"""
    Returns the pair (P',Q') under the uncrowding map for sequence seq and pair (P,Q).

        INPUT::
            seq - a sequence of integers to be inserted. This should insert to a hook shape.
            P - a semistandard Young tableau; empty if none initialized
            Q - a flagged increasing tableau with same shape as P; empty if none initialized

        OUTPUT::
            P' - a semistandard Young tableau
            Q' - a flagged increasing tableau with same shape as P'

        EXAMPLES::
            sage: P = SemistandardTableau([[4],[6],[7]])
            sage: Q = Tableau([['X'],['X'],[1]])
            sage: uncrowding_insertion([3,3,9,8],P,Q)
            ([[3, 3, 8], [4, 9], [6], [7]], [['X', 'X', 'X'], ['X', 1], ['X'], [1]])

            sage: uncrowding_insertion([1,1,4,5,4,3,2])
            ([[1, 1, 2, 4], [3], [4], [5]], [['X', 'X', 'X', 'X'], [1], [2], [3]])
    """
    Pp = P
    Qq = Q
    if P is None:
        Pp = SemistandardTableau([])
    if Q is None:
        Qq = Tableau([])
    if not isinstance(Pp,SemistandardTableau):
        raise ValueError("P should be instance of SemistandardTableau")
    if not isinstance(Qq,Tableau):
        raise ValueError("Q should be instance of Tableau")
    if Pp.shape()!=Qq.shape():
        raise ValueError("P and Q must be of same shape")        

    for x in seq:
        Pp = Pp.bump(x)
    shape = Pp.shape()
    temp = Tableau([ ['X']*shape[0] ]+Qq.to_list())
    cells = [pair for pair in Pp.cells() if pair not in temp.cells()]
    for cell in cells: 
        temp = temp.add_entry(cell,cell[0])
    Qq = temp
    return Pp, Qq

def uncrowding_map(T):
    r"""
    Returns the image of semistandard set-valued tableau T under the uncrowding map.

        EXAMPLE::
            sage: T = Tableau([ [ [1],[1,2,3] ],[ [2,3] ] ])
            sage: uncrowding_map(T)
            ([[1, 1], [2, 2], [3, 3]], [['X', 'X'], ['X', 1], [1, 2]])

            sage: T = Tableau([ [ [1],[1,2],[2] ],[ [2,3],[3,4,5] ],[ [4] ] ])
            sage: uncrowding_map(T)
            ([[1, 1, 2], [2, 2], [3, 3], [4, 4], [5]], [['X', 'X', 'X'], ['X', 'X'], ['X', 1], [2, 3], [3]])
    """
    P = SemistandardTableau([])
    Q = Tableau([])
    sequences=insertion_sequence(T.to_list(),sorted=True)
    for seq in sequences:
        P,Q = uncrowding_insertion(seq,P,Q)
    return P,Q