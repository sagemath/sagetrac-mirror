r"""
Multiset tableaux

AUTHORS:

- Wencin Poh, Jeremy Meza, Harrison Chapman (2019): initial version
"""

#*****************************************************************************
#       Copyright (C) 2019 Wencin Poh
#                          Jeremy Meza
#                          Harrison Chapman <hchaps at gmail dot com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#****************************************************************************

from collections import defaultdict

from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.combinat.tableau import (SemistandardTableau, SemistandardTableaux,
                                   Tableau, Tableaux)
from sage.combinat.partition import Partition, Partitions
from sage.combinat.permutation import Arrangements
from sage.combinat.vector_partition import VectorPartitions

class SemistandardMultisetTableau(Tableau):
    """
    A class to model a semistandard multiset tableau.
    """

    @staticmethod
    def __classcall_private__(self, t, order=None, key=None):
        r"""
        Ensures that a SMT is only ever created from a class of SMT
        """
        if isinstance(t, SemistandardMultisetTableau):
            # TODO do we want to try to change the order fn of t?
            return t
        # TODO implement this
        #elif t in SemistandardMultisetTableaux(key=key_fn):
        #    return SemistandardTableaux_all().element_class(SemistandardTableaux_all(), t)

        # t is not a semistandard tableau so we give an appropriate error message
        if t not in Tableaux():
            raise ValueError('%s is not a tableau' % t)

        if not all(isinstance(c, (int, Integer)) and c > 0
                   for row in t for sset in row for c in sset):
            raise ValueError("entries must be finite iterables of positive integers"%t)

        if any(key_fn(row[ci]) > key_fn(row[ci+1]) for row in t for ci in range(len(row)-1)):
            raise ValueError("The rows of %s are not weakly increasing"%t)

        # If we're still here ``t`` cannot be column strict
        # TODO are there other edge cases for SSMT?
        raise ValueError('%s is not a column strict tableau' % t)

    def check(self):
        """
        Check that ``self`` is a valid semistandard multiset tableau, for
        multiset order provided
        """
        super(SemistandardMultisetTableau, self).check()

        # Tableau() has checked that t is tableau, so it remains to check that
        # the entries of t are positive integers which are weakly increasing
        # along rows
        from sage.sets.positive_integers import PositiveIntegers
        PI = PositiveIntegers()

        for row in self:
            if any(self.key(row[c]) > self.key(row[c+1])
                   for c in range(len(row)-1)):
                raise ValueError("the entries in each row of a semistandard"
                                 "multiset tableau must be weakly increasing")

        # and strictly increasing down columns
        if self:
            for row, next in zip(self, self[1:]):
                if not all(self.key(row[c]) < self.key(next[c])
                           for c in range(len(next))):
                    raise ValueError(
                        "the entries of each column of a semistandard"
                        "multiset tableau must be strictly increasing")

class SemistandardMultisetTableaux(Tableaux):
    r"""
    Class of semistandard multiset tableaux.

    EXAMPLES:

        SemistandardMultisetTableaux(shape, content, order=my_compare)
    """
    @staticmethod
    def __classcall_private__(cls, *args, **kwargs):
        """
        """
        from sage.combinat.partition import _Partitions

        if args:
            n = args[0]
        elif 'n' in kwargs:
            n = kwargs['n']
        else:
            n = None

        weight = kwargs.get('weight', None)
        order = kwargs.get('order', 'grlex')

        # neither shape nor size given
        if n is None:
            return SemistandardMultisetTableaux_all(order)

        # weight not given
        if weight is None:
            raise NotImplementedError('class of all semistandard multiset '
                                      'tableaux of given weight not yet '
                                      'implemented')

        # shape given
        elif n in _Partitions:
            return SemistandardMultisetTableaux_shape_weight(_Partitions(n), weight, order)

        # size given
        if not isinstance(n, (int, Integer)) or n < 0:
            raise ValueError(
                "the argument must be a non-negative integer or a partition")

        return SemistandardMultisetTableaux_size_weight(n, weight, order)

    Element = SemistandardMultisetTableau
    def __init__(self, order='last_letter', key=None, **kwds):
        """
        Initialize ``self`` given an order.

        EXAMPLES::

            sage: S = SemistandardMultisetTableaux()
            sage: TestSuite(S).run()
        """
        if callable(order):
            self.order = order

            def order_to_key(order_fn):
                class K(object):
                    __slots__ = ['obj']
                    def __init__(self, obj, *args):
                        self.obj = obj
                    def __eq__(self, other):
                        return self.obj == other.obj
                    def __ne__(self, other):
                        return not self == other
                    def __lt__(self, other):
                        return order_fn(x, y) and not self == other
                    def __gt__(self, other):
                        return (not self <= other)
                    def __le__(self, other):
                        return self < other or self == other
                    def __ge__(self, other):
                        return not self < other
                    def __hash__(self):
                        raise TypeError('hash not implemented')
                return K

            self.key = order_to_key(order)

        elif order == 'last_letter':
            self.order = 'last_letter'
            self.key = self.last_letter_key

        elif order == 'grlex':
            self.order = 'grlex'
            self.key = self.grlex_key

        else:
            raise ValueError("An order should be given")

        if 'max_entry' in kwds:
            self.max_entry = kwds['max_entry']
            kwds.pop('max_entry')
        else:
            self.max_entry = None
        Tableaux.__init__(self, **kwds)

    @staticmethod
    def last_letter_key(x):
        return -1 if not x else max(x)

    @staticmethod
    def grlex_key(x):
        return (len(x), sorted(x))

    def __contains__(self,x):
        """
        Checks if x is an element of self.

        TESTS::

            sage: T = sage.combinat.tableau.SemistandardMultisetTableaux_all()
            sage: [[[1,2],[1,2]],[[2]]] in T
            True
            sage: [] in T
            True
            sage: Tableau([[[1]]]) in T
            True
            sage: StandardMultisetTableau([[[1]]]) in T
            True

            sage: [[[1],[2]],[[1]]] in T
            False
            sage: [[[1],[1]],[[5]]] in T
            True
            sage: [[[1],[3],[2]]] in T
            False
        """
        # TODO self.order == x.order will return False unless they are a reference
        # to exactly the same order function (i.e. (lambda x: x) != (lambda x: x) since
        # they are different references to 'the same function definition')
        if isinstance(x,SemistandardMultisetTableau) and x.order==self.order:
            return True
        elif isinstance(x,Tableau):
            order = self.order
            for row in x:
                if not all(all(c>0 for c in C) for C in row):
                    return False
                elif not all(order(row[i],row[i+1]) for i in range(len(row)-1)):
                    return False
                elif not all(order):
                    pass

class SemistandardMultisetTableaux_all(Tableaux):

    def __init__(self, order):
        pass

    def _repr_(self):
        pass

    def an_element(self):
        pass

class SemistandardMultisetTableaux_size_weight(Tableaux):
    def __init__(self, size, weight, order):
        pass

    def _repr_(self):
        pass

class SemistandardMultisetTableaux_shape_weight(SemistandardMultisetTableaux):

    def __init__(self, shape, weight, order='last_letter', **kwds):
        super(SemistandardMultisetTableaux_shape_weight, self).__init__(order=order, **kwds)
        self.shape = shape
        self.weight = weight

    def _repr_(self):
        """
        TESTS::

            sage: SemistandardMultisetTableaux([2,1,1], max_entry=4)
            Semistandard multiset tableaux of shape [2, 1, 1] and max entry 4 with order last_letter
        """
        return ("Semistandard multiset tableaux of shape {0} and max entry {1} "
                "with order {2}".format(self.shape, self.max_entry, self.order))

    def __iter__(self):
        def weight_part_to_mset_maxlen(vp, N):
            mset = []

            if len(vp) > N:
                # Weight partition is too long
                return None

            for part in vp:
                mset_part = [i+1 for i,n in enumerate(part) for _ in range(n)]
                mset.append(mset_part)

            # Pad by empty sets
            mset.extend([] for _ in range(N-len(vp)))
            return mset

        # MSET is list of multiset-partitions of content to fill shape of tableau
        # a typical element of MSET will be list of lists (some repeated and some empty) 
        # with total content equal to self.content
        weight_vec_p = VectorPartitions(self.weight)
        N = sum(self.shape)
        mset_partitions = [
            sorted(wp, key=self.key, reverse=True)
            for wp in (weight_part_to_mset_maxlen(vp, N) for vp in weight_vec_p)
            if wp is not None
        ]

        # for each multiset partition, group into equivalence classes and totally order
        def order_equivalent(x, y):
            kx, ky = self.key(x), self.key(y)
            return kx <= ky and ky <= kx

        for mset_part in mset_partitions:
            tab_entries = defaultdict(list)
            current_val = 1

            tab_entries[current_val].append(mset_part[0])

            for current_p, next_p in zip(mset_part, mset_part[1:]):
                if order_equivalent(current_p, next_p):
                    tab_entries[current_val].append(next_p)
                else:
                    current_val += 1
                    tab_entries[current_val].append(next_p)

        # run over SSYT and replace i's with multiset of things in ith equiv class
        wt = [len(eqclass) for eqclass in tab_entries.values()]
        print(tab_entries)
        print(wt)
        print()
        for t in SemistandardTableaux(self.shape, weight=wt):
            mt = [list(r) for r in t]
            print("sst:", mt)
            for i, eq_class in tab_entries.items():
                n_icells = len(t.cells_containing(i))
                print("eqc:", eq_class, "  nic:", n_icells)
                for assignment in Arrangements([tuple(x) for x in eq_class], n_icells):
                    print("assg:", assignment)
                    for (row, col), assign_set in zip(t.cells_containing(i), assignment):
                        mt[row][col] = tuple(assign_set)
            yield mt#self.element_class(self, mt, order=self.order, key=self.key)

        return
