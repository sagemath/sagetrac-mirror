r"""
SkewTableax, SemistandardSkewTableaux, StandardSkewTableaux, and related
Parent classes.

See :class:`SkewTableau` and its subclasses for the corresponding Element
classes.

AUTHORS:

- Mike Hansen: Initial version
- Travis Scrimshaw, Arthur Lubovsky (2013-02-11): Factored out
  ``CombinatorialClass``
- Josh Swanson (and others) (2015): tableaux refactoring/cleanup
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
#                     2013 Travis Scrimshaw <tscrim at ucdavis.edu>
#                     2013 Arthur Lubovsky
#                     2015 Jan Keitel
#                     2015 Josh Swanson
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
#*****************************************************************************

import six
from itertools                                 import ifilter
from sage.categories.infinite_enumerated_sets  import InfiniteEnumeratedSets
from sage.categories.finite_enumerated_sets    import FiniteEnumeratedSets
from sage.combinat.integer_vector              import IntegerVectors
from sage.combinat.partition                   import Partitions
from sage.combinat.skew_partition              import SkewPartition, SkewPartitions
from sage.matrix.all                           import zero_matrix
from sage.rings.all                            import Integer, ZZ, QQ, NN
from sage.rings.arith                          import factorial
from sage.structure.parent                     import Parent

from sage.combinat.tableaux.tableaux_options   import TableauOptions
from sage.combinat.tableaux.bad_shape_tableaux import BadShapeTableaux
from sage.combinat.tableaux.skew_tableau       import SkewTableau
from sage.combinat.tableaux.skew_tableau       import SkewTableauFactory
from sage.combinat.tableaux.skew_tableau       import SemistandardSkewTableau
from sage.combinat.tableaux.skew_tableau       import SemistandardSkewTableauFactory
from sage.combinat.tableaux.skew_tableau       import StandardSkewTableau
from sage.combinat.tableaux.skew_tableau       import StandardSkewTableauFactory

class SkewTableaux(BadShapeTableaux):
    r"""
    Parent class of all skew tableaux.

    See meth:`_element_constructor_` for Element class construction
    options and further examples. See :class:`SkewTableau` for the Element
    class.

    EXAMPLES::

        sage: SkewTableau([[None, 1, 2], [None, 3], [1]]).shape()
        [3, 2, 1] / [1, 1]

    TESTS::

        sage: S = SkewTableaux()
        sage: TestSuite(S).run()
    """
    Element = SkewTableau
    global_options = TableauOptions

    def _element_constructor_(self, st=0, expr=None, shape_word=None,
                                    dct=None, check=True):
        if st is not 0:
            return self.from_st(st, check)

        if expr is not None:
            return self.from_expr(expr, check)

        if shape_word is not None:
            shape, word = shape_word
            return self.from_shape_and_word(shape, word, check)

        if dct is not None:
            return self.from_dict(dct, check)

        return self._new_element(())
    # For user convenience, documentation for this method has
    #    been placed on the skew tableau factory function.
    _element_constructor_.__doc__ = SkewTableauFactory.__doc__

    def _an_element_(self):
        r"""
        Return a typical element of ``self``.

        TESTS::

            sage: SkewTableaux().an_element() # indirect doctest
            [[None, None, None, 3, 1, 2], [None, 2, 3, 3, 0], [1, -1]]
        """
        return self([[None, None, None, 3, 1, 2],
                     [None, 2, 3, 3, 0],
                     [1, -1]])

    def _coerce_map_from_(self, S):
        r"""
        Determine whether or not we can coerce from `S` to the current
        Parent.

        Subclasses are not (currently) allowed to coerce by default, but they
        should typically just work. Skew tableaux are often thought of as lists
        of lists, so allow that format as well.

        TESTS::

            sage: SkewTableaux().has_coerce_map_from(SemistandardSkewTableaux())
            True
            sage: SkewTableaux().has_coerce_map_from(list)
            True
        """
        if isinstance(S, SkewTableaux):
            return True
        elif S in (list,tuple):
            return True
        return False

    def __contains__(self, other):
        r"""
        Determine whether or not ``other`` is contained in ``self``.
        
        The default implementation from :class:`Parent` is too restrictive
        and requires `other == self(other)`. This does not work for
        skew tableaux, which are frequently thought of as lists of lists.

        TESTS::

            sage: SkewTableau([[None, 1, 2], [3]]) in SkewTableaux()
            True
            sage: [[None, 1, 2], [3]] in SkewTableaux()
            True
            sage: [[None, 1, 2], [3, None]] in SkewTableaux()
            False
            sage: StandardSkewTableau([[None, 1, 2], [None, 3, 4]]) in SkewTableaux()
            True
            sage: [[None, 2], [1, 3]] in SemistandardSkewTableaux()
            True
            sage: [[None, 2], [1, 3]] in StandardSkewTableaux()
            True
            sage: [[None, 2], [2, 4]] in StandardSkewTableaux()
            False
        """
        if isinstance(other, self.Element):
            return True
        try:
            self(other)
        except Exception:
            return False
        return True

    def from_st(self, st, check=True):
        r"""
        Construct a skew tableaux, optionally validating input.

        INPUT:
        - ``st``    -- an iterable of rows from top to bottom in English
          notation, where each row is an iterable of entries from left
          to right but where ``None``'s indicate the cells of the inner
          shape
        - ``check`` -- boolean, ``True`` if input should be validated

        EXAMPLE::

            sage: SkewTableaux().from_st([[None, 1, 2], [3, 4]])
            [[None, 1, 2], [3, 4]]

        TESTS::

            sage: SkewTableaux().from_st([[None, 1, None], [3, 4]])
            Traceback (most recent call last):
            ...
            ValueError: Row (None, 1, None) has a None after a non-None
            sage: SkewTableaux().from_st(([[None, 1, 2], [3, 4, 5, 6]]))
            Traceback (most recent call last):
            ...
            ValueError: Input must be of skew shape
        """
        try:
            # Remove empty rows, normalize input
            st = tuple( ifilter(None,
                                (tuple(row) for row in st)) )
        except TypeError:
            raise TypeError("each element of the skew tableau must be an iterable")

        if check:
            # Make sure rows begin with blocks of nones
            inner = [0]*len(st)
            outer = [0]*len(st)
            for i in range(len(st)):
                row = st[i]

                finishedNones = False
                noneCount = 0
                for val in row:
                    if val is None:
                        if finishedNones:
                            raise ValueError("Row %s has a None after a non-None"%str(row))
                        noneCount += 1
                    else:
                        finishedNones = True

                inner[i] = noneCount
                outer[i] = len(row)

            # Make sure inner and outer form a skew partition
            if (outer, inner) not in SkewPartitions():
                raise ValueError("Input must be of skew shape")

        return self._new_element(st)

    def from_expr(self, expr, check=True):
        r"""
        Return a :class:`SkewTableau` from a MuPAD-Combinat ``expr``
        for a skew tableau.

        Provided primarily for compatibility with MuPAD-Combinat.

        INPUT:

        - ``expr`` -- a pair (``inner``, ``rows``) where ``inner`` is
          the inner partition and ``rows`` is an iterable of rows from
          bottom to top in English notation, where each row is an iterable
          of the entries in that row from left to right.
        - ``check`` -- boolean, ``True`` if input should be validated

        EXAMPLES::

            sage: SkewTableaux().from_expr([[1,1], [[5],[3,4],[1,2]]])
            [[None, 1, 2], [None, 3, 4], [5]]
            sage: SkewTableaux().from_expr([[1,1,1], [[5],[3,4]]])
            [[None, 3, 4], [None, 5], [None]]

        TESTS::

            sage: SkewTableaux().from_expr([[1,2,1], [[5],[3,4]]])
            Traceback (most recent call last):
            ...
            ValueError: shape must form a skew partition
        """
        inner, rows = expr

        # Make inner and rows the same length by padding
        li = len(inner)
        lr = len(rows)
        if lr > li:
            inner += [0]*(lr-li)
            li=lr
        elif lr < li:
            rows = [[]]*(li-lr) + rows
            lr=li

        if check:
            # Make sure shape is skew
            outer = [inner[i] + len(rows[-(i+1)]) for i in range(lr)]
            if (outer, inner) not in SkewPartitions():
                raise ValueError("shape must form a skew partition")

        # Create a normalized st representation
        st = tuple( (None,)*(inner[i]) + tuple(rows[-(i+1)])
                    for i in range(lr) )

        return self._new_element(st)

    def from_shape_and_word(self, shape, word, check=True):
        r"""
        Return the skew tableau corresponding to the skew partition ``shape``
        and the word ``word`` obtained from the row reading.

        INPUT:

        - ``shape`` -- a skew partition
        - ``word`` -- the row reading word
        - ``check`` -- boolean, ``True`` if input should be validated

        EXAMPLES::

            sage: t = SkewTableau([[None, 1, 3], [None, 2], [4]])
            sage: shape = t.shape(); shape
            [3, 2, 1] / [1, 1]
            sage: word  = t.to_word(); word
            word: 4213
            sage: SkewTableaux().from_shape_and_word(shape, word)
            [[None, 1, 3], [None, 2], [4]]

        TESTS::

            sage: SkewTableaux().from_shape_and_word([[3, 2], [1]], [1, 2, 3])
            Traceback (most recent call last):
            ...
            ValueError: Shape and word must be the same length
            sage: SkewTableaux().from_shape_and_word([[2, 3], [1]], [1, 2, 3, 3])
            Traceback (most recent call last):
            ...
            ValueError: Input must be of skew shape
        """
        outer, inner = shape

        if check:
            if sum(outer) - sum(inner) != len(word):
                raise ValueError("Shape and word must be the same length")

        # Construct the default input format
        st = [ [None]*row_length for row_length in outer ]
        w_count = 0
        for i in reversed(range(len(outer))):
            for j in range(outer[i]):
                if i >= len(inner) or j >= inner[i]:
                    st[i][j] = word[w_count]
                    w_count += 1

        # from_st performs additional validation and normalization
        return self.from_st(st, check)

    def from_dict(self, dct, check=True):
        r"""
        Return the skew tableau obtained from the given dictionary.

        The shape of the skew tableau is not uniquely specified by the
        dictionary in general, so the smallest outer shape is used.

        INPUT:

        - ``dct`` -- a dictionary whose keys are pairs of non-negative
          integers
        - ``check`` -- boolean, ``True`` if input should be validated

        EXAMPLES::

            sage: SkewTableaux().from_dict({(0, 1): 1, (0, 2): 2, (1, 1): 2})
            [[None, 1, 2], [None, 2]]
            sage: s1 = SkewTableau(dct={(0, 1): 1}); s1
            [[None, 1]]
            sage: s2 = SkewTableau([[None, 1], [None]]); s2
            [[None, 1], [None]]
            sage: s1.dict() == s2.dict()
            True
            sage: s1.shape() == s2.shape()
            False

        TESTS::

            sage: SkewTableaux().from_dict({'cat': 1, (0, 2): 3})
            Traceback (most recent call last):
            ...
            ValueError: Keys must be pairs of non-negative integers
            sage: SkewTableaux().from_dict({(1, 2): 3, (2, 1): 4, (15, 3): 5})
            Traceback (most recent call last):
            ...
            ValueError: Outer shape must form a partition
        """
        # Can assume non-emptiness
        if len(dct)==0:
            return self._new_element(())

        if check:
            # Make sure indexes are pairs of non-negative integers
            try:
                for v in six.iterkeys(dct):
                    if v[0] not in NN or v[1] not in NN:
                        raise ValueError
            except (TypeError, ValueError):
                raise ValueError("Keys must be pairs of non-negative integers")

        # Figure out the outer shape
        max_row_index = max(six.iterkeys(dct), key=lambda k: k[0])[0]
        outer = [0]*(max_row_index+1)
        for r, c in six.iterkeys(dct):
            if outer[r] < c+1: outer[r] = c+1
        for i in range(1,len(outer)+1):
            if outer[-i] == 0:
                outer[-i] = c
            else:
                c = outer[-i]

        if check:
            if outer not in Partitions():
                raise ValueError("Outer shape must form a partition")

        # Compute the st representation
        st = [[None]*l for l in outer]
        for k, v in six.iteritems(dct):
            st[k[0]][k[1]] = v

        # Additional normalization and validation is performed in from_st
        return self.from_st(st, check)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: SkewTableaux()
            Skew tableaux
        """
        return "Skew tableaux"

class SemistandardSkewTableaux(SkewTableaux):
    r"""
    Parent class of all semistandard skew tableaux, possibly with a given
    upper bound on entries.

    See meth:`_element_constructor_` for Element class construction
    options and further examples. See :class:`SemistandardSkewTableau`
    for the Element class itself.

    EXAMPLES::

        sage: SemistandardSkewTableau([[None, 1, 2], [None, 3], [1]]).shape()
        [3, 2, 1] / [1, 1]
        sage: [[None, 1, 2], [None, 3], [1]] in SemistandardSkewTableaux(max_entry=4)
        True
        sage: [[None, 1, 2], [None, 3], [1]] in SemistandardSkewTableaux(max_entry=2)
        False
    """
    Element = SemistandardSkewTableau

    def __init__(self, max_entry=None):
        r"""
        Initialize ``self``.

        INPUT:
        
        - ``max_entry`` -- an upper bound on entries.

        TESTS::

            sage: S = SemistandardSkewTableaux()
            sage: TestSuite(S).run()

            sage: S = SemistandardSkewTableaux(max_entry=3)
            sage: TestSuite(S).run()
        """
        Parent.__init__(self, category=InfiniteEnumeratedSets())
        if max_entry is not None:
            self.max_entry = max_entry

    def _element_constructor_(self, *args, **kwds):
        ret = super(SemistandardSkewTableaux, self)._element_constructor_(*args, **kwds)

        if kwds.get("check", True):
            if not ret.is_semistandard():
                raise ValueError("Input is not semistandard")

            if hasattr(self, "max_entry") and ret.size() > 0:
                if max(ret.iter_entries()) > self.max_entry:
                    raise ValueError("Entries must be at most %i"%self.max_entry)

        return ret
    # For user convenience, documentation for this method has
    #    been placed on the semistandard skew tableau factory function.
    _element_constructor_.__doc__ = SemistandardSkewTableauFactory.__doc__

    def _an_element_(self):
        r"""
        Return a typical element of ``self``.

        TESTS::

            sage: SemistandardSkewTableaux().an_element() in SemistandardSkewTableaux()
            True
        """
        # Use Parent's default implementation, since we have an iterator available
        return super(Parent, self)._an_element_()

    def _repr_(self):
        r"""
        Return the string representation of ``self``.

        TESTS::

            sage: SemistandardSkewTableaux()
            Semistandard skew tableaux
            sage: SemistandardSkewTableaux(max_entry=5)
            Semistandard skew tableaux with maximum entry 5
        """
        if not hasattr(self, "max_entry"):
            return "Semistandard skew tableaux"
        return "Semistandard skew tableaux with maximum entry {}".format(self.max_entry)

    def __iter__(self):
        r"""
        Iterate over the tableux in ``self``.
        
        We restrict to reduced skew shapes in the sense of
        :class:`SkewPartitions`.

        EXAMPLES::

            sage: it = iter(SemistandardSkewTableaux(max_entry=5))
            sage: [next(it) for x in range(12)]
            [[],
             [[1]],
             [[2]],
             [[3]],
             [[4]],
             [[5]],
             [[1, 1]],
             [[1, 2]],
             [[1, 3]],
             [[1, 4]],
             [[1, 5]],
             [[2, 2]]]

        If no max entry is specified, the iteration goes over all
        semistandard skew tableaux of size `n` with max entry `n`,
        for all `n`::

            sage: it = iter(SemistandardSkewTableaux())
            sage: [next(it) for x in range(10)]
            [[],
             [[1]],
             [[1, 1]],
             [[1, 2]],
             [[2, 2]],
             [[1], [2]],
             [[None, 1], [1]],
             [[None, 2], [1]],
             [[None, 1], [2]],
             [[None, 2], [2]]]
        """
        if not hasattr(self, "max_entry"):
            n = 0
            while True:
                for ssst in SemistandardSkewTableaux_size(n, n):
                    yield self(ssst, check=False)
                n += 1
        else:
            n = 0
            while True:
                for ssst in SemistandardSkewTableaux_size(n, self.max_entry):
                    yield self(ssst, check=False)
                n += 1

class SemistandardSkewTableaux_size(SemistandardSkewTableaux):
    r"""
    Parent class of all semistandard skew tableaux of a fixed size `n`,
    possibly with an upper bound on entries which defaults to `n`.

    EXAMPLES::

        sage: [[None, 1, 2, 3], [1, 2]] in SemistandardSkewTableaux(5)
        True
        sage: [[None, 1, 2, 3], [1, 2]] in SemistandardSkewTableaux(4)
        False
    """
    def __init__(self, n, max_entry=None):
        r"""
        Initialize ``self``.

        INPUT:
        - ``n``         -- the number of cells in these tableaux
        - ``max_entry`` -- an upper bound on the entries in these
          tableaux

        EXAMPLES::

            sage: S = SemistandardSkewTableaux(3)
            sage: TestSuite(S).run()
        """
        Parent.__init__(self, category=FiniteEnumeratedSets())
        self.n = n
        if max_entry is None:
            self.max_entry = n
        else:
            self.max_entry = max_entry

    def _element_constructor_(self, *args, **kwds):
        r"""
        Return a new Element of ``self``.
        
        See :class:`SkewTableaux` for further details.

        TESTS::

            sage: S = SemistandardSkewTableaux(p=4, max_entry=3)
            sage: [[None, 1, 3], [2, 3]] in S
            True
            sage: [[None, 1, 3], [2, 4]] in S
            False
            sage: [[None, 1, 3], [2]] in S
            False
        """
        ret = super(SemistandardSkewTableaux_size, self)._element_constructor_(*args, **kwds)

        if kwds.get("check", True):
            s=ret.size()
            if s != self.n:
                raise ValueError("Input must have size %i"%self.n)
            if s > 0 and max(ret.iter_entries()) > self.max_entry:
                raise ValueError("Entries must be at most %i"%self.max_entry)

        return ret

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: SemistandardSkewTableaux(3)
            Semistandard skew tableaux of size 3 and maximum entry 3
            sage: SemistandardSkewTableaux(3, max_entry=8)
            Semistandard skew tableaux of size 3 and maximum entry 8
        """
        return ("Semistandard skew tableaux of size %s and maximum entry %s"
                %(repr(self.n), repr(self.max_entry)))

    def cardinality(self):
        r"""
        Return the cardinality of ``self``.

        EXAMPLES::

            sage: SemistandardSkewTableaux(2).cardinality()
            8
        """
        return sum(SemistandardSkewTableaux_shape(p, self.max_entry).cardinality()
                   for p in SkewPartitions(self.n))

    def __iter__(self):
        r"""
        Iterate over the tableaux in ``self``.

        EXAMPLES::

            sage: SemistandardSkewTableaux(2).list()
            [[[1, 1]],
             [[1, 2]],
             [[2, 2]],
             [[1], [2]],
             [[None, 1], [1]],
             [[None, 2], [1]],
             [[None, 1], [2]],
             [[None, 2], [2]]]
        """
        for p in SkewPartitions(self.n):
            for ssst in SemistandardSkewTableaux_shape(p, self.max_entry):
                yield self(ssst, check=False)

class SemistandardSkewTableaux_size_weight(SemistandardSkewTableaux):
    r"""
    Parent class of all semistandard skew tableaux of a fixed size and
    weight.

    To ensure finiteness, we restrict to reduced skew shapes in the
    sense of :class:`SkewPartitions`.

    EXAMPLES::

        sage: [[None, 1, 2, 3], [1, 2]] in SemistandardSkewTableaux(5, (2, 2, 1))
        True
        sage: [[None, 1, 2, 3], [1, 2]] in SemistandardSkewTableaux(5, (2, 1, 2))
        False
    """
    def __init__(self, n, weight):
        r"""
        Initialize ``self``.

        INPUT:
        - ``n``      -- the number of cells of these tableaux
        - ``weight`` -- a tuple of integers giving the weight of these tableaux

        TESTS::

            sage: S = SemistandardSkewTableaux(3,[2,1])
            sage: TestSuite(S).run()
        """
        Parent.__init__(self, category=FiniteEnumeratedSets())
        self.n      = n
        self.weight = weight

    def _element_constructor_(self, *args, **kwds):
        r"""
        Return a new Element of ``self``.
        
        See :class:`SkewTableaux` for further details.

        TESTS::

            sage: S = SemistandardSkewTableaux(p=4, max_entry=3)
            sage: [[None, 1, 3], [2, 3]] in S
            True
            sage: [[None, 1, 3], [2, 4]] in S
            False
            sage: [[None, 1, 3], [2]] in S
            False
        """
        ret = super(SemistandardSkewTableaux_size_weight, self)._element_constructor_(*args, **kwds)

        if kwds.get("check", True):
            if ret.size() != self.n:
                raise ValueError("Input must have size %i"%self.n)
            if tuple(ret.weight()) != self.weight:
                raise ValueError("Input must have weight %s"%repr(self.weight))

        return ret

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: SemistandardSkewTableaux(3,[2,1])
            Semistandard skew tableaux of size 3 and weight [2, 1]
        """
        return ("Semistandard skew tableaux of size %s and weight %s"
                %(repr(self.n), list(self.weight)))

    def cardinality(self):
        r"""
        Return the cardinality of ``self``.

        EXAMPLES::

            sage: SemistandardSkewTableaux(2,[1,1]).cardinality()
            4
        """
        return sum(SemistandardSkewTableaux_shape_weight(p, self.weight).cardinality()
                   for p in SkewPartitions(self.n))

    def __iter__(self):
        r"""
        Iterate over the tableaux in ``self``.

        We restrict to reduced skew shapes in the
        sense of :class:`SkewPartitions`.

        EXAMPLES::

            sage: list(SemistandardSkewTableaux(2,[1,1]))
            [[[1, 2]], [[1], [2]], [[None, 2], [1]], [[None, 1], [2]]]
        """
        for p in SkewPartitions(self.n):
            for ssst in SemistandardSkewTableaux_shape_weight(p, self.weight):
                yield self(ssst, check=False)

class SemistandardSkewTableaux_shape(SemistandardSkewTableaux):
    r"""
    Parent class of semistandard skew tableaux of a fixed skew shape
    with a given upper bound on entries.

    EXAMPLES::

        sage: [[None, 1, 2, 3], [1, 2]] in SemistandardSkewTableaux([[4,2], [1]])
        True
        sage: [[None, 1, 2, 3], [1, 2]] in SemistandardSkewTableaux([[4,2], [2]])
        False
    """
    def __init__(self, p, max_entry=None):
        r"""
        Initialize ``self``.

        INPUT:

        - ``p``         -- a skew partition
        - ``max_entry`` -- an upper bound on entries; defaults to the
          size of ``p``.

        TESTS::

            sage: S = SemistandardSkewTableaux([[2,1],[]])
            sage: S == loads(dumps(S))
            True
            sage: TestSuite(S).run()
        """
        Parent.__init__(self, category=FiniteEnumeratedSets())

        if max_entry is None:
            max_entry = sum(p[0])-sum(p[1])
        self.p = p
        self.max_entry = max_entry

    def _element_constructor_(self, *args, **kwds):
        r"""
        Return a new Element of ``self``.
        
        See :class:`SkewTableaux` for further details.

        TESTS::

            sage: S = SemistandardSkewTableaux(p=[[3, 1], [1]])
            sage: [[None, 1, 3], [2]] in S
            True
            sage: [[None, 1, 3], [2, 4]] in S
            False
            sage: [[None, 1, 4], [2]] in S
            False
        """
        ret = super(SemistandardSkewTableaux_shape, self)._element_constructor_(*args, **kwds)

        if kwds.get("check", True):
            if ret.shape() != self.p:
                raise ValueError("Input must have shape %i"%repr(self.p))
            if ret.size() > 0 and max(ret.iter_entries()) > self.max_entry:
                raise ValueError("Entries must be at most %i"%self.max_entry)

        return ret

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: SemistandardSkewTableaux([[2,1],[]])
            Semistandard skew tableaux of shape [2, 1] / [] and maximum entry 3
        """
        return ("Semistandard skew tableaux of shape %s and maximum entry %s"
                %(repr(self.p), repr(self.max_entry)))

    def cardinality(self):
        r"""
        Return the cardinality of ``self``.

        EXAMPLES::

            sage: SemistandardSkewTableaux([[2,1],[]]).cardinality()
            8
            sage: SemistandardSkewTableaux([[2,1],[]], max_entry=2).cardinality()
            2
        """
        return sum(SemistandardSkewTableaux_shape_weight(self.p, tuple(weight)).cardinality()
                   for weight in IntegerVectors(self.p.size(), self.max_entry))

    def __iter__(self):
        r"""
        Iterate over the tableaux in ``self``.

        EXAMPLES::

            sage: SemistandardSkewTableaux([[2,1],[]]).list()
            [[[1, 1], [2]],
             [[1, 1], [3]],
             [[1, 2], [2]],
             [[1, 3], [2]],
             [[1, 2], [3]],
             [[1, 3], [3]],
             [[2, 2], [3]],
             [[2, 3], [3]]]
            sage: from sage.combinat.skew_tableau import SemistandardSkewTableaux_shape
            sage: SemistandardSkewTableaux_shape([[2,1],[]], max_entry=2).list()
            [[[1, 1], [2]], [[1, 2], [2]]]
        """
        for weight in IntegerVectors(self.p.size(), self.max_entry):
            for ssst in SemistandardSkewTableaux_shape_weight(self.p, tuple(weight)):
                yield self(ssst, check=False)

class SemistandardSkewTableaux_shape_weight(SemistandardSkewTableaux):
    r"""
    Parent class of semistandard skew tableaux of a fixed skew shape
    and weight.

    EXAMPLES::

        sage: [[None, 1, 2, 3], [1, 2]] in SemistandardSkewTableaux(p=[[4,2], [1]], weight=(2, 2, 1))
        True
        sage: [[None, 1, 2, 3], [1, 1]] in SemistandardSkewTableaux(p=[[4,2], [1]], weight=(2, 2, 1))
        False
    """
    def __init__(self, p, weight):
        r"""
        INPUT:

        - ``p``      -- a skew partition giving the shape of these tableaux
        - ``weight`` -- a tuple of integers giving the weight of these tableaux

        TESTS::

            sage: S = SemistandardSkewTableaux([[2,1],[]],[2,1])
            sage: S == loads(dumps(S))
            True
            sage: TestSuite(S).run()
        """
        Parent.__init__(self, category=FiniteEnumeratedSets())
        self.p      = p
        self.weight = weight


    def _element_constructor_(self, *args, **kwds):
        r"""
        Return a new Element of ``self``.
        
        See :class:`SkewTableaux` for further details.

        TESTS::

            sage: S = SemistandardSkewTableaux(p=[[3, 1], [1]], weight=(2, 1))
            sage: [[None, 1, 2], [1]] in S
            True
            sage: [[None, 1, 3], [2]] in S
            False
            sage: [[1, 1], [2]] in S
            False
        """
        ret = super(SemistandardSkewTableaux_shape_weight, self)._element_constructor_(*args, **kwds)

        if kwds.get("check", True):
            if ret.shape() != self.p:
                raise ValueError("Input must have shape %i"%repr(self.p))
            if tuple(ret.weight()) != self.weight:
                raise ValueError("Input must have weight %s"%repr(self.weight))

        return ret

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: SemistandardSkewTableaux([[2,1],[]],[2,1])
            Semistandard skew tableaux of shape [2, 1] / [] and weight [2, 1]
        """
        return ( "Semistandard skew tableaux of shape %s and weight %s"
                 %(repr(self.p), list(self.weight)) )

    def __iter__(self):
        r"""
        Iterate over ``self``.

        EXAMPLES::

            sage: list(SemistandardSkewTableaux([[2,1],[]],[2,1]))
            [[[1, 1], [2]]]
        """
        from sage.combinat.ribbon_tableau import RibbonTableaux_shape_weight_length
        for x in RibbonTableaux_shape_weight_length(self.p, tuple(self.weight), 1):
            yield self(x, check=False)

class StandardSkewTableaux(SemistandardSkewTableaux):
    r"""
    Parent class of standard skew tableaux.

    EXAMPLES::

        sage: [[None, 1, 2, 4], [3, 5]] in StandardSkewTableaux()
        True
        sage: [[None, 1, 2, 3], [1, 1]] in StandardSkewTableaux()
        False
    """
    Element = StandardSkewTableau

    def _element_constructor_(self, *args, **kwds):
        ret = super(StandardSkewTableaux, self)._element_constructor_(*args, **kwds)

        if kwds.get("check", True):
            # Semistandardness has already been checked
            if not ret.has_standard_entries():
                raise ValueError("Input is not standard")

        return ret
    # For user convenience, documentation for this method has
    #    been placed on the standard skew tableau factory function.
    _element_constructor_.__doc__ = StandardSkewTableauFactory.__doc__

    def __init__(self):
        r"""
        Initialize ``self``.

        TESTS::

            sage: S = StandardSkewTableaux()
            sage: TestSuite(S).run()
        """
        Parent.__init__(self, category=InfiniteEnumeratedSets())

    def _repr_(self):
        r"""
        Return a string representatino of ``self``.

        TESTS::

            sage: StandardSkewTableaux()
            Standard skew tableaux
        """
        return "Standard skew tableaux"

    def __iter__(self):
        r"""
        Iterate through the tableaux in ``self``.

        We restrict to reduced skew shapes in the
        sense of :class:`SkewPartitions`.

        EXAMPLES::

            sage: it = StandardSkewTableaux().__iter__()
            sage: [next(it) for x in range(10)]
            [[],
             [[1]],
             [[1, 2]], [[1], [2]], [[None, 1], [2]], [[None, 2], [1]],
             [[1, 2, 3]], [[1, 2], [3]], [[1, 3], [2]],
             [[None, 1, 2], [3]]]
        """
        n = 0
        while True:
            for st in StandardSkewTableaux_size(n):
                yield self(st, check=False)
            n += 1

class StandardSkewTableaux_size(StandardSkewTableaux):
    r"""
    Parent class of standard skew tableaux of a fixed size.

    EXAMPLES::

        sage: [[None, 1, 2, 4], [3, 5]] in StandardSkewTableaux(5)
        True
        sage: [[None, 1, 2, 3], [4]] in StandardSkewTableaux(5)
        False
    """
    def __init__(self, n):
        r"""
        Initialize ``self``.

        INPUT:
        - ``n`` -- the number of cells in these tableaux

        TESTS::

            sage: S = StandardSkewTableaux(3)
            sage: TestSuite(S).run()
        """
        Parent.__init__(self, category=FiniteEnumeratedSets())
        self.n = n

    def _element_constructor_(self, *args, **kwds):
        r"""
        Return a new Element of ``self``.

        See :class:`SkewTableaux` for further details.

        TESTS::

            sage: S = StandardSkewTableaux(3)
            sage: [[None, 1, 2], [3]] in S
            True
            sage: [[None, 1, 2], [2]] in S
            False
            sage: [[1, 2], [2, 4]] in S
            False
        """
        ret = super(StandardSkewTableaux_size, self)._element_constructor_(*args, **kwds)

        if kwds.get("check", True):
            if ret.size() != self.n:
                raise ValueError("Input must have size %i"%repr(self.n))

        return ret

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: StandardSkewTableaux(3)
            Standard skew tableaux of size 3
        """
        return "Standard skew tableaux of size %s"%self.n

    def cardinality(self):
        r"""
        Return the cardinality of ``self``.

        We restrict to reduced skew shapes in the
        sense of :class:`SkewPartitions` (otherwise
        the count would be infinite).

        EXAMPLES::

            sage: StandardSkewTableaux(1).cardinality()
            1
            sage: StandardSkewTableaux(2).cardinality()
            4
            sage: StandardSkewTableaux(3).cardinality()
            24
            sage: StandardSkewTableaux(4).cardinality()
            194
        """
        return sum(StandardSkewTableaux_shape(skp).cardinality()
                   for skp in SkewPartitions(self.n))

    def __iter__(self):
        r"""
        Iterate through the tableaux in ``self``.

        We restrict to reduced skew shapes in the
        sense of :class:`SkewPartitions`.

        EXAMPLES::

            sage: StandardSkewTableaux(2).list()
            [[[1, 2]], [[1], [2]], [[None, 1], [2]], [[None, 2], [1]]]

            sage: StandardSkewTableaux(3).list()
            [[[1, 2, 3]],
             [[1, 2], [3]], [[1, 3], [2]],
             [[None, 1, 2], [3]], [[None, 1, 3], [2]],
             [[None, 2, 3], [1]],
             [[None, 1], [2, 3]], [[None, 2], [1, 3]],
             [[None, None, 1], [2, 3]], [[None, None, 2], [1, 3]], [[None, None, 3], [1, 2]],
             [[1], [2], [3]],
             [[None, 1], [None, 2], [3]], [[None, 1], [None, 3], [2]], [[None, 2], [None, 3], [1]],
             [[None, 1], [2], [3]], [[None, 2], [1], [3]], [[None, 3], [1], [2]],
             [[None, None, 1], [None, 2], [3]], [[None, None, 1], [None, 3], [2]],
             [[None, None, 2], [None, 1], [3]], [[None, None, 3], [None, 1], [2]],
             [[None, None, 2], [None, 3], [1]], [[None, None, 3], [None, 2], [1]]]
        """
        for skp in SkewPartitions(self.n):
            for sst in StandardSkewTableaux_shape(skp):
                yield self(sst, check=False)

class StandardSkewTableaux_shape(StandardSkewTableaux):
    r"""
    Parent class of standard skew tableaux of a fixed skew shape.

    EXAMPLES::

        sage: [[None, 1, 2, 4], [3, 5]] in StandardSkewTableaux([[4, 2], [1]])
        True
        sage: [[None, 1, 2, 3], [4]] in StandardSkewTableaux([[4, 2], [1]])
        False
    """
    def __init__(self, skp):
        r"""
        Initialize ``self``.

        INPUT:
        - ``skp`` -- a :class:`SkewPartition`

        TESTS::

            sage: S = StandardSkewTableaux([[3, 2, 1], [1, 1]])
            sage: TestSuite(S).run()
        """
        Parent.__init__(self, category=FiniteEnumeratedSets())
        self.skp = skp

    def _element_constructor_(self, *args, **kwds):
        r"""
        Return a new Element of ``self``.

        See :class:`SkewTableaux` for further details.

        TESTS::

            sage: S = StandardSkewTableaux([[3, 1], [1]])
            sage: [[None, 1, 2], [3]] in S
            True
            sage: [[None, 1, 2], [2]] in S
            False
            sage: [[1, 2], [2, 4]] in S
            False
        """
        ret = super(StandardSkewTableaux_shape, self)._element_constructor_(*args, **kwds)

        if kwds.get("check", True):
            if ret.shape() != self.skp:
                raise ValueError("Input must have shape %i"%repr(self.p))

        return ret

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: StandardSkewTableaux([[3, 2, 1], [1, 1]])
            Standard skew tableaux of shape [3, 2, 1] / [1, 1]
        """
        return "Standard skew tableaux of shape %s"%repr(self.skp)

    def cardinality(self):
        r"""
        Return the cardinality of ``self``.

        This uses a formula due to Aitken
        (see Cor. 7.16.3 of [Sta-EC2]_).

        EXAMPLES::

            sage: StandardSkewTableaux([[3, 2, 1], [1, 1]]).cardinality()
            8
        """
        outer, inner = self.skp
        m = len(outer)
        n = sum(outer) - sum(inner)
        outer = list(outer)
        inner = list(inner) + [0]*(m-len(inner))
        a = zero_matrix(QQ, m)
        for i in range(m):
            for j in range(m):
                v = outer[i] - inner[j] - i + j
                if v < 0:
                    a[i,j] = 0
                else:
                    a[i,j] = 1/factorial(v)
        return ZZ(factorial(n) * a.det())

    def __iter__(self):
        r"""
        Iterate over the tableaux in ``self``.

        The resulting standard skew tableaux are ordered
        lexicographically by the word obtained from their
        row reading.

        EXAMPLES::

            sage: [st for st in StandardSkewTableaux([[3, 2, 1], [1, 1]])]
            [[[None, 1, 2], [None, 3], [4]],
             [[None, 1, 2], [None, 4], [3]],
             [[None, 1, 3], [None, 2], [4]],
             [[None, 1, 4], [None, 2], [3]],
             [[None, 1, 3], [None, 4], [2]],
             [[None, 1, 4], [None, 3], [2]],
             [[None, 2, 3], [None, 4], [1]],
             [[None, 2, 4], [None, 3], [1]]]
        """
        dag = self.skp.to_dag(format="tuple")
        le_list = list(dag.topological_sort_generator())

        empty = [[None]*row_length for row_length in self.skp.outer()]

        for le in le_list:
            i = 1
            skew = [list(row) for row in empty]
            for row, column in le:
                    skew[row][column] = i
                    i += 1

            yield self(skew, check=False)

# Factory methods for constructing various Parent classes of tableaux
def SemistandardSkewTableauxFactory(p=None, weight=None, max_entry=None):
    r"""
    Factory function to create various Parent classes of certain types of
    semistandard skew tableaux.

    INPUT:
    - ``p`` -- either an integer specifying the number of cells in the
      tableaux, or a :class:`SkewPartition` specifying the skew shape,
      or `None` if unrestricted
    - ``weight`` -- an iterable specifying the weight of the tableaux
    - ``max_entry`` -- an upper bound for entries (the lower bound of 1
      being implicit), or `None` if unrestricted

    If neither the weight nor the maximum entry is specified, the
    maximum entry defaults to the size of the tableau. Note that
    "maximum entry" does not literally mean the highest entry; instead
    it is just an upper bound that no entry is allowed to surpass.

    EXAMPLES:

    The (infinite) class of all semistandard skew tableaux::

        sage: SemistandardSkewTableaux()
        Semistandard skew tableaux

    The (still infinite) class of all semistandard skew tableaux
    with maximum entry `2`::

        sage: SemistandardSkewTableaux(max_entry=2)
        Semistandard skew tableaux with maximum entry 2

    The class of all semistandard skew tableaux of given size `3`
    and maximum entry `3`::

        sage: SemistandardSkewTableaux(3)
        Semistandard skew tableaux of size 3 and maximum entry 3

    To set a different maximum entry::

        sage: SemistandardSkewTableaux(3, max_entry=7)
        Semistandard skew tableaux of size 3 and maximum entry 7

    Specifying a shape::

        sage: SemistandardSkewTableaux([[2,1],[]])
        Semistandard skew tableaux of shape [2, 1] / [] and maximum entry 3

    Specifying both a shape and a maximum entry::

        sage: S = SemistandardSkewTableaux([[2,1],[1]], max_entry = 3); S
        Semistandard skew tableaux of shape [2, 1] / [1] and maximum entry 3
        sage: S.list()[0:3]
        [[[None, 1], [1]], [[None, 2], [1]], [[None, 1], [2]]]

        sage: for n in range(5):
        ....:     print n, len(SemistandardSkewTableaux([[2,2,1],[1]], max_entry = n))
        0 0
        1 0
        2 1
        3 9
        4 35

    Specifying a shape and a weight::

        sage: SemistandardSkewTableaux([[2,1],[]],[2,1])
        Semistandard skew tableaux of shape [2, 1] / [] and weight [2, 1]

    (The maximum entry is redundant in this case and thus is ignored.)

    Specifying a size and a weight::

        sage: SemistandardSkewTableaux(3, [2,1])
        Semistandard skew tableaux of size 3 and weight [2, 1]

    TESTS::

        sage: S1 = SemistandardSkewTableaux(3, [2,1])
        sage: S2 = SemistandardSkewTableaux(int(3), (2,1))
        sage: S1 is S2
        True
        sage: S1 = StandardSkewTableaux([[3, 2, 1], [1, 1]])
        sage: S2 = StandardSkewTableaux(SkewPartition([[3, 2, 1], [1, 1]]))
        sage: S1 is S2
        True

    .. WARNING::

        If the shape is not specified, the iterator for the resulting
        Parents yields only those skew tableaux whose shape is reduced, in
        the sense of :class:`SkewPartitions`, meaning there are no empty rows
        before the last nonempty row, and there are no empty columns before
        the last nonempty column. (This is typically what one wants in
        applications.)
    """
    # Input validation and normalization for these Parent classes
    #   should be done here rather than in their `__init__` methods;
    #   analogously, input validation and normalization for Element
    #   classes should be done by their Parents.
    if p is None:
        if weight is None:
            return SemistandardSkewTableaux(max_entry)
        raise ValueError("With a weight, you must also specify a size or shape")

    if isinstance(p, (int, Integer)):
        if weight is None:
            return SemistandardSkewTableaux_size(p, max_entry)
        else:
            return SemistandardSkewTableaux_size_weight(p, tuple(weight))

    if p in SkewPartitions():
        if weight is None:
            return SemistandardSkewTableaux_shape(SkewPartition(p), max_entry)
        else:
            return SemistandardSkewTableaux_shape_weight(SkewPartition(p), tuple(weight))

    raise ValueError("Invalid input")

def StandardSkewTableauxFactory(skp=None):
    r"""
    Factory function to create various Parent classes of certain types of
    standard tableaux.

    INPUT:
    - ``skp` -- either an integer specifying the number of cells in the
      tableaux, or a :class:`SkewPartition` specifying the skew shape,
      or `None` if unrestricted

    EXAMPLES::

    The (infinite) class of all standard skew tableaux::

        sage: S = StandardSkewTableaux(); S
        Standard skew tableaux
        sage: S.cardinality()
        +Infinity

    The class of all semistandard skew tableaux of size `3`

        sage: S = StandardSkewTableaux(3); S
        Standard skew tableaux of size 3
        sage: S.cardinality()
        24

    Specifying a shape::

        sage: StandardSkewTableaux([[2,1],[]])
        Standard skew tableaux of shape [2, 1] / []
        sage: StandardSkewTableaux([[3, 2, 1], [1, 1]]).list()
        [[[None, 1, 2], [None, 3], [4]],
         [[None, 1, 2], [None, 4], [3]],
         [[None, 1, 3], [None, 2], [4]],
         [[None, 1, 4], [None, 2], [3]],
         [[None, 1, 3], [None, 4], [2]],
         [[None, 1, 4], [None, 3], [2]],
         [[None, 2, 3], [None, 4], [1]],
         [[None, 2, 4], [None, 3], [1]]]
    """
    if skp is None:
        return StandardSkewTableaux()
    elif isinstance(skp, (int, Integer)):
        return StandardSkewTableaux_size(skp)
    elif skp in SkewPartitions():
        return StandardSkewTableaux_shape(SkewPartition(skp))
    else:
        raise TypeError("Invalid argument")

# Allow unpickling old pickles using deprecated class names
from sage.structure.sage_object import register_unpickle_override
register_unpickle_override('sage.combinat.skew_tableau', 'StandardSkewTableaux_n',  StandardSkewTableaux_size)
register_unpickle_override('sage.combinat.skew_tableau', 'SemistandardSkewTableaux_n',  SemistandardSkewTableaux_size)
register_unpickle_override('sage.combinat.skew_tableau', 'SemistandardSkewTableaux_nmu',  SemistandardSkewTableaux_size_weight)
register_unpickle_override('sage.combinat.skew_tableau', 'SemistandardSkewTableaux_p',  SemistandardSkewTableaux_shape)
register_unpickle_override('sage.combinat.skew_tableau', 'SemistandardSkewTableaux_pmu',  SemistandardSkewTableaux_shape_weight)
register_unpickle_override('sage.combinat.skew_tableau', 'StandardSkewTableaux_skewpartition',  StandardSkewTableaux_shape)
#register_unpickle_override('sage.combinat.skew_tableau', 'SkewTableau_class',  SkewTableau_class)
