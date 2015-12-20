r"""
Fast word datatype using an array of unsigned char.

By convention, if stop is SIZE_T_MAX the word is considered to be infinite.
"""
#*****************************************************************************
#       Copyright (C) 2014 Vincent Delecroix <20100.delecroix@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include 'sage/ext/interrupt.pxi'
include 'sage/ext/stdsage.pxi'
include "sage/data_structures/bitset.pxi"

cimport cython

from sage.ext.memory cimport check_malloc, check_realloc

from sage.rings.integer cimport Integer, smallInteger
from sage.rings.rational cimport Rational
from sage.rings.infinity import Infinity
from libc.string cimport memcpy, memcmp, memset


from cpython.number cimport PyIndex_Check, PyNumber_Check
from cpython.sequence cimport PySequence_Check
from cpython.slice cimport PySlice_Check, PySlice_GetIndicesEx

import itertools

# the maximum value of a size_t
cdef size_t SIZE_T_MAX = -(<size_t> 1)

cdef extern from "Python.h":
    Py_ssize_t PY_SSIZE_T_MAX

def reversed_word_iterator(WordDatatype_char_finite w):
    r"""
    This function exists only because it is not possible to use yield in the
    special method ``__reversed__``.

    EXAMPLES::

        sage: W = Words([0,1,2])
        sage: w = W([0,1,0,0,1,2])
        sage: list(reversed(w)) # indirect doctest
        [2, 1, 0, 0, 1, 0]
    """
    cdef ssize_t i
    cdef unsigned char * p = w._data[0]
    for i in range((<Py_ssize_t>w._stop)-1, (<Py_ssize_t>w._start)-1, -1):
        yield p[i]

cdef class WordDatatype_char(WordDatatype):
    r"""
    A base class for :class:`WordDatatype_char_finite` and
    :class:`WordDatatype_char_infinite`.

    Currently, only handles letters in [0,255].
    """
    def __dealloc__(self):
        r"""
        Deallocate memory only if self uses it own memory.

        Note that ``sage_free`` will not deallocate memory if self is the
        master of another word.
        """
        # it is strictly forbidden to access _master here! (because Cython
        # deallocate Python attributes before calling this method, _master is
        # likely to be None in any case)
        if self._is_slice == 0:
            sage_free(self._data[0])  # free the space used by the word
            sage_free(self._data)     # free the pointer

    def is_empty(self):
        r"""
        Return whether the word is empty.

        EXAMPLES::

            sage: W = Words([0,1,2])
            sage: W([0,1,2,2]).is_empty()
            False
            sage: W([]).is_empty()
            True
        """
        return not self

    cdef _update_word_up_to(self, size_t n):
        pass

cdef class WordDatatype_char_finite(WordDatatype_char):
    def __init__(self, parent, data):
        r"""
        Constructor

        TESTS::

            sage: W = Words([0,1,2,3])
            sage: W([0,1,2,3])
            word: 0123
            sage: W(iter([0,1,2,3]))
            word: 0123
        """
#        print "WordDatatype_char_finite.__init__({}, {}, {}".format(
#                id(self), parent, data)
        self._parent = parent

        if not PySequence_Check(data):
            data = list(data)
        self._set_data(data)

    @cython.boundscheck(False) # assume that indexing will not cause any IndexErrors
    @cython.wraparound(False)  # not check not correctly handle negative indices
    cdef _set_data(self, data):
        r"""
        set the attribute ._data and ._stop from the sequence data
        (usually data is a word, a tuple or a list)
        """
        cdef size_t i
        self._stop = len(data)

        self._stop = len(data)
        self._data = <unsigned char **> check_malloc(sizeof(unsigned char *))
        self._data[0] = <unsigned char *> check_malloc((self._stop-self._start) * sizeof(unsigned char))
        for i in range(self._stop):
            self._data[0][i] = data[i]

    def __nonzero__(self):
        r"""
        Test whether the word is not empty.

        EXAMPLES::

            sage: W = Words([0,3,5])
            sage: bool(W([0,3,3,5]))
            True
            sage: bool(W([]))
            False
        """
        return self._start != self._stop

    def __len__(self):
        r"""
        Return the length of the word as a Python integer.

        TESTS::

            sage: W = Words([0,1,2,3])
            sage: w = W([0,1,2,0,3,2,1])
            sage: len(w)
            7
            sage: type(len(w))
            <type 'int'>
        """
        return self._stop - self._start

    def length(self):
        r"""
        Return the length of the word as a Sage integer.

        EXAMPLES::

            sage: W = Words([0,1,2,3,4])
            sage: w = W([0,1,2,0,3,2,1])
            sage: w.length()
            7
            sage: type(w.length())
            <type 'sage.rings.integer.Integer'>
            sage: type(len(w))
            <type 'int'>
        """
        return smallInteger(self._stop - self._start)

    def __hash__(self):
        r"""
        Return the hash value.

        EXAMPLES::

            sage: W = Words([0,1,2])
            sage: hash(W([0,1,0,1,0,0,0], datatype='list'))
            102060647
            sage: hash(W([0,1,0,1,0,0,0], datatype='char'))
            102060647
        """
        cdef int res = 5381
        cdef size_t i
        cdef unsigned char * data = self._data[0]
        if self._hash is None:
            for i in range(self._start, min(self._start+1024, self._stop)):
                res = ((res << 5) + res) + data[i]
            self._hash = res
        return self._hash

    def letters(self):
        r"""
        Return the list of letters that appear in this word, listed in the
        order of first appearance.

        EXAMPLES::

            sage: W = Words(5)
            sage: W([1,3,1,2,2,3,1]).letters()
            [1, 3, 2]
            sage: W().letters()
            []
        """
        cdef bitset_t seen
        bitset_init(seen, 256) # allocation + initialization to 0

        cdef size_t i
        cdef list res = []
        cdef unsigned char letter
        cdef unsigned char * data = self._data[0]
        for i in range(self._start, self._stop):
            letter = data[i]
            if not bitset_in(seen, letter):
                bitset_add(seen, letter)
                res.append(letter)
        bitset_free(seen)
        return res

    cdef WordDatatype_char_finite _new_c(self, unsigned char ** data, start, stop, WordDatatype_char_finite master):
        r"""
        TO DISCUSS: in Integer (sage.rings.integer) this method is actually an
        external function. But we might want to have several possible inheritance.
        """
        cdef type t = type(self)
        cdef WordDatatype_char_finite other = WordDatatype_char_finite.__new__(t)

        if data == NULL:
            other._data = <unsigned char **> check_malloc(sizeof(unsigned char *))
            other._data[0] = NULL
        else:
            other._data = data

        other._start = start
        other._stop = stop
        other._parent = self._parent

        if master is None:
            other._master = master
            other._is_slice = 0
        elif master._master is None:
            other._master = master
            other._is_slice = 1
        else:
            other._master = master._master
            other._is_slice = 1

        return other

    def __richcmp__(self, other, op):
        r"""
        INPUT:

        - ``other`` -- a word (WordDatatype_char_finite)
        - ``op`` -- int, from 0 to 5

        TESTS::

            sage: W = Words(range(100))
            sage: w = W(range(10) * 2)
            sage: w == w
            True
            sage: w != w
            False
            sage: w[:-1] != w[1:]
            True
            sage: w < w[1:] and w[1:] > w
            True
            sage: w > w[1:] or w[1:] < w
            False

            sage: W() == W([0])
            False
            sage: W([0]) == W()
            False
            sage: W() == W()
            True
        """
        # 0: <
        # 1: <=
        # 2: ==
        # 3: !=
        # 4: >
        # 5: >=
        if not isinstance(other, WordDatatype_char_finite):
            return NotImplemented

        # word of different lengths are not equal!
        cdef size_t self_len = (<WordDatatype_char_finite> self)._stop - (<WordDatatype_char_finite> self)._start
        cdef size_t other_len = (<WordDatatype_char_finite> other)._stop - (<WordDatatype_char_finite> other)._start
        if (op == 2 or op == 3) and self_len != other_len:
            return op == 3

        cdef int test = (<WordDatatype_char_finite> self)._lexico_cmp(other)
        if test < 0:
            return op < 2 or op == 3
        elif test > 0:
            return op > 3
        else:
            return op == 1 or op == 2 or op == 5

    def __cmp__(self, other):
        r"""
        INPUT:

        - ``other`` -- a word (WordDatatype_char_finite)

        TESTS::

            sage: W = Words([1,3,5])
            sage: cmp(W([3,1,3]), W([3,1,3]))
            0
            sage: cmp(W([1,5,3,3]), W([1,5,5]))
            -1
            sage: cmp(W(), W())
            0
        """
        if not isinstance(other, WordDatatype_char_finite):
            return NotImplemented

        cdef int test = self._lexico_cmp(other)
        if test:
            return test
        return (<Py_ssize_t> (self._stop - self._start)) - \
               (<Py_ssize_t> ((<WordDatatype_char_finite> other)._stop - (<WordDatatype_char_finite>other)._start))

    cdef int _lexico_cmp(self, WordDatatype_char_finite other) except -2:
        r"""
        Lexicographic comparison of self and other up to
        the letter at position min(len(self),len(other))
        """
        cdef size_t l = min(self._stop - self._start, other._stop - other._start)

        sig_on()
        cdef int test = memcmp(
                      <void *> (<WordDatatype_char_finite> self)._data[0] + self._start,
                      <void *> (<WordDatatype_char_finite> other)._data[0] + other._start,
                      l * sizeof(unsigned char))
        sig_off()

        return (test > 0) - (test < 0)

    def __getitem__(self, key):
        r"""
        INPUT:

        - ``key`` -- index

        TESTS::

            sage: W = Words([0,1,2,3])
            sage: w = W([0,1,0,2,0,3,1,2,3])
            sage: w[0]
            0
            sage: w[2]
            0
            sage: w[1:]
            word: 10203123
            sage: w[5::-2]
            word: 321

            sage: w[1:8][2:4] == w[3:5]
            True
            sage: w[1:7][5:0:-1] == w[2:7][::-1]
            True

            sage: w = W([randint(0,3) for _ in range(20)])
            sage: list(w) == [w[i] for i in range(len(w))]
            True

            sage: w['foo':'bar']
            Traceback (most recent call last):
            ...
            TypeError: slice indices must be integers or None or have an __index__ method

        Check a weird behavior of PySlice_GetIndicesEx (:trac:`17056`)::

            sage: w[1:0]
            word:
        """
        cdef Py_ssize_t i, start, stop, step, slicelength
        cdef unsigned char ** data
        cdef size_t j,k
        if PySlice_Check(key):
            # here the key is a slice
            PySlice_GetIndicesEx(key,
                    (self._stop - self._start),
                    &start, &stop, &step,
                    &slicelength)
            if slicelength == 0:
                return self._new_c(NULL, 0, 0, None)
            if step == 1:
                return self._new_c(self._data, self._start+start, self._start+stop, self)
            data = <unsigned char **> check_malloc(sizeof(unsigned char *))
            data[0] = <unsigned char *> check_malloc(slicelength * sizeof(unsigned char))
            j = 0
            for k in range(start,stop,step):
                data[0][j] = self._data[0][self._start+k]
                j += 1
            return self._new_c(data, 0, slicelength, None)

        elif PyIndex_Check(key):
            # here the key is an int
            i = key    # cast key into a size_t
            if i < 0:
                i += (self._stop - self._start)
            if i < 0 or i >= (self._stop - self._start):
                raise IndexError("word index out of range")
            return self._data[0][self._start+i]

        raise TypeError("word indices must be integers")

    def __iter__(self):
        r"""
        Iterator over the letter of self

        EXAMPLES::

            sage: W = Words([0,1,2,3])
            sage: for i in W([0,0,1,0]):  # indirect doctest
            ....:     print i,
            0 0 1 0
            sage: for i in W():
            ....:     print i
        """
        cdef size_t i
        for i in range(self._start, self._stop):
            yield self._data[0][i]

    def __reversed__(self):
        r"""
        Reversed iterator over the letter of self

        EXAMPLES::

            sage: W = Words([0,1,2,3])
            sage: w = [0,0,1,2,0]
            sage: list(reversed(w)) # indirect doctest
            [0, 2, 1, 0, 0]
            sage: list(reversed(w[1:4])) # indirect doctest
            [2, 1, 0]

        TESTS::

            sage: list(reversed(W([])))
            []
            sage: list(reversed(W([1])))
            [1]
        """
        return reversed_word_iterator(self)

    cdef _concatenate(self, WordDatatype_char_finite other):
        cdef unsigned char ** data
        cdef size_t ls = self._stop - self._start
        cdef size_t lo = other._stop - other._start

        if ls == 0:
            return other
        elif lo == 0:
            return self

        data = <unsigned char **> check_malloc(sizeof(unsigned char *))
        data[0] = <unsigned char *> check_malloc((ls + lo) * sizeof(unsigned char))

        sig_on()
        memcpy(data[0], self._data[0]+self._start, ls * sizeof(unsigned char))
        memcpy(data[0] + ls, other._data[0]+other._start, lo * sizeof(unsigned char))
        sig_off()

        return self._new_c(data, 0, ls+lo, None)

    def __mul__(self, other):
        r"""
        Concatenation of ``self`` and ``other``.

        TESTS:

            sage: W = Words(IntegerRange(1,18))
            sage: W([3,1]) * W([2,2])
            word: 3122

        The result is automatically converted to a WordDatatype_char_finite. Currently we can
        even do::

            sage: w = W([1,2,3,1])
            sage: w * [4,5,4,5]
            word: 12314545
        """
        cdef WordDatatype_char_finite w

        if isinstance(other, WordDatatype_char_finite):
            return (<WordDatatype_char_finite> self)._concatenate(other)

        elif PySequence_Check(other):
            # we convert other to a WordDatatype_char and perform the concatenation
            return (<WordDatatype_char_finite> self)._concatenate(self._parent(other))

        raise TypeError("not able to initialize a word from {}".format(other))

    def __pow__(self, exp, mod):
        r"""
        Power

        INPUT:

        -  ``exp``  - an integer, a rational, a float number or plus infinity.

        TESTS::

            sage: W = Words(range(20))
            sage: w = W([0,1,2,3])
            sage: w
            word: 0123
            sage: w ** (1/2)
            word: 01
            sage: w ** 2
            word: 01230123
            sage: w ** 3
            word: 012301230123
            sage: w ** (7/2)
            word: 01230123012301
            sage: len(((w ** 2) ** 3) ** 5) == len(w) * 2 * 3 * 5
            True

        Infinite exponents::

            sage: W([0,1]) ** Infinity
            word: 0101010101010101010101010101010101010101...
        """
        if not PyNumber_Check(exp):
            raise ValueError("the exponent must be a number or infinity")
        if mod is not None:
            raise ValueError("a word can not be taken modulo")

        if exp == float('inf'):
            fcn = lambda n: self[n % self.length()]
            return self._parent.shift()(fcn, datatype='callable')

        if exp < 0:
            raise ValueError("can not take negative power of a word")

        cdef WordDatatype_char_finite w = self
        cdef size_t l = w._stop - w._start
        cdef size_t i, rest

        if type(exp) is Rational:
            if l % exp.denominator():
                raise ValueError("undefined")
            i = exp.floor()
            rest = (exp - exp.floor()) * l
        else:
            i = exp
            rest = 0

        cdef unsigned char ** data

        # first handle the cases (i*length + rest) <= length and return the
        # corresponding prefix of self
        if i == 1 and rest == 0:
            return self
        if l == 0:
            return w._new_c(NULL, 0, 0, None)
        if i == 0:
            if rest == 0:
                return w._new_c(NULL, 0, 0, None)
            else:
                return w._new_c(w._data, w._start, w._start+rest, w)

        # now consider non trivial powers
        if l > SIZE_T_MAX / (i+1):
            raise OverflowError("the length of the result is too large")
        cdef size_t new_length = l * i + rest

        data = <unsigned char **> check_malloc(sizeof(unsigned char *))
        data[0] = <unsigned char *> check_malloc(new_length * sizeof(unsigned char))

        cdef Py_ssize_t j = l
        memcpy(data[0], w._data[0], j * sizeof(unsigned char))
        while 2*j < new_length:
            memcpy(data[0] + j, data[0], j * sizeof(unsigned char))
            j *= 2
        memcpy(data[0] + j, data[0], (new_length - j) * sizeof(unsigned char))

        return w._new_c(data, 0, new_length, None)

    def is_square(self):
        r"""
        Returns True if self is a square, and False otherwise.

        EXAMPLES::

            sage: w = Word([n % 4 for n in range(48)], alphabet=[0,1,2,3])
            sage: w.is_square()
            True

        ::

            sage: w = Word([n % 4 for n in range(49)], alphabet=[0,1,2,3])
            sage: w.is_square()
            False
            sage: (w*w).is_square()
            True

        TESTS:

        The above tests correspond to the present class (char)::

            sage: type(w)
            <class 'sage.combinat.words.word.FiniteWord_char'>

        ::

            sage: Word([], alphabet=[0,1]).is_square()
            True
            sage: Word([0], alphabet=[0,1]).is_square()
            False
            sage: Word([0,0], alphabet=[0,1]).is_square()
            True
        """
        cdef size_t length = self._stop - self._start
        if length % 2 != 0:
            return False
        else:
            length //= 2
            return memcmp(self._data[0],
                          self._data[0] + length,
                          length * sizeof(unsigned char)) == 0

    def has_prefix(self, other):
        r"""
        Test whether ``other`` is a prefix of ``self``.

        INPUT:

        - ``other`` -- a word or a sequence (e.g. tuple, list)

        EXAMPLES::

            sage: W = Words([0,1,2])
            sage: w = W([0,1,1,0,1,2,0])
            sage: w.has_prefix([0,1,1])
            True
            sage: w.has_prefix([0,1,2])
            False
            sage: w.has_prefix(w)
            True
            sage: w.has_prefix(w[:-1])
            True
            sage: w.has_prefix(w[1:])
            False

        TESTS:

        :trac:`19322`::

            sage: W = Words([0,1,2])
            sage: w = W([0,1,0,2])
            sage: w.has_prefix(words.FibonacciWord())
            False

            sage: w.has_prefix([0,1,0,2,0])
            False
            sage: w.has_prefix([0,1,0,2])
            True
            sage: w.has_prefix([0,1,0])
            True
        """
        cdef size_t i
        cdef WordDatatype_char_finite w
        cdef size_t ls = self._stop - self._start
        cdef size_t lo

        if isinstance(other, WordDatatype_char_finite):
            # C level
            w = <WordDatatype_char_finite> other
            lo = w._stop - w._start

            if lo == 0:
                return True
            elif (lo > ls):
                return False
            else:
                return (memcmp(self._data[0]+self._start,
                          w._data[0]+w._start,
                          lo * sizeof(unsigned char)) == 0)

        elif PySequence_Check(other):
            # python level
            from sage.combinat.words.infinite_word import InfiniteWord_class
            if isinstance(other, InfiniteWord_class):
                return False
            lo = len(other)
            if lo > ls:
                return False
            self._update_word_up_to(lo)
            for i in range(lo):
                if other[i] != self._data[0][self._start+i]:
                    return False
            return True

        raise TypeError("not able to initialize a word from {}".format(other))

    #TODO: move into Word_datatype_char
    def longest_common_prefix(self, other):
        r"""
        Return the longest common prefix of this word and ``other``.

        EXAMPLES::

            sage: W = Words([0,1,2])
            sage: W([0,1,0,2]).longest_common_prefix([0,1])
            word: 01
            sage: W([0,1,0,2]).longest_common_prefix([0])
            word: 0
            sage: W([0,1,0,2]).longest_common_prefix([])
            word:

            sage: u = W([0,1,0,0,1,1,1])
            sage: v = W([0,1,0,2,0,1,0,0,1])
            sage: u.longest_common_prefix(v)
            word: 010
            sage: v.longest_common_prefix(u)
            word: 010
            sage: u[:5].longest_common_prefix(v[4:9])
            word: 01001

        Using infinite words is also possible (and the return type is also a
        of the same type as ``self``)::

            sage: W([0,1,0,0]).longest_common_prefix(words.FibonacciWord())
            word: 0100
            sage: type(_)
            <class 'sage.combinat.words.word.FiniteWord_char'>

        An example of an intensive usage::

            sage: W = Words([0,1])
            sage: w = words.FibonacciWord()
            sage: w = W(list(w[:5000]))
            sage: L = [[len(w[n:].longest_common_prefix(w[n+fibonacci(i):]))
            ....:      for i in range(5,15)] for n in range(1,1000)]
            sage: for n,l in enumerate(L):
            ....:     if l.count(0) > 4: print n+1,l
            375 [0, 13, 0, 34, 0, 89, 0, 233, 0, 233]
            376 [0, 12, 0, 33, 0, 88, 0, 232, 0, 232]
            608 [8, 0, 21, 0, 55, 0, 144, 0, 377, 0]
            609 [7, 0, 20, 0, 54, 0, 143, 0, 376, 0]
            985 [0, 13, 0, 34, 0, 89, 0, 233, 0, 610]
            986 [0, 12, 0, 33, 0, 88, 0, 232, 0, 609]

        TESTS::

            sage: W = Words([0,1,2])
            sage: w = W([0,2,1,0,0,1])
            sage: w.longest_common_prefix(0)
            Traceback (most recent call last):
            ...
            TypeError: unsupported input 0
        """
        cdef WordDatatype_char_finite w
        cdef size_t i = SIZE_T_MAX
        cdef size_t ls
        cdef size_t lo

        if isinstance(other, WordDatatype_char_finite):
            # C level
            # (this can be much faster if we allow to compare larger memory
            # zones)
            w = <WordDatatype_char_finite> other
            ls = self._stop - self._start
            lo = w._stop - w._start
            for i in range(min(ls,lo)):
                if self._data[0][self._start+i] != w._data[0][w._start+i]:
                    break
            else:
                if ls <= lo:
                    return self
                else:
                    return other

            return self._new_c(self._data, self._start, self._start+i, self)

        elif PySequence_Check(other):
            # Python level
            # we avoid to call len(other) since it might be an infinite word
            for i,a in enumerate(itertools.islice(other, self._stop-self._start)):
                if self._data[0][self._start+i] != a:
                    break
            else:
                if i == SIZE_T_MAX:
                    # empty word
                    return self._new_c(NULL, 0, 0, None)
                else:
                    i += 1

            return self._new_c(self._data, self._start, self._start+i, self)

        raise TypeError("unsupported input {}".format(other))

    def longest_common_suffix(self, other):
        r"""
        Return the longest common suffix between this word and ``other``.

        EXAMPLES::

            sage: W = Words([0,1,2])
            sage: W([0,1,0,2]).longest_common_suffix([2,0,2])
            word: 02
            sage: u = W([0,1,0,0,1])
            sage: v = W([1,0,0,1,2,0,0,1])
            sage: u.longest_common_suffix(v)
            word: 001
            sage: v.longest_common_suffix(u)
            word: 001
            sage: u[1:].longest_common_suffix(v[1:4])
            word: 001

        TESTS::

            sage: W = Words([0,1,2])
            sage: w = W([0,2,1,0,0,1])
            sage: w.longest_common_suffix(0)
            Traceback (most recent call last):
            ...
            TypeError: unsupported input 0
        """
        cdef WordDatatype_char_finite w
        cdef size_t i = 0
        cdef size_t ls = self._stop - self._start
        cdef size_t lo

        if isinstance(other, WordDatatype_char_finite):
            # C level
            # (this can be much faster if we could compare larger memory
            # zones)
            w = other
            lo = w._stop - w._start
            for i in range(min(ls,lo)):
                if self._data[0][self._stop-i-1] != w._data[0][w._stop-i-1]:
                    break
            else:
                if ls <= lo:
                    return self
                else:
                    return other

            return self._new_c(self._data, self._stop-i, self._stop, self)

        elif PySequence_Check(other):
            # Python level
            lo = len(other)
            for i in range(min(ls,lo)):
                if self._data[0][self._stop-i-1] != other[lo-i-1]:
                    break
            else:
                if ls <= lo:
                    return self
                elif i == 0:
                    # empty word
                    return self._new_c(NULL, 0, 0, None)
                else:
                    i += 1

            return self._new_c(self._data, self._stop-i, self._stop, self)

        raise TypeError("unsupported input {}".format(other))

cdef class Slicer(object):
    r"""
    A simple class to handle non standard slices of infinite words.

    This class is mainly a helper for ``__getitem__`` of infinite words. We made
    it a class to allows pickling support. Let ``w`` be an infinite words. When
    considering a slice, if step is ``1`` (like ``w[4::1]``) then the infinite
    words will share their memory. But if it is not (like ``w[5::2]``) then we
    need to create a new word that will read inside the first one. This class is
    built to handle this second case.
    """
    cdef WordDatatype_char_infinite word
    cdef Py_ssize_t start
    cdef Py_ssize_t step

    def __init__(self, word, start, step):
        r"""
        INPUT:

        - ``word`` - an infinite word

        - ``start`` - the starting point

        - ``step`` - step

        TESTS::

            sage: from sage.combinat.words.word_char import Slicer
            sage: W = InfiniteWords([0,1,2,3,4,5,6,7])
            sage: w = W(lambda n: n%5)
            sage: s = Slicer(w, 3, 4)
            sage: all(s(i) == w[3 + 4*i] for i in range(10))
            True
        """
        self.word = word
        self.start = start
        self.step = step

    def __reduce__(self):
        r"""
        Pickling support

        TESTS::

            sage: from sage.combinat.words.word_char import Slicer
            sage: W = InfiniteWords([0,1,2,3,4,5,6,7])
            sage: w = W(lambda n: n%5)
            sage: s = Slicer(w, 3, 4)
            sage: t = loads(dumps(s))
            sage: all(t(i) == s(i) for i in range(10))
            True
        """
        return Slicer, (self.word, self.start, self.step)

    def __call__(self, n):
        r"""
        TESTS::

            sage: from sage.combinat.words.word_char import Slicer
            sage: W = InfiniteWords([0,1,2,3,4,5,6,7])
            sage: w = W(lambda n: n%5)
            sage: s = Slicer(w, 3, 4)
            sage: all(s(i) == w[3 + 4*i] for i in range(10))
            True
        """
        self.word._update_word_up_to(self.start + self.step*n + 1)
        return self.word._data[0][self.start + self.step*n]

def slice_unpickle(data, start, stop, step):
    r"""
    Return ``data[start:stop:step]``.

    This method is used for unpickling slice of infinite words (see the method
    ``__reduce__`` of :class:`WordDatatype_char_infinite`.

    TESTS::

        sage: from sage.combinat.words.word_char import slice_unpickle
        sage: slice_unpickle(range(10), 2, 7, 3)
        [2, 5]
    """
    return data[start:stop:step]

cdef class WordDatatype_char_infinite(WordDatatype_char):
    r"""
    In order to build an infinite word with underlying C datatype being an array
    of C ``unsigned char`` one needs to inherit from this class and either
    implement ``_new_slice`` (for Python classes) or ``_extend_word`` (for
    Cython extension classes).

    EXAMPLES:

    To create your own class of words in Python, you need:

    - to inherit from this class and `InfiniteWord_class`

    - to call the method ``_init_c_data(alloc_size)`` inside the constructor
    ``__init__``. The argument ``alloc_size`` is the amount of memory
    preallocated for the word.

    - to implement a method ``_new_slice()`` that return a new chunk of the word
    each time

    As an example, we create the Kolakoski word::

        sage: from sage.combinat.words.word_char import WordDatatype_char_infinite
        sage: from sage.combinat.words.infinite_word import InfiniteWord_class
        sage: class MyKolakoskiWord(WordDatatype_char_infinite, InfiniteWord_class):
        ....:     def __init__(self):
        ....:         self._init_c_data(13)   # initialize the C data structure
        ....:         self._parent = InfiniteWords([1,2])
        ....:         self._letter = None     # current letter
        ....:         self._pos = None        # current position
        ....:     def _new_slice(self):
        ....:         if self._pos is None:
        ....:             self._pos = 2
        ....:             self._letter = 1
        ....:             return [1,2,2]
        ....:         else:
        ....:             n = self[self._pos]
        ....:             self._pos += 1
        ....:             letter = self._letter
        ....:             self._letter = 1 if letter == 2 else 2
        ....:             return [letter]*n

    Then::

        sage: w = MyKolakoskiWord()
        sage: w
        word: 1221121221221121122121121221121121221221...
        sage: w.parent()
        Infinite words over {1, 2}

    For examples in Cython, you might have a look at the code of
    :class:`WordDatatype_callable_char`, :class:`WordDatatype_iter_char` or
    :class:`WordDatatype_substitutive_char`.
    """
    cpdef int _init_c_data(self, size_t alloc_size) except -1:
        r"""
        Any words should call this method at initialization.

        This method initialize the following attributes:

        - ``self._data``
        - ``self._alloc_size``

        It must be call in the constructor (i.e. ``__init__``) *before* any
        access to the letters.

        INPUT:

        - ``alloc_size`` - the size of the memory to reserve for the word
          (beyond that limit there will be dynamic reallocation)

        TESTS::

            sage: W = InfiniteWords([0,1])
            sage: w = W(lambda n: 0)          # indirect doctest
        """
        self._data = <unsigned char **> check_malloc(sizeof(unsigned char *))
        self._alloc_size = alloc_size
        self._data[0] = <unsigned char *> check_malloc(alloc_size * sizeof(unsigned char))
        return 0

    cdef int _realloc(self, size_t new_size) except -1:
        r"""
        Allocate more memory for the word if needed.

        INPUT:

        - ``new_size`` - the minimum size of the new allocated memory
        """
        if new_size >= self._alloc_size:
            self._alloc_size = max(new_size, <size_t> ((1.3) * self._alloc_size
))
            self._data[0] = <unsigned char *> check_realloc(self._data[0], self._alloc_size * sizeof(unsigned char))

        return 0

    def __reduce__(self):
        r"""
        TESTS::

            sage: w = WordMorphism({1:[4], 4:[4,1]}).fixed_point(4)
            sage: u = w[10:][4:]
            sage: u
            word: 1441414414414144144141441414414414144141...
            sage: loads(dumps(u))
            word: 1441414414414144144141441414414414144141...
        """
        if self._is_slice:
            return slice_unpickle, (self._master, self._start, None, 1)
        else:
            raise NotImplementedError

    def __hash__(self):
        r"""
        A hash function that takes care of the first 1024 letters

        EXAMPLES::

            sage: w = WordMorphism({0:[0,1],1:[1]}).fixed_point(0)
            sage: hash(w)
            -1533913820
        """
        cdef int res = 5381
        cdef int i
        cdef unsigned char * data
        if (<WordDatatype> self)._hash is None:
            self._update_word_up_to(1025)
            data = self._data[0]
            for i in range(1024):
                res = ((res << 5) + res) + data[i]
            (<WordDatatype> self)._hash = res
        return (<WordDatatype> self)._hash

    def __iter__(self):
        r"""
        TESTS::

            sage: W = InfiniteWords([0,1,2])
            sage: w = W(lambda n: Integer(n).popcount() % 3)
            sage: it = iter(w)      # indirect doctest
            sage: for i in range(20):
            ....:     assert it.next() == w[i]
        """
        cdef size_t i = self._start
        while True:
            while i < self._current_stop:
                yield self._data[0][i]
                i += 1
            self._extend_word()

    cpdef WordDatatype_char_infinite _infinite_slice(self, size_t start, size_t step):
        r"""
        Return an infinite slice

        TESTS::

            sage: W = InfiniteWords([0,1])
            sage: w = W(lambda n: Integer(n).popcount()%2)
            sage: w[1:]    # indirect doctest
            word: 1101001100101101001011001101001100101100...
            sage: type(_)
            <class 'sage.combinat.words.word.InfiniteWord_char'>

            sage: w[::3]   # indirect doctest
            word: 0000000100000010000000010001110100000001...
            sage: type(_)
            <class 'sage.combinat.words.word.InfiniteWord_callable_char'>
        """
        cdef WordDatatype_char_infinite w
        P = self._parent

        start += self._start

        if step == 1:
            if start == 0:
                return self
            else:
                w = WordDatatype_char_infinite.__new__(P._element_classes['char'])
                w._parent = P
                w._start = start
                w._data = self._data
                w._is_slice = 1
                w._master = self._master if self._is_slice else self
                return w
        else:
            if self._is_slice:
                return P(Slicer(self._master, start, step))
            else:
                return P(Slicer(self, start, step))

    cpdef WordDatatype_char_finite _finite_slice(self, Py_ssize_t start, Py_ssize_t stop, Py_ssize_t step, Py_ssize_t slicelength):
        r"""
        Return a finite slice

        There is no check on input. This method is mostly a helper for
        ``__getitem__``.

        INPUT:

        - ``start`` - where does the slice starts (relative to this word)

        - ``stop`` - where does the slice stops (relative to this word)

        - ``step`` - the step

        - ``slicelength`` -- the length of the slize

        TESTS::

            sage: W = InfiniteWords([0,1,2])
            sage: w = W(lambda n: Integer(n).popcount() % 3)
            sage: l = [Integer(n).popcount()%3 for n in range(20)]
            sage: for (start,stop,step) in [(None,7,None), (2,9,None),
            ....:         (None,14,2), (3,17,2), (10,None,-1),
            ....:         (14,10,-1), (1,0,None), (0,1,-1),
            ....:         (1,5,2),  (2,5,2), (1,5,3), (2,5,3),
            ....:         (5,1,-2), (5,2,-2), (5,1,-3), (5,2,-3)]:
            ....:    assert list(w[start:stop:step]) == l[start:stop:step]
        """
        cdef WordDatatype_char_finite w
        cdef unsigned char ** data
        cdef size_t j
        cdef Py_ssize_t k

        if step > 0:
            self._update_word_up_to(stop)
        else:
            self._update_word_up_to(start)

        P = self._parent.factors()
        w = WordDatatype_char_finite.__new__(P._element_classes['char'])
        w._parent = P

        start += self._start
        stop += self._start

        if slicelength == 0:
            w._is_slice = 0
            w._data = <unsigned char **> check_malloc(sizeof(unsigned char *))
            w._data[0] = NULL
        elif step == 1:
            w._start = start
            w._stop = stop
            w._is_slice = 1
            w._master = self
            w._data = self._data
        else:
            data = <unsigned char **> check_malloc(sizeof(unsigned char *))
            data[0] = <unsigned char *> check_malloc(slicelength * sizeof(unsigned char))
            j = 0
            for k in range(start,stop,step):
                data[0][j] = self._data[0][self._start+k]
                j += 1
            w._is_slice = 0
            w._master = None
            w._data = data
            w._start = 0
            w._stop = slicelength

        return w

    def __getitem__(self, key):
        r"""
        INPUT:

        - ``key`` -- index

        TESTS::

            sage: W = Words([0,1,2,3])
            sage: w = W([0,1,0,2,0,3,1,2,3])
            sage: w[0]
            0
            sage: w[2]
            0
            sage: w[1:]
            word: 10203123
            sage: w[5::-2]
            word: 321

            sage: w[1:8][2:4] == w[3:5]
            True
            sage: w[1:7][5:0:-1] == w[2:7][::-1]
            True

            sage: w = W([randint(0,3) for _ in range(20)])
            sage: list(w) == [w[i] for i in range(len(w))]
            True

            sage: w['foo':'bar']
            Traceback (most recent call last):
            ...
            TypeError: slice indices must be integers or None or have an __index__ method

        Check a weird behavior of PySlice_GetIndicesEx (:trac:`17056`)::

            sage: w[1:0]
            word:


        TESTS::

            sage: W = InfiniteWords([0,1])
            sage: w = W(lambda n: Integer(n).popcount()%2)
            sage: w1 = w[1:]    # indirect doctest
            sage: w1[250] == w[251] and w1[251] == w[252] and w1[252] == w[253]
            True
            sage: w1[2046:2051] == w[2047:2052]
            True

            sage: w4 = w1[3:]
            sage: w4[250] == w[254] and w4[251] == w[255] and w4[252] == w[256]
            True
            sage: w4[1037:1055] == w[1041:1059]
            True

            sage: u2 = w[1::2]
            sage: all(u2[i] == w[1+2*i] for i in range(10,30))
            True
            sage: u2[2:18:3] == w[5:37:6]
            True

            sage: w[-3]
            Traceback (most recent call last):
            ...
            IndexError: infinite word index must be positive
            sage: w[1:-1]
            Traceback (most recent call last):
            ...
            IndexError: infinite word index must be positive
        """
        cdef Py_ssize_t i, start, stop, step, slicelength
        cdef slice s

        if PySlice_Check(key):
            s = key
            if s.step is None:
                step = 1
            else:
                step = s.step

            if not step:
                raise ValueError("slice step can not be 0")

            if s.stop is None and step > 0:
                # build a new infinite word
                if s.start is None:
                    start = 0
                else:
                    start = s.start

                return self._infinite_slice(start, step)

            else:
                # build a new finite word
                # we avoid to use PySlice_GetIndicesEx because of negative
                # indices are invalid here
                if s.start is None:
                    if step < 0:
                        raise IndexError("if step is negative, start must be specified")
                    else:
                        start = 0
                else:
                    start = s.start
                    if start < 0:
                        raise IndexError("infinite word index must be positive")

                if s.stop is None:
                    stop = -1
                else:
                    stop = s.stop
                    if stop < 0:
                        raise IndexError("infinite word index must be positive")

                if step > 0:
                    if stop <= start:
                        slicelength = 0
                    else:
                        slicelength = (stop - start + step - 1) // step
                else:
                    if start <= stop:
                        slicelength = 0
                    else:
                        slicelength = (stop - start + step + 1) // step
                return self._finite_slice(start, stop, step, slicelength)

        elif PyIndex_Check(key):
            # here the key is an int
            i = key
            if i < 0:
                raise IndexError("infinite word index must be positive")
            self._update_word_up_to(self._start+i+1)
            return self._data[0][self._start+i]

        raise TypeError("word indices must be integers")

    def __nonzero__(self):
        r"""
        An infinite word is always non-empty.

        EXAMPLES::

            sage: F = InfiniteWords([0,1,2])
            sage: w = F(lambda n: Integer(n).popcount() % 3)
            sage: bool(w)
            True
        """
        return True

    def is_empty(self):
        r"""
        An infinite word is always non-empty.

        EXAMPLES::

            sage: W = InfiniteWords([0,1,2])
            sage: W(lambda n: n%2).is_empty()
            False
        """
        return False

    cdef int _extend_word(self) except -1:
        r"""
        Extend the word once.

        It is up to subclasses to determine how to handle this method. Cython
        class might just implement it and Python class should implement
        `_new_slice`. See examples in :class:`WordDatatype_callable_char`,
        :class:`WordDatatype_iter_char` and
        :class:`WordDatatype_substitutive_char`.
        """
        if self._is_slice:
            if (<WordDatatype_char_infinite> self._master)._current_stop == \
                    self._current_stop:
                (<WordDatatype_char_infinite> self._master)._extend_word()
            self._current_stop = (<WordDatatype_char_infinite> self._master)._current_stop
        else:
            for letter in self._new_slice():
                if self._alloc_size == self._current_stop:
                    self._realloc(self._alloc_size)
                self._data[0][self._current_stop] = letter
                self._current_stop += 1
        return 0

    cdef _update_word_up_to(self, size_t n):
        r"""
        Make sure that the prefix of length ``n`` is computed.

        This method just call repetitively ``_extend_word``.
        """
        if self._is_slice:
            self._master._update_word_up_to(self._start + n)
            self._current_stop = \
                (<WordDatatype_char_infinite> self._master)._current_stop
        else:
            while self._current_stop < n:
                self._extend_word()

cdef class WordDatatype_callable_char(WordDatatype_char_infinite):
    r"""
    Base class for words built from a function.

    EXAMPLES::

        sage: W = InfiniteWords([0,1])
        sage: w = W(lambda n: Integer(n).popcount()%2)
        sage: w
        word: 0110100110010110100101100110100110010110...
        sage: type(w)
        <class 'sage.combinat.words.word.InfiniteWord_callable_char'>
        sage: w.length()
        +Infinity
        sage: w[113]
        0
        sage: w[10:25]
        word: 010110100101100
    """
    def __init__(self, parent, f, length=None, slice_length=23):
        r"""
        INPUT:

        - ``parent`` - the parent of the infinite word

        - ``f`` - a function that to the integer ``n`` associates the letter at
          position ``n`` in the word

        - ``length`` - ignored (compatibility with other class of infinite
          words)

        - ``slice_length`` - determine the number of letters that is added to
          the word at each update. Having a large number lowers the number of
          calls. But if ``f`` is complicated enough that might be better to set
          it to ``1``. Default is ``23``.

        TESTS::

            sage: W = InfiniteWords([2,4])
            sage: w = W(lambda n: 2 + 2*(Integer(n).popcount()%2))
            sage: w
            word: 2442422442242442422424422442422442242442...
            sage: w[:100] == w[:200:2]
            True
            sage: ldw = loads(dumps(w))
            sage: ldw.parent() == w.parent()
            True
            sage: ldw[:200] == w[:200]
            True
        """
        self._parent = parent
        self._init_c_data(13)
        self._f = f
        assert slice_length > 0
        self._slice_length = slice_length

    cdef int _extend_word(self) except -1:
        cdef int j
        self._realloc(self._current_stop + self._slice_length)
        for j in range(self._slice_length):
            self._data[0][self._current_stop] = self._f(self._current_stop)
            self._current_stop += 1
        return 0

    def __reduce__(self):
        r"""
        Pickling support.

        TESTS::

            sage: W = InfiniteWords([0,1])
            sage: w = W(is_prime)
            sage: w
            word: 0011010100010100010100010000010100000100...
            sage: loads(dumps(w))
            word: 0011010100010100010100010000010100000100...
        """
        from sage.misc.fpickle import pickle_function
        try:
            s = pickle_function(self._f)
        except Exception:
            return self._parent, (self._f, 'callable', True)
        else:
            return self._parent, (s, 'pickled_function', True)

cdef class WordDatatype_iter_char(WordDatatype_char_infinite):
    r"""
    Base class for infinite words built from an iterator.

    EXAMPLES::

        sage: from itertools import count
        sage: W = InfiniteWords([0,1])
        sage: w = W(n%2 for n in count())
        sage: w
        word: 0101010101010101010101010101010101010101...
        sage: w.length()
        +Infinity
        sage: type(w)
        <class 'sage.combinat.words.word.InfiniteWord_iter_char'>
        sage: w[10]
        0
        sage: w[11]
        1
        sage: w[10:2:-1]
        word: 01010101
    """
    def __init__(self, parent, iterator, length=None, slice_length=23):
        r"""
        INPUT:

        - ``parent`` - the parent of this infinite word

        - ``iterator`` - an iterator

        - ``length`` - ignored (compatibility with other words)

        - ``slice_length`` - the number of letters that is added to the cache at
          each time. If larger then the number of calls might be significatively
          lower. Though, if the iterator is relatively complicated, it might be
          better to set it to ``1``.

        TESTS::

            sage: W = InfiniteWords([0,1])
            sage: from itertools import count
            sage: w = W(Integer(i).is_prime() for i in count())
            sage: w[::2]
            word: 0100000000000000000000000000000000000000...
            sage: w[10:3:-1]
            word: 0001010
        """
        self._parent = parent
        self._init_c_data(13)
        self._iterator = iterator
        assert slice_length > 0
        self._slice_length = slice_length

    cdef int _extend_word(self) except -1:
        cdef int j
        self._realloc(self._current_stop + self._slice_length)
        for j in range(self._slice_length):
            self._data[0][self._current_stop] = self._iterator.next()
            self._current_stop += 1
        return 0

cdef class WordDatatype_substitutive_char(WordDatatype_char_infinite):
    r"""
    Base class for fixed point of substitution.

    EXAMPLES::

        sage: s = WordMorphism({0: [0,1], 1: [0]})
        sage: w = s.fixed_point(0)
        sage: type(w)
        <class 'sage.combinat.words.word.SubstitutiveWord_char'>
        sage: w
        word: 0100101001001010010100100101001001010010...
        sage: w[:20]
        word: 01001010010010100101

        sage: WordMorphism({1: [4], 4: [4,7], 7: [1]}).fixed_point(4)
        word: 4714474714714471447471447471471447471471...
        sage: type(_)
        <class 'sage.combinat.words.word.SubstitutiveWord_char'>
    """
    def __init__(self, parent, s, letter):
        r"""
        INPUT:

        - ``parent`` - the parent for this infinite word

        - ``s`` - the substitution

        - ``letter`` - the first letter of the word

        TESTS::

            sage: s = WordMorphism({0: [0,1], 1: [0]})
            sage: w = s.fixed_point(0)
            sage: w[::2]
            word: 0011001100010001000110011000100010001100...
        """
        assert s.is_prolongable(letter) and s.is_growing(letter)
        assert parent == s.domain().shift() == s.codomain().shift()

        self._parent = parent
        alphabet = parent.alphabet()
        d = max(alphabet)+1
        py_images = [(a,s.image(a)) for a in alphabet]
        cdef size_t M = sum(len(w) for _,w in py_images)

        self._init_c_data(2*M)

        self._images = <unsigned char **> check_malloc(d * sizeof(unsigned char *))
        self._lengths = <int *> check_malloc(d * sizeof(int))
        memset(self._lengths, 0, d * sizeof(int))
        cdef unsigned char * ptr =  <unsigned char *> check_malloc(M * sizeof(unsigned char))
        self._images[0] = ptr
        for a,w in py_images:
            self._lengths[a] = len(w)
            self._images[a] = ptr
            for j in w:
                ptr[0] = j
                ptr += 1

        self._i = 0
        memcpy((<WordDatatype_char>self)._data[0], self._images[letter], self._lengths[letter] * sizeof(unsigned char))
        self._current_stop = self._lengths[letter]

    def __dealloc__(self):
        # NOTE: the __dealloc__ of WordDatatype_char_infinite is automatically
        # called. There is no need to care about self._data here.
        sage_free(self._images[0])
        sage_free(self._images)
        sage_free(self._lengths)

    cdef int _extend_word(self) except -1:
        # add a new image of a letter
        self._i += 1
        cdef unsigned char letter = self._data[0][self._i]
        self._realloc(self._current_stop + self._lengths[letter])

        memcpy((<WordDatatype_char>self)._data[0] + self._current_stop,
                self._images[letter],
                self._lengths[letter] * sizeof(unsigned char))
        self._current_stop += self._lengths[letter]
        return 0

    def morphism(self):
        r"""
        Return the morphism from which this infinite word was built as a fixed point.

        EXAMPLES::

            sage: w = WordMorphism({1:[4], 4:[4,1]}).fixed_point(4)
            sage: w.morphism()
            WordMorphism: 1->4, 4->41
        """
        from sage.combinat.words.morphism import WordMorphism
        cdef int i,j
        W = self._parent.factors()
        d = {i: W([self._images[i][j] for j in range(self._lengths[i])]) for i in self._parent.alphabet()}
        return WordMorphism(d, domain=W)

    def __reduce__(self):
        r"""
        Pickling support.

        EXAMPLES::

            sage: w = WordMorphism({1:[4], 4:[4,1]}).fixed_point(4)
            sage: w
            word: 4144141441441414414144144141441441414414...
            sage: loads(dumps(w))
            word: 4144141441441414414144144141441441414414...
            sage: parent(_)
            Infinite words over {1, 4}
        """
        s = self.morphism()
        return s.fixed_point, (self[0],)
