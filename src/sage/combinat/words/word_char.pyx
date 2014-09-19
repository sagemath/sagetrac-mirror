r"""
TESTS::

    sage: from sage.combinat.words.word import FiniteWord_char
    sage: W = Words(IntegerRange(0,256))
    sage: w = FiniteWord_char(W, range(1,10)*2)
    sage: w
    word: 123456789123456789

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

    sage: list(w) == [w[i] for i in range(len(w))]
    True

    sage: type(len(w))
    <type 'int'>
    sage: type(w.length())
    <type 'sage.rings.integer.Integer'>

    sage: w.has_prefix([1,2,3,4])
    True
    sage: w.has_prefix([1,2,4,4])
    False
    sage: w.has_prefix(FiniteWord_char(W,[1,2,3,4]))
    True
    sage: w.has_prefix(FiniteWord_char(W,[1,2,4,4]))
    False

    sage: w.is_palindrome()
    False
    sage: (w*w[::-1]).is_palindrome()
    True
    sage: (w[:-1:]*w[::-1]).is_palindrome()
    True

    sage: w.is_lyndon()
    False
    sage: FiniteWord_char(W, range(10)+[10,10]).is_lyndon()
    True

    sage: w.is_square()
    True
    sage: w[:-2].is_square()
    False
    sage: w.is_square_free()
    False
    sage: w[:-1].is_square_free()
    True
    sage: u = FiniteWord_char(W, [randint(0,255) for i in range(10)])
    sage: (u*u).is_square()
    True
    sage: (u*u*u).is_cube()
    True

    sage: len(w.factor_set())
    127
    sage: w.rauzy_graph(5)
    Looped digraph on 9 vertices


    sage: u = FiniteWord_char(W,[1,2,3])
    sage: u.first_pos_in(w)
    0
    sage: u.first_pos_in(w[1:])
    8
"""

include 'sage/ext/interrupt.pxi'
include 'sage/ext/stdsage.pxi'

cimport cython
from sage.rings.integer cimport Integer, smallInteger
from libc.string cimport memcpy, memcmp

cdef extern from "Python.h":
    ctypedef struct PyObject:
        pass
    ctypedef struct PySliceObject:
        pass

    # check functions
    int PyIndex_Check(PyObject *o)
    int PySlice_Check(PyObject *o)
    int PySequence_Check(PyObject *o)


    int PySlice_GetIndicesEx(PySliceObject *slice, Py_ssize_t length,
            Py_ssize_t *start, Py_ssize_t *stop, Py_ssize_t *step,
            Py_ssize_t *slicelength)
    Py_ssize_t PyNumber_AsSsize_t(PyObject *o, PyObject *exc)

    # Error handling
    PyObject * PyExc_IndexError
    PyObject* PyErr_Occurred()


cdef class Word_char(object):
    cdef unsigned char * _data
    cdef size_t _length
    cdef Word_char _master
    cdef public _parent
    cdef int _hash

    def __init__(self, parent, data):
        self._parent = parent
        self._master = None
        self._data = NULL
        self._length = 0

        if data:
            if not PySequence_Check(<PyObject *>data):
                raise TypeError("not able to initialize a word from {}".format(data))

            self._set_data(data)

    @cython.boundscheck(False) # assume that indexing will not cause any IndexErrors
    @cython.wraparound(False)  # not check not correctly handle negative indices
    cdef _set_data(self, data):
        r"""
        set the attribute ._data and ._length from the sequence data
        (usually data is a word, a tuple or a list)
        """
        cdef size_t i
        self._length = len(data)
        self._data = <unsigned char *> sage_malloc(self._length * sizeof(unsigned char))
        if self._data == NULL:
            raise MemoryError

        for i in range(self._length):
            self._data[i] = data[i]

    def __dealloc__(self):
        if self._master is None:
            sage_free(self._data)

    # TO DISCUSS: in Integer (sage.rings.integer) this method is actually an
    # external function
    cdef _new_c(self, unsigned char * data, size_t length, Word_char master):
        cdef Word_char other = PY_NEW_SAME_TYPE(self)
        if HAS_DICTIONARY(self):
            other.__class__ = self.__class__
        other._data = data
        other._master = master # can be None
        other._length = length
        other._parent = self._parent

        # TODO: hack to copy the Python attributes in inherited classes
        return other

    def __hash__(self):
        r"""
        Returns the hash for this word.

        TESTS::

             sage: h = hash(Word('abc'))    # indirect test
             sage: Word('abc').__hash__() == Word('abc').__hash__()
             True
        """
        cdef int res = 5381
        cdef size_t i
        # if the hash is zero we compute it all the time!!
        if self._hash == 0:
            for i in range((<Word_char> self)._length):
                res = ((res << 5) + res) + (<Word_char> self)._data[i]
            self._hash = res
        return self._hash

    def __richcmp__(self, other, op):
        # 0: <
        # 1: <=
        # 2: ==
        # 3: !=
        # 4: >
        # 5: >=

        if not PY_TYPE_CHECK(other, Word_char):
            return NotImplemented

        # word of different lengths are not equal!
        if (op == 2 or op == 3) and (<Word_char> self)._length != (<Word_char> other)._length:
            return op == 3

        cdef int test = (<Word_char> self)._lexico_cmp(<Word_char> other)
        if test < 0:
            return op < 2 or op == 3
        elif test > 0:
            return op > 3
        else:
            return op == 1 or op == 2 or op == 5

    cdef int _lexico_cmp(self, Word_char other) except -2:
        r"""
        Lexicographic comparison of self and other up to
        the letter at position min(len(self),len(other))
        """
        cdef size_t l = min(self._length, other._length)

        sig_on()
        cdef int test = memcmp(<void *> (<Word_char> self)._data,
                      <void *> (<Word_char> other)._data,
                      l * sizeof(unsigned char))
        sig_off()

        if test == 0:
            return 0
        if test < 0:
            return -1
        else:
            return 1

    def _get_info(self):
        r"""
        Temporary function that print useful debug information
        """
        if self._master:
            print "master at %u"%(id(self._master))
        print "data at %u"%(<size_t>self._data)

    def __getitem__(self, key):
        r"""
        note: the implementation is mostly borrowed from the implementation of
        __getitem__ of lists. The main problem is that it seems hard to
        propagate the error we might get from PySlice_GetIndicesEx or
        PyNumber_asSize_t.
        """
        cdef Py_ssize_t i, start, stop, step, slicelength
        cdef unsigned char * data
        cdef size_t j,k
        cdef int res
        if PySlice_Check(<PyObject *>key):
            # here the key is a slice
            res = PySlice_GetIndicesEx(<PySliceObject *>key,
                    self._length,
                    &start, &stop, &step,
                    &slicelength)
            if res < 0:
                print "HAAAAAAAAAAA... there will soon be a crash"
            if step == 1:
                return self._new_c(self._data+start, stop-start, self)
            data = <unsigned char *> sage_malloc(slicelength * sizeof(unsigned char))
            j = 0
            for k in range(start,stop,step):
                data[j] = self._data[k]
                j += 1
            return self._new_c(data, slicelength, None)

        elif PyIndex_Check(<PyObject *>key):
            # here the key is an int
            i = PyNumber_AsSsize_t(<PyObject *>key, PyExc_IndexError)
            if i == -1 and PyErr_Occurred():
                print "HAAAAAAAAAAA... there will soon be a crash"
            if i < 0:
                i += self._length;
            if i < 0 or i >= self._length:
                raise IndexError("word index out of range")
            return self._data[i]

        raise TypeError("word indices must be integers")

    def __iter__(self):
        cdef size_t i
        for i in range(self._length):
            yield self._data[i]

    cdef _concatenate(self, Word_char other):
        cdef unsigned char * data
        data = <unsigned char *> sage_malloc((self._length + other._length) * sizeof(unsigned char))
        if data == NULL:
            raise MemoryError

        memcpy(data, self._data, self._length * sizeof(unsigned char))
        memcpy(data+self._length, other._data, other._length * sizeof(unsigned char))

        return self._new_c(data, self._length + other._length, None)

    def __mul__(self, other):
        r"""
        Concatenation of ``self`` and ``other``.

        The result is automatically converted to a Word_char.
        """
        cdef Word_char w

        if PY_TYPE_CHECK(other, Word_char):
            return (<Word_char> self)._concatenate(other)

        elif PySequence_Check(<PyObject *>other):
            # we convert other to a Word_char and perform the concatenation
            w = (<Word_char> self)._new_c(NULL, 0, None)
            w._set_data(other)
            return (<Word_char> self)._concatenate(w)

        raise TypeError("not able to initialize a word from {}".format(other))

    def __len__(self):
        return self._length

    def length(self):
        return smallInteger(self._length)

    @cython.boundscheck(False)
    def has_prefix(self, other):
        r"""
        Test whether ``other`` is a prefix of ``self``.

        INPUT:

        - ``other`` -- a word or a sequence (e.g. tuple, list)
        """
        cdef size_t i
        cdef int test
        cdef unsigned char * data

        if PY_TYPE_CHECK(other, Word_char):
            # C level
            if (<Word_char> other)._length > self._length:
                return False
            return memcmp((<Word_char> self)._data,
                          (<Word_char> other)._data,
                          (<Word_char> other)._length) == 0

        elif PySequence_Check(<PyObject *>other):
            # python level
            if len(other) > self._length:
                return False

            for i in range(len(other)):
                if other[i] != self._data[i]:
                    return False
            return True

        raise TypeError("not able to initialize a word from {}".format(other))
