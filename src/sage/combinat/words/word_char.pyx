include 'sage/ext/interrupt.pxi'
include 'sage/ext/stdsage.pxi'

from sage.rings.integer cimport Integer, smallInteger

from libc.string cimport memcpy, memcmp

cdef class Word_char(object):
    cdef unsigned char * _data
    cdef size_t _length
    cdef Word_char _master
    cdef public _parent

    def __init__(self, parent, data):
        self._parent = parent
        self._master = None
        self._data = NULL
        self._length = 0

        if data:
            self._set_data(data)

    cdef _set_data(self, data):
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

    def __richcmp__(self, other, op):
        if not PY_TYPE_CHECK(other, Word_char):
            return NotImplemented

        cdef int test = (<Word_char> self)._cmp_c_impl(<Word_char> other)
        # 0: <
        # 1: <=
        # 2: ==
        # 3: !=
        # 4: >
        # 5: >=
        if test < 0:
            return op < 4
        elif test > 0:
            return op > 3
        else:
            return op == 1 or op == 2 or op == 5

    cdef int _cmp_c_impl(self, Word_char other) except -2451:

        cdef size_t l = min(self._length, other._length)

        test = memcmp(<void *> (<Word_char> self)._data,
                      <void *> (<Word_char> other)._data,
                      l * sizeof(unsigned char))

        if test == 0:
            if self._length > other._length:
                return 1
            elif self._length < other._length:
                return -1
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
        cdef long i, start, stop, step
        cdef unsigned char * data
        cdef size_t j,k
        if isinstance(key, slice):
            start, stop, step = key.indices(self._length)
            if step == 1:
                return self._new_c(self._data+start, stop-start, self)
            elif step > 0:
                length = (stop-start)/step
            else:
                length = (start-stop-step-1)/(-step) # ceil (stop-start)/step
                if length < 0:
                    length = 0
            data = <unsigned char *> sage_malloc(length * sizeof(unsigned char))
            j = 0
            for k in range(start,stop,step):
                data[j] = self._data[k]
                j += 1
            return self._new_c(data, length, None)

        else:
            # here the key is an int
            i = key.__index__()
            if 0 <= i < self._length:
                return self._data[i]
            elif 0 < -i <= self._length:
                return self._data[self._length+i]
            raise IndexError("word index out of range")

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
        else:
            # we convert other to a Word_char and perform the concatenation
            w = (<Word_char> self)._new_c(NULL, 0, None)
            w._set_data(other)
            return (<Word_char> self)._concatenate(w)

    def __len__(self):
        return self._length

    def length(self):
        return smallInteger(self._length)
    
    def index(self, unsigned char letter, size_t start=0):
        cdef size_t i
        for i in range(start, self._length):
            if self._data[i] == letter:
                return i
        return -1

    def has_prefix(self, other):
        r"""
        Test whether ``other`` is a prefix of ``self``.
        """
        cdef size_t i
        cdef unsigned char * data

        if PY_TYPE_CHECK(other, Word_char):
            # C level
            if (<Word_char> other)._length > self._length:
                return False
            data = (<Word_char> other)._data
            for i in range((<Word_char> other)._length):
                if data[i] != self._data[i]:
                    return False
            return True

        else:
            # python level
            if len(other) > self._length:
                return False

            for i in range(len(other)):
                if other[i] != self._data[i]:
                    return False
            return True

