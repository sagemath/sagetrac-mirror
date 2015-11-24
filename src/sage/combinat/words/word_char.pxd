from sage.combinat.words.word_datatypes cimport WordDatatype

cdef class WordDatatype_char(WordDatatype):
    # the begining of the word is always _data[0] + _start
    cdef unsigned char ** _data
    cdef size_t _start

    # _master is a just a reference to another Python object in case the finite
    # word is just a slice of another one. But because Cython takes care of
    # Python attributes *before* the call to __dealloc__ we need to duplicate
    # the information.
    cdef int _is_slice
    cdef WordDatatype_char _master

    # for words that are lazy (like infinite ones), this method guarantee that
    # the word is computed to at least length n.
    cdef _update_word_up_to(self, size_t n)

cdef class WordDatatype_char_finite(WordDatatype_char):
    cdef size_t _stop    # end of the word

    cdef _set_data(self, data)
    cdef int _lexico_cmp(self, WordDatatype_char_finite other) except -2

    cdef WordDatatype_char_finite _new_c(self, unsigned char ** data, start, stop, WordDatatype_char_finite master)
    cdef _concatenate(self, WordDatatype_char_finite other)

cdef class WordDatatype_char_infinite(WordDatatype_char):
    cdef size_t _current_stop  # current position up to which the word is known
    cdef size_t _alloc_size

    cdef int _init_c_data(self, size_t alloc_size) except -1
    cdef int _realloc(self, size_t new_size) except -1
    cdef int _extend_word(self) except -1
    cdef _update_word_up_to(self, size_t n)

    cpdef WordDatatype_char_finite _finite_slice(self, Py_ssize_t start, Py_ssize_t stop, Py_ssize_t step, Py_ssize_t slicelength)
    cpdef WordDatatype_char_infinite _infinite_slice(self, size_t start, size_t step)

cdef class WordDatatype_callable_char(WordDatatype_char_infinite):
    cdef _f
    cdef size_t _slice_length

cdef class WordDatatype_iter_char(WordDatatype_char_infinite):
    cdef _iterator
    cdef size_t _slice_length

cdef class WordDatatype_substitutive_char(WordDatatype_char_infinite):
    cdef unsigned char ** _images    # images as pointers to char*
    cdef int * _lengths              # length of the images
    cdef size_t _i                   # current index

