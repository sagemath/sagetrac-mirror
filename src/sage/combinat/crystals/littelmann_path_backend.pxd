cdef class LittelmannPath(object):
    cdef public list value

    cpdef copy(self)
    cpdef list endpoint(self)
    cdef compress(self)
    cdef split_step(self, int which_step, r)
    cpdef reflect_step(self, int which_step, int i, list root)
    cdef list _string_data(self, int i)
    cpdef LittelmannPath dualize(self, inplace=*)
    cpdef int epsilon(self, int i)
    cpdef int phi(self, int i)
    cpdef e(self, int i, list root, int power=*, to_string_end=*)
    cpdef f(self, int i, list root, int power=*, to_string_end=*)
    cpdef LittelmannPath s(self, int i, list root)

cdef add_inplace_lists(list l, list a)
cdef inplace_axpy(a, list x, list y)
cdef positively_parallel_weights(list v, list w)

