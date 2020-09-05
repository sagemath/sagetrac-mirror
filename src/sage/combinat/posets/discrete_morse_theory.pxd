
cdef list minimal_skipped_intervals(list L, Py_ssize_t c, bint pure)
cdef bint is_msi(list L, Py_ssize_t c, Py_ssize_t s, Py_ssize_t t, bint pure)
cdef bint interval_removed(list A, list B)
cdef list J_intervals(list L, list M)

