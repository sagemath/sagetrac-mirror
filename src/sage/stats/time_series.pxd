cdef class TimeSeries:
    cdef double* _values
    cdef Py_ssize_t _length
    cdef rescale(self, double s)
    cdef double sum(self)
