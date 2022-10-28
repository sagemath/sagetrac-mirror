cdef class PolyDict:
    cdef dict __repn
    cdef object __zero

cdef class ETuple:
    cdef size_t _length
    cdef size_t _nonzero
    cdef int *_data

    cdef size_t unweighted_degree(self)
    cdef size_t weighted_degree(self, tuple w)
    cdef size_t unweighted_quotient_degree(self, ETuple other)
    cdef size_t weighted_quotient_degree(self, ETuple other, tuple w)
    cdef ETuple eadd(ETuple self, ETuple other)
    cdef ETuple esub(ETuple self, ETuple other)
    cdef ETuple emul(ETuple self, int factor)
    cdef ETuple emin(ETuple self, ETuple other)
    cdef ETuple emax(ETuple self, ETuple other)
    cdef ETuple eadd_p(ETuple self, int other, int pos)
    cdef ETuple eadd_scaled(ETuple self, ETuple other, int scalar)
    cdef int dotprod(ETuple self, ETuple other)
    cdef ETuple escalar_div(ETuple self, int n)
    cdef ETuple divide_by_gcd(self, ETuple other)
    cdef ETuple divide_by_var(self, size_t index)
    cdef bint divides(self, ETuple other)
    cdef bint is_constant(ETuple self)
    cdef bint is_multiple_of(ETuple self, int n)
    cdef list nonzero_positions(ETuple self, bint sort=*)
    cdef common_nonzero_positions(ETuple self, ETuple other, bint sort=*)
    cdef list nonzero_values(ETuple self, bint sort=*)
    cdef ETuple reversed(ETuple self)
    cdef ETuple _new(ETuple self)
    cdef size_t get_exp(ETuple self, int i)

