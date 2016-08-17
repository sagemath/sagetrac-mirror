from skew_polynomial_element cimport SkewPolynomial_generic_dense
from sage.matrix.matrix_dense cimport Matrix_dense

cdef class SkewPolynomial_finite_field_dense (SkewPolynomial_generic_dense):

    cdef SkewPolynomial_finite_field_dense _rgcd(self,SkewPolynomial_finite_field_dense other)
    cdef void _inplace_lrem(self, SkewPolynomial_finite_field_dense other)
    cdef void _inplace_rrem(self, SkewPolynomial_finite_field_dense other)
    cdef void _inplace_lfloordiv(self, SkewPolynomial_finite_field_dense other)
    cdef void _inplace_rfloordiv(self, SkewPolynomial_finite_field_dense other)
    cdef void _inplace_lmonic(self)
    cdef void _inplace_rmonic(self)
    cdef void _inplace_rgcd(self,SkewPolynomial_finite_field_dense other)
    cdef Py_ssize_t _val_inplace_unit(self)
    cdef SkewPolynomial_finite_field_dense _rquo_inplace_rem(self, SkewPolynomial_finite_field_dense other)

    cdef Matrix_dense _matmul_c(self)

    # Karatsuba
    #cpdef RingElement _mul_karatsuba(self, RingElement other, cutoff=*)
    cpdef SkewPolynomial_finite_field_dense _mul_central(self, SkewPolynomial_finite_field_dense right)
    cpdef RingElement _mul_(self, RingElement right)
    cpdef rquo_rem_karatsuba(self, RingElement other, cutoff=*)

 cdef class SkewPolynomial_finite_field_karatsuba:
    cdef _parent
    cdef Py_ssize_t _order
    cdef Py_ssize_t _cutoff
    cdef RingElement _zero
    cdef _twist
    cdef char _algo_matrix
    cdef RingElement _t
    cdef Matrix_dense _T, _Tinv

    cdef list mul_step (self, list x, list y)
    cdef list mul_step_matrix(self, list x, list y)
    cdef list mul_iter(self, list x, list y, char flag)
    cdef list _twinv
    cdef list div_step(self, list a, Py_ssize_t ia, Py_ssize_t da, list b, Py_ssize_t ib, Py_ssize_t db)
    cdef list div_iter(self, list a, Py_ssize_t ia, Py_ssize_t da, list b, Py_ssize_t ib, Py_ssize_t db)
