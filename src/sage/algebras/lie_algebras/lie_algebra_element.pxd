from sage.structure.element cimport Element
from sage.structure.element_wrapper cimport ElementWrapper

cdef class LieAlgebraElementWrapper(ElementWrapper):
    cpdef _add_(self, right)
    cpdef _sub_(self, right)

cdef class LieAlgebraMatrixWrapper(LieAlgebraElementWrapper):
    pass

cdef class StructureCoefficientsElement(LieAlgebraMatrixWrapper):
    cpdef bracket(self, right)
    cpdef _bracket_(self, right)
    cpdef to_vector(self)
    cpdef dict monomial_coefficients(self, bint copy=*)
    #cpdef lift(self)

cdef class UntwistedAffineLieAlgebraElement(Element):
    cdef dict _t_dict
    cdef _c_coeff
    cdef _d_coeff
    cdef long _hash

    cpdef dict t_dict(self)
    cpdef c_coefficient(self)
    cpdef d_coefficient(self)

    cpdef bracket(self, y)
    cpdef _bracket_(self, y)
    cpdef canonical_derivation(self)
    cpdef monomial_coefficients(self, bint copy=*)

