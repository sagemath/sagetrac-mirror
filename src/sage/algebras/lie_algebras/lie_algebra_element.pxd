from sage.structure.element cimport Element
from sage.structure.element_wrapper cimport ElementWrapper
from sage.structure.sage_object cimport SageObject
from sage.modules.with_basis.indexed_element cimport IndexedFreeModuleElement

cdef class LieAlgebraElement(IndexedFreeModuleElement):
    cdef lift(self)

cdef class LieAlgebraElementWrapper(ElementWrapper):
    cdef _add_(self, right)
    cdef _sub_(self, right)

cdef class LieAlgebraMatrixWrapper(LieAlgebraElementWrapper):
    pass

cdef class LieSubalgebraElementWrapper(LieAlgebraElementWrapper):
    cdef dict _monomial_coefficients
    cdef dict monomial_coefficients(self, bint copy=*)

cdef class StructureCoefficientsElement(LieAlgebraMatrixWrapper):
    cdef bracket(self, right)
    cdef _bracket_(self, right)
    cdef to_vector(self, bint sparse=*)
    cdef dict monomial_coefficients(self, bint copy=*)
    #cdef lift(self)

cdef class UntwistedAffineLieAlgebraElement(Element):
    cdef dict _t_dict
    cdef _c_coeff
    cdef _d_coeff
    cdef long _hash

    cdef _add_(self, other)
    cdef _sub_(self, other)
    cdef _neg_(self)

    cdef dict t_dict(self)
    cdef c_coefficient(self)
    cdef d_coefficient(self)

    cdef bracket(self, y)
    cdef _bracket_(self, y)
    cdef canonical_derivation(self)
    cdef monomial_coefficients(self, bint copy=*)

cdef class LieObject(SageObject):
    cdef tuple _word
    cdef public tuple _index_word
    cdef tuple to_word(self)

cdef class LieGenerator(LieObject):
    cdef public str _name

cdef class LieBracket(LieObject):
    cdef public LieObject _left
    cdef public LieObject _right
    cdef long _hash

    cdef lift(self, dict UEA_gens_dict)

cdef class GradedLieBracket(LieBracket):
    cdef public _grade

cdef class LyndonBracket(GradedLieBracket):
    pass
