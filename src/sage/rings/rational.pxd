from sage.libs.gmp.types cimport mpq_t

cimport sage.structure.element
cimport sage.rings.integer as integer

cdef rational_power_parts(a, Rational b, factor_limit=?)

cdef class Rational(sage.structure.element.FieldElement):
    cdef mpq_t value

    cdef _add_(self, other)
    cdef _mul_(self, other)
    cdef _pow_(self, other)
    cdef __set_value(self, x, unsigned int base)
    cdef void set_from_mpq(Rational self, mpq_t value)
    cdef _lshift(self, long int exp)
    cdef _rshift(self, long int exp)

    cdef _val_unit(self, integer.Integer p)
