cdef inline type element_class_zero = pAdicLazyElement_zero
cdef inline type element_class_one = pAdicLazyElement_one
cdef inline type element_class_bound = pAdicLazyElement_bound
cdef inline type element_class_value = pAdicLazyElement_value
cdef inline type element_class_random = pAdicLazyElement_random
cdef inline type element_class_slice = pAdicLazyElement_slice
cdef inline type element_class_add = pAdicLazyElement_add
cdef inline type element_class_sub = pAdicLazyElement_sub
cdef inline type element_class_mul = pAdicLazyElement_mul
cdef inline type element_class_muldigit = pAdicLazyElement_muldigit
cdef inline type element_class_div = pAdicLazyElement_div
cdef inline type element_class_sqrt = pAdicLazyElement_sqrt
cdef inline type element_class_teichmuller = pAdicLazyElement_teichmuller
cdef inline type element_class_unknown = pAdicLazyElement_unknown

include "sage/libs/linkages/padics/lazy/flint.pxi"
include "lazy_template.pxi"
