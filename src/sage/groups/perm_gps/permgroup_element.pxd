from sage.structure.element cimport MultiplicativeGroupElement, MonoidElement, Element
from sage.structure.list_clone cimport ClonableIntArray
from sage.rings.polynomial.polydict cimport ETuple
from sage.libs.gap.element cimport GapElement

cdef class PermutationGroupElement(MultiplicativeGroupElement):
    cdef int* perm
    cdef int n
    cdef int perm_buf[15] # to avoid malloc for small elements
    cdef GapElement _libgap
    cdef PermutationGroupElement _new_c(self)
    cdef _alloc(self, int)
    cdef _set_identity(self)
    cdef _set_list_images(self, v, bint convert)
    cdef _set_libgap(self, GapElement p)
    cdef _set_list_cycles(self, c, bint convert)
    cdef _set_string(self, str s)
    cdef _set_permutation_group_element(self, PermutationGroupElement p, bint convert)

    cdef _mul_(self, other)
    cdef PermutationGroupElement _generate_new(self, list new_list)
    cdef PermutationGroupElement _generate_new_GAP(self, old)
    cdef _gap_list(self)
    cdef domain(self)
    cdef public __custom_name
    cdef list _act_on_list_on_position(self, list x)
    cdef ClonableIntArray _act_on_array_on_position(self, ClonableIntArray x)
    cdef ETuple _act_on_etuple_on_position(self, ETuple x)
