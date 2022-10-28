from sage.structure.list_clone cimport ClonableIntArray

cdef list all_children(ClonableIntArray v, int max_part)
cdef int lex_cmp_partial(ClonableIntArray t1, ClonableIntArray t2, int step)
cdef int lex_cmp(ClonableIntArray t1, ClonableIntArray t2)
cdef bint is_canonical(list sgs, ClonableIntArray v) except -1
cdef ClonableIntArray canonical_representative_of_orbit_of(list sgs, ClonableIntArray v)
cdef list canonical_children(list sgs, ClonableIntArray v, int max_part)
cdef set orbit(list sgs, ClonableIntArray v)
