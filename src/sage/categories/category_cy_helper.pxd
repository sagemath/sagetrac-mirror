cdef tuple _sort_uniq(categories)
cdef class AxiomContainer(dict):
     pass
cdef tuple canonicalize_axioms(AxiomContainer all_axioms, axioms)
from sage.misc.classcall_metaclass cimport ClasscallMetaclass
cdef tuple _flatten_categories(categories, ClasscallMetaclass JoinCategory)
cdef tuple join_as_tuple(tuple categories, tuple axioms, tuple ignore_axioms)
