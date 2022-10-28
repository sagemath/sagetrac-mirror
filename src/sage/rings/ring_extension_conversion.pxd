from sage.rings.ring_extension cimport RingExtension_generic


cdef backend_parent(R)
cdef from_backend_parent(R, RingExtension_generic E)

cdef backend_element(x)
cdef from_backend_element(x, RingExtension_generic E)

cdef _backend_morphism(f)
cdef backend_morphism(f, forget=*)
cdef from_backend_morphism(f, RingExtension_generic E)

cdef to_backend(arg)
cdef from_backend(arg, E)


