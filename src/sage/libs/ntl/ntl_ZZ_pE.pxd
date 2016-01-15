from .types cimport ZZ_pE_c
from .ntl_ZZ_pEContext cimport ntl_ZZ_pEContext_class
from .ntl_ZZ_pX cimport ntl_ZZ_pX

cdef class ntl_ZZ_pE(object):
    cdef ZZ_pE_c x
    cdef ntl_ZZ_pEContext_class c
    cdef public ntl_ZZ_pX get_as_ZZ_pX(ntl_ZZ_pE self)
    cdef public void set_from_ZZ_pX(ntl_ZZ_pE self, ntl_ZZ_pX value)
    cdef ntl_ZZ_pE _new(self)
