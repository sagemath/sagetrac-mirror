from .types cimport mat_ZZ_p_c
from .ntl_ZZ_pContext cimport ntl_ZZ_pContext_class
from .ntl_ZZ_p cimport ntl_ZZ_p

cdef class ntl_mat_ZZ_p(object):
    cdef mat_ZZ_p_c x
    cdef ntl_ZZ_pContext_class c
    cdef ntl_ZZ_p _new_element(self)
    cdef long __nrows, __ncols
    cdef ntl_mat_ZZ_p _new(self)
