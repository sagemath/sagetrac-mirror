from .types cimport GEN
from .pari_instance cimport PariInstance

cdef void _pari_init_error_handling(PariInstance pari_instance)
cdef int _pari_err_handle(GEN E) except 0
cdef void _pari_err_recover(long errnum)
