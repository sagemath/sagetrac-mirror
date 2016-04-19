from sage.libs.pari.gen cimport gen
from .pari_instance cimport PariInstance

cpdef gen objtoclosure(PariInstance pari, f)
