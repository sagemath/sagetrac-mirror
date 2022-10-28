from cypari2.gen cimport Gen
from sage.rings.real_double cimport RealDoubleElement

cdef Gen new_gen_from_real_double_element(RealDoubleElement self)
