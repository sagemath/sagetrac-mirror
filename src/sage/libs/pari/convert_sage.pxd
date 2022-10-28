from cypari2.gen cimport Gen
from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational

cdef gen_to_sage(Gen z, locals=*)

cdef set_integer_from_gen(Integer self, Gen x)
cdef Gen new_gen_from_integer(Integer self)
cdef set_rational_from_gen(Rational self, Gen x)
cdef Gen new_gen_from_rational(Rational self)

cdef pari_is_prime(Integer p)
cdef pari_is_prime_power(Integer q, bint get_data)
cdef unsigned long pari_maxprime()
cdef list pari_prime_range(long c_start, long c_stop, bint py_ints=*)
