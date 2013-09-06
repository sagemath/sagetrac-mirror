from common cimport *
from integer cimport *

cdef extern from "numerix/modulus.hpp" namespace "mmx":
     cdef cppclass modulus_int_preinverse[n]:
          pass
     cdef cppclass modulus[C]:
          pass

cdef extern from "numerix/modular.hpp" namespace "mmx":
     cdef cppclass modular[C]:
          pass

cdef extern from "algebramix/series.hpp" namespace "mmx":
     cdef cppclass series[C, V]:
          series()
          series(const C& c)

cdef extern from "algebramix/series.hpp" namespace "mmx":
     cdef cppclass series_carry_variant_helper[C]:
          pass
     cdef cppclass series_carry_p_adic[T]:
          pass
     cdef cppclass series_carry_naive:
          pass

ctypedef series[modular[modulus[int] ], series_carry_p_adic[series_carry_naive] ] mmx_padic

cdef extern from "algebramix/p_adic.hpp" namespace "mmx":
     cdef void mmx_set_modulus "mmx::modular<mmx::modulus<int> >::set_modulus" (int)
     mmx_padic* create_padic "new mmx::series<mmx::modular<mmx::modulus<int> >, mmx::series_carry_p_adic<mmx::series_carry_naive> >" (const mmx_padic& C)
     mmx_padic itsc "mmx::integer_to_series_carry<mmx::modular<mmx::modulus<int> >, mmx::series_carry_p_adic<mmx::series_carry_naive> >" (integer)
     mmx_padic operator + (const mmx_padic& f, const mmx_padic& g)
     mmx_padic operator - (const mmx_padic& f, const mmx_padic& g)
     mmx_padic operator * (const mmx_padic& f, const mmx_padic& g)
     mmx_padic operator / (const mmx_padic& f, const mmx_padic& g)

     syntactic flatten (const mmx_padic& f)
     port operator << (const port& out, const mmx_padic& x)

cdef class MMXpadic:
     cdef mmx_padic* x
     cdef int p