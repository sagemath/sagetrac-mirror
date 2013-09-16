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
     cdef cppclass series_carry_relaxed[C]:
          pass

ctypedef series[modular[modulus[int] ], series_carry_p_adic[series_carry_naive] ] mmx_padic_naive
ctypedef series[modular[modulus[int] ], series_carry_p_adic[series_carry_relaxed[series_carry_naive]] ] mmx_padic_relaxed
#ctypedef mmx_padic_naive mmx_padic_relaxed

# Common stuff
cdef extern from "numerix/modular.hpp" namespace "mmx":
     # Cython does not support templated functions
     cdef void mmx_set_modulus "mmx::modular<mmx::modulus<int> >::set_modulus" (int)

# Naive implementation
cdef extern from "algebramix/p_adic.hpp" namespace "mmx":
     # Cython does not support templated functions
     mmx_padic_naive itsc_naive "mmx::integer_to_series_carry<mmx::modular<mmx::modulus<int> >, mmx::series_carry_p_adic<mmx::series_carry_naive> >" (integer)

     mmx_padic_naive operator + (const mmx_padic_naive& f, const mmx_padic_naive& g)
     mmx_padic_naive operator - (const mmx_padic_naive& f, const mmx_padic_naive& g)
     mmx_padic_naive operator * (const mmx_padic_naive& f, const mmx_padic_naive& g)
     mmx_padic_naive operator / (const mmx_padic_naive& f, const mmx_padic_naive& g)

     syntactic flatten_naive "mmx::flatten" (const mmx_padic_naive& f)
     port operator << (const port& out, const mmx_padic_naive& x)

# Relaxed implementation
cdef extern from "algebramix/p_adic.hpp" namespace "mmx":
     # Cython does not support templated functions
     mmx_padic_relaxed itsc_relaxed "mmx::integer_to_series_carry<mmx::modular<mmx::modulus<int> >, mmx::series_carry_p_adic<mmx::series_carry_relaxed<mmx::series_carry_naive> > >" (integer)

     mmx_padic_relaxed operator + (const mmx_padic_relaxed& f, const mmx_padic_relaxed& g)
     mmx_padic_relaxed operator - (const mmx_padic_relaxed& f, const mmx_padic_relaxed& g)
     mmx_padic_relaxed operator * (const mmx_padic_relaxed& f, const mmx_padic_relaxed& g)
     mmx_padic_relaxed operator / (const mmx_padic_relaxed& f, const mmx_padic_relaxed& g)

     syntactic flatten_relaxed "mmx::flatten" (const mmx_padic_relaxed& f)
     port operator << (const port& out, const mmx_padic_relaxed& x)

cdef class MMXpadic_naive:
     cdef mmx_padic_naive* x
     cdef int p

cdef class MMXpadic_relaxed:
     cdef mmx_padic_relaxed* x
     cdef int p
