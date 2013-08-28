# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
SUBDIRS = \
    finite_rings function_field number_field padics polynomial \
    universal_cyclotomic_field semirings
include $(top_srcdir)/common.mk
ringsdir = $(pythondir)/sage/rings

rings_PYTHON = \
    misc.py infinity.py power_series_ring.py integral_domain_element.py __init__.py \
    invariant_theory.py complex_field.py multi_power_series_ring_element.py \
    complex_interval_field.py rational_field.py euclidean_domain_element.py \
    tests.py ring_element.py principal_ideal_domain.py ideal_monoid.py \
    integral_domain.py big_oh.py arith.py qqbar.py field.py pari_ring.py \
    all.py laurent_series_ring.py fraction_field.py \
    commutative_algebra.py quotient_ring_element.py commutative_ring.py \
    commutative_algebra_element.py multi_power_series_ring.py monomials.py ideal.py \
    integer_ring_python.py quotient_ring.py dedekind_domain.py \
    dedekind_domain_element.py homset.py euclidean_domain.py contfrac.py \
    field_element.py principal_ideal_domain_element.py commutative_ring_element.py \
    real_interval_field.py

AM_DEFAULT_SOURCE_EXT = .pyx
rings_LTLIBRARIES = \
    bernmm.la bernoulli_mod_p.la complex_double.la \
    complex_interval.la complex_mpc.la complex_number.la factorint.la fast_arith.la \
    fraction_field_FpT.la fraction_field_element.la integer.la integer_ring.la \
    laurent_series_ring_element.la morphism.la noncommutative_ideals.la \
    power_series_mpoly.la power_series_poly.la power_series_ring_element.la \
    rational.la real_double.la real_interval_absolute.la real_lazy.la real_mpfi.la \
    real_mpfr.la residue_field.la ring.la

bernmm_la_SOURCES = bernmm.pyxx bernmm/bern_modp.cpp bernmm/bern_modp_util.cpp bernmm/bern_rat.cpp
bernmm_la_LIBADD = -lntl @GMP_OR_MPIR_LIBS@ -lpthread -lcsage
bernmm_la_CPPFLAGS = -DUSE_THREADS=1 -DTHREAD_STACK_SIZE=4096 $(AM_CPPFLAGS)

bernoulli_mod_p_la_LIBADD = -lntl
bernoulli_mod_p_la_SOURCES = bernoulli_mod_p.pyxx
bernoulli_mod_p_la_CPPFLAGS = -I@SAGE_DIR@/libs/ntl/ $(AM_CPPFLAGS)

complex_double_la_LIBADD = -lgsl -l@BLAS@ -l@BLAS2@ -lpari @GMP_OR_MPIR_LIBS@ -lm -lcsage
complex_double_la_CFLAGS = -std=c99 -D_XPG6 $(AM_CFLAGS)

complex_interval_la_LIBADD = -lmpfi -lmpfr @GMP_OR_MPIR_LIBS@ -lcsage

complex_mpc_la_LIBADD = -lmpc -lmpfr @GMP_OR_MPIR_LIBS@ -lcsage

complex_number_la_LIBADD = -lmpfr @GMP_OR_MPIR_LIBS@ -lcsage
factorint_la_LIBADD = @GMP_OR_MPIR_LIBS@ -lcsage
fast_arith_la_LIBADD = -lpari @GMP_OR_MPIR_LIBS@ -lcsage

# fixme?
fraction_field_FpT_la_LIBADD = -lcsage -lflint @GMP_OR_MPIR_LIBS@ @GMP_OR_MPIR_XXLIBS@ -lntl -lzn_poly
fraction_field_FpT_la_SOURCES = fraction_field_FpT.pyxx

integer_la_LIBADD = -lntl -lpari -lflint @GMP_OR_MPIR_LIBS@ -lcsage
integer_la_CPPFLAGS = $(PYTHON_INCLUDES) $(AM_CPPFLAGS)

integer_ring_la_LIBADD = -lntl @GMP_OR_MPIR_LIBS@ -lcsage

laurent_series_ring_element_la_LIBADD = -lcsage

power_series_mpoly_la_LIBADD = -lcsage
power_series_poly_la_LIBADD = -lcsage
power_series_ring_element_la_LIBADD = -lcsage

residue_field_la_LIBADD = -lcsage

ring_la_LIBADD = -lntl @GMP_OR_MPIR_LIBS@ -lcsage

morphism_la_LIBADD = -lcsage
rational_la_LIBADD = -lntl @GMP_OR_MPIR_LIBS@ -lcsage

real_lazy_la_LIBADD = -lcsage

real_double_la_LIBADD = -lgsl @GMP_OR_MPIR_LIBS@ -l@BLAS@ -l@BLAS2@ -lcsage
real_double_la_CPPFLAGS = -DGSL_DISABLE_DEPRECATED=1 $(AM_CPPFLAGS)

real_interval_absolute_la_LIBADD = @GMP_OR_MPIR_LIBS@

real_mpfi_la_LIBADD = -lmpfi -lmpfr @GMP_OR_MPIR_LIBS@ -lcsage

real_mpfr_la_LIBADD = -lmpfr -lpari @GMP_OR_MPIR_LIBS@ -lcsage

noinst_HEADERS = mpfi.pxi mpc.pxi
