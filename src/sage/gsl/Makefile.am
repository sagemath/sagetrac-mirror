# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
include $(top_srcdir)/common.mk
gsldir = $(pythondir)/sage/gsl
gsl_PYTHON = __init__.py dft.py all.py

noinst_HEADERS = \
    gsl_airy.pxi gsl_bessel.pxi gsl_blas.pxi gsl_blas_types.pxi gsl_block.pxi \
    gsl_chebyshev.pxi gsl_clausen.pxi gsl_combination.pxi gsl_complex.pxi \
    gsl_coulomb.pxi gsl_coupling.pxi gsl_dawson.pxi gsl_debye.pxi gsl_diff.pxi \
    gsl_dilog.pxi gsl_eigen.pxi gsl_elementary.pxi gsl_ellint.pxi gsl_elljac.pxi \
    gsl_erf.pxi gsl_errno.pxi gsl_expint.pxi gsl_exp.pxi gsl_fermi_dirac.pxi \
    gsl_fft.pxi gsl_fit.pxi gsl_gamma.pxi gsl_gegenbauer.pxi gsl_histogram.pxi \
    gsl_hyperg.pxi gsl_integration.pxi gsl_interp.pxi gsl_laguerre.pxi \
    gsl_lambert.pxi gsl_legendre.pxi gsl_linalg.pxi gsl_log.pxi gsl_math.pxi \
    gsl_matrix_complex.pxi gsl_matrix.pxi gsl_min.pxi gsl_mode.pxi gsl_monte.pxi \
    gsl_ntuple.pxi gsl_odeiv.pxi gsl_permutation.pxi gsl_poly.pxi gsl_pow_int.pxi \
    gsl_psi.pxi gsl.pxi gsl_qrng.pxi gsl_random.pxi gsl_rng.pxi gsl_roots.pxi \
    gsl_sf_result.pxi gsl_sort.pxi gsl_statistics.pxi gsl_sum.pxi \
    gsl_synchrotron.pxi gsl_transport.pxi gsl_trig.pxi gsl_vector_complex.pxi \
    gsl_vector.pxi gsl_wavelet.pxi gsl_zeta.pxi math.pxi

AM_DEFAULT_SOURCE_EXT = .pyx
gsl_LTLIBRARIES = \
	 callback.la dwt.la fft.la gsl_array.la integration.la interpolation.la \
	 ode.la probability_distribution.la

callback_la_LIBADD = -lgsl -l@BLAS@ -l@BLAS2@
callback_la_CPPFLAGS = -DGSL_DISABLE_DEPRECATED=1 $(AM_CPPFLAGS)
dwt_la_LIBADD = -lgsl -l@BLAS@
dwt_la_CPPFLAGS = -DGSL_DISABLE_DEPRECATED=1 $(AM_CPPFLAGS)
fft_la_LIBADD = -lgsl -l@BLAS@ -l@BLAS2@ -lcsage
fft_la_CPPFLAGS = -DGSL_DISABLE_DEPRECATED=1 $(AM_CPPFLAGS)
gsl_array_la_LIBADD = -lgsl -l@BLAS@ -l@BLAS2@ -lcsage
gsl_array_la_CPPFLAGS = -DGSL_DISABLE_DEPRECATED=1 $(AM_CPPFLAGS)
integration_la_LIBADD = -lgsl -l@BLAS@ -l@BLAS2@ -lcsage
integration_la_CPPFLAGS = -DGSL_DISABLE_DEPRECATED=1 $(AM_CPPFLAGS)
interpolation_la_LIBADD = -lgsl -l@BLAS@ -l@BLAS2@ -lcsage
interpolation_la_CPPFLAGS = -DGSL_DISABLE_DEPRECATED=1 $(AM_CPPFLAGS)
ode_la_LIBADD = -lgsl -l@BLAS@ -l@BLAS2@ -lcsage
ode_la_CPPFLAGS = -DGSL_DISABLE_DEPRECATED=1 $(AM_CPPFLAGS)
probability_distribution_la_LIBADD = -lgsl -l@BLAS@ -l@BLAS2@ -lcsage
probability_distribution_la_CPPFLAGS = -DGSL_DISABLE_DEPRECATED=1 $(AM_CPPFLAGS)
