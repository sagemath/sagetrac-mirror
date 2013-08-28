SUBDIRS =
include $(top_srcdir)/common.mk
modular_modformdir = $(pythondir)/sage/modular/modform
modular_modform_PYTHON = \
space.py cuspidal_submodule.py __init__.py hecke_operator_on_qexp.py \
find_generators.py defaults.py vm_basis.py j_invariant.py tests.py \
ambient_g0.py test.py ambient.py submodule.py ambient_g1.py element.py \
constructor.py all.py theta.py half_integral.py eisenstein_submodule.py \
eis_series.py notes.py ambient_eps.py periods.py ambient_R.py numerical.py

AM_DEFAULT_SOURCE_EXT = .pyx
modular_modform_LTLIBRARIES = eis_series_cython.la
eis_series_cython_la_LIBADD = @GMP_OR_MPIR_LIBS@ -lflint -lcsage
eis_series_cython_la_CFLAGS = -std=c99 $(AM_CFLAGS)
