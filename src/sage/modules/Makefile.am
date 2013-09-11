# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
SUBDIRS = fg_pid
include $(top_srcdir)/common.mk
modulesdir = $(pythondir)/sage/modules
modules_PYTHON = \
    misc.py __init__.py real_double_vector.py free_module_morphism.py \
    vector_space_morphism.py vector_callable_symbolic_dense.py \
    vector_symbolic_dense.py vector_space_homspace.py complex_double_vector.py \
    all.py matrix_morphism.py free_module.py quotient_module.py module_element.py \
    free_quadratic_module.py free_module_homspace.py

noinst_HEADERS = \
	 binary_search.pxi vector_integer_sparse_c.pxi vector_integer_sparse_h.pxi \
    vector_modn_sparse_c.pxi vector_modn_sparse_h.pxi vector_rational_sparse_c.pxi \
    vector_rational_sparse_h.pxi

AM_DEFAULT_SOURCE_EXT = .pyx
modules_LTLIBRARIES = \
    finite_submodule_iter.la free_module_element.la module.la \
    vector_complex_double_dense.la vector_double_dense.la vector_integer_dense.la \
    vector_mod2_dense.la vector_modn_dense.la vector_rational_dense.la \
    vector_real_double_dense.la

vector_complex_double_dense_la_LIBADD = -l@BLAS@ -l@BLAS2@ -lcsage
vector_complex_double_dense_la_CPPFLAGS = @NUMPY_INCLUDES@ $(AM_CPPFLAGS)
vector_double_dense_la_LIBADD = -l@BLAS@ -l@BLAS2@ -lcsage
vector_double_dense_la_CPPFLAGS = @NUMPY_INCLUDES@ $(AM_CPPFLAGS)
vector_integer_dense_la_LIBADD = @GMP_OR_MPIR_LIBS@ -lcsage
vector_mod2_dense_la_LIBADD = @GMP_OR_MPIR_LIBS@ -lm4ri -lpng12 -lgd -lcsage
vector_mod2_dense_la_CFLAGS = @M4RI_CFLAGS@ $(AM_CFLAGS)

vector_modn_dense_la_CFLAGS = $(AM_CFLAGS) -std=c99
vector_modn_dense_la_LIBADD = -lcsage

vector_rational_dense_la_LIBADD = @GMP_OR_MPIR_LIBS@ -lcsage

vector_real_double_dense_la_LIBADD = -l@BLAS@ -l@BLAS2@ -lcsage
vector_real_double_dense_la_CPPFLAGS = @NUMPY_INCLUDES@ $(AM_CPPFLAGS)

free_module_element_la_LIBADD = -lcsage
