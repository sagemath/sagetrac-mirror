# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
include $(top_srcdir)/common.mk
AM_CPPFLAGS+= $(PYTHON_INCLUDES) @LIBCSAGE_INCLUDES@
matrixdir = $(pythondir)/sage/matrix
matrix_PYTHON = \
    all.py benchmark.py berlekamp_massey.py constructor.py docs.py __init__.py \
    matrix_integer_dense_hnf.py matrix_integer_dense_saturation.py matrix_misc.py \
    matrix_space.py operation_table.py symplectic_basis.py tests.py

AM_DEFAULT_SOURCE_EXT = .pyx
matrix_LTLIBRARIES = \
    action.la change_ring.la matrix.la matrix0.la matrix1.la matrix2.la \
    matrix_complex_double_dense.la matrix_cyclo_dense.la matrix_dense.la \
    matrix_double_dense.la matrix_generic_dense.la matrix_generic_sparse.la \
    matrix_integer_2x2.la matrix_integer_dense.la matrix_integer_sparse.la \
    matrix_mod2_dense.la matrix_mod2e_dense.la matrix_modn_dense.la \
    matrix_modn_dense_double.la matrix_modn_dense_float.la matrix_modn_sparse.la \
    matrix_mpolynomial_dense.la matrix_rational_dense.la matrix_rational_sparse.la \
    matrix_real_double_dense.la matrix_sparse.la matrix_symbolic_dense.la \
    matrix_window.la matrix_window_modn_dense.la misc.la strassen.la

action_la_CPPFLAGS = $(PYTHON_INCLUDES) -I$(C_LIB_INCLUDE) $(AM_CPPFLAGS)
change_ring_la_LIBADD = -l@BLAS@ -l@BLAS2@ @GMP_OR_MPIR_LIBS@
change_ring_la_CPPFLAGS = @NUMPY_INCLUDES@ $(PYTHON_INCLUDES) $(AM_CPPFLAGS)

matrix_la_CPPFLAGS = $(PYTHON_INCLUDES) $(AM_CPPFLAGS)
matrix_la_LIBADD = -lcsage

matrix0_la_LIBADD = -lcsage
matrix1_la_LIBADD = -lcsage
matrix2_la_LIBADD = -lcsage

matrix_dense_la_LIBADD = -lcsage

matrix_complex_double_dense_la_LIBADD = -l@BLAS@ -l@BLAS2@ -lcsage
matrix_complex_double_dense_la_CPPFLAGS = @NUMPY_INCLUDES@ $(AM_CPPFLAGS)

matrix_cyclo_dense_la_LIBADD = -lntl @GMP_OR_MPIR_LIBS@ -lcsage
matrix_cyclo_dense_la_SOURCES = matrix_cyclo_dense.pyxx
matrix_cyclo_dense.cc: @srcdir@/matrix_cyclo_dense.pyxx

matrix_double_dense_la_LIBADD = -l@BLAS@ -l@BLAS2@
matrix_double_dense_la_CPPFLAGS = @NUMPY_INCLUDES@ $(AM_CPPFLAGS)
matrix_generic_dense_la_LIBADD = -lcsage
matrix_integer_2x2_la_LIBADD = @GMP_OR_MPIR_LIBS@ -lcsage
matrix_integer_dense_la_LIBADD = -liml -lpari @GMP_OR_MPIR_LIBS@ -lm -l@BLAS@ -l@BLAS2@ -lcsage
matrix_integer_dense_la_CFLAGS = @M4RI_CFLAGS@ $(AM_CFLAGS)
matrix_integer_dense_la_CPPFLAGS = -I@SAGE_DIR@/ext $(AM_CPPFLAGS)
matrix_integer_sparse_la_LIBADD = @GMP_OR_MPIR_LIBS@ -lcsage
matrix_mod2_dense_la_LIBADD = @GMP_OR_MPIR_LIBS@ -lm4ri -lgd -lpng12 -lz -lcsage
matrix_mod2_dense_la_CFLAGS = @M4RI_CFLAGS@ $(AM_CFLAGS)
matrix_mod2e_dense_la_LIBADD = -lm4rie -lm4ri -lm -lcsage
matrix_mod2e_dense_la_CFLAGS = @M4RI_CFLAGS@ $(AM_CFLAGS)
matrix_mod2e_dense_la_CPPFLAGS = @LIBM4RIE_INCLUDES@ $(AM_CPPFLAGS)
matrix_modn_dense_la_LIBADD = @GMP_OR_MPIR_LIBS@ -lcsage

matrix_modn_dense_double_la_SOURCES = matrix_modn_dense_double.pyxx
matrix_modn_dense_double_la_LIBADD = -llinbox -lgivaro -lmpfr @GMP_OR_MPIR_XXLIBS@ @GMP_OR_MPIR_LIBS@ -l@BLAS@ -l@BLAS2@ -lcsage
matrix_modn_dense_double.cc: @srcdir@/matrix_modn_dense_double.pyxx
matrix_modn_dense_double_la_CPPFLAGS = -DDISABLE_COMMENTATOR -D__STDC_LIMIT_MACROS $(AM_CXXFLAGS) \
                                       -DINTEGER_MOD_INT32_LIMIT=46341 \
                                       -DINTEGER_MOD_INT64_LIMIT=2147483647 $(AM_CPPFLAGS)

matrix_modn_dense_float_la_SOURCES = matrix_modn_dense_float.pyxx
matrix_modn_dense_float_la_LIBADD = -llinbox -lgivaro -lmpfr @GMP_OR_MPIR_XXLIBS@ @GMP_OR_MPIR_LIBS@ -l@BLAS@ -l@BLAS2@ -lcsage
matrix_modn_dense_float.cc: @srcdir@/matrix_modn_dense_float.pyxx
matrix_modn_dense_float_la_CPPFLAGS = -DDISABLE_COMMENTATOR -D__STDC_LIMIT_MACROS $(AM_CPPFLAGS)
matrix_modn_sparse_la_LIBADD = @GMP_OR_MPIR_LIBS@ -lcsage

matrix_mpolynomial_dense_la_SOURCES = matrix_mpolynomial_dense.pyxx
matrix_mpolynomial_dense_la_LIBADD = -lm -lreadline -lsingular -lgivaro -lntl @GMP_OR_MPIR_XXLIBS@ @GMP_OR_MPIR_LIBS@ -lcsage
matrix_mpolynomial_dense_la_CPPFLAGS = $(AM_CPPFLAGS) @SINGULAR_INCLUDES@ \
                                       -D__STDC_LIMIT_MACROS
matrix_mpolynomial_dense.cc: @srcdir@/matrix_mpolynomial_dense.pyxx

matrix_rational_dense_la_LIBADD = -lpari @GMP_OR_MPIR_LIBS@ -lcsage
matrix_rational_sparse_la_LIBADD = @GMP_OR_MPIR_LIBS@ -lcsage
matrix_real_double_dense_la_LIBADD = -l@BLAS@ -l@BLAS2@ -lcsage
matrix_real_double_dense_la_CPPFLAGS = @NUMPY_INCLUDES@ $(AM_CPPFLAGS)

matrix_sparse_la_LIBADD = -lcsage -lcsage

matrix_window_modn_dense_la_LIBADD = -lcsage
matrix_window_la_LIBADD = -lcsage

misc_la_LIBADD = -lmpfr @GMP_OR_MPIR_LIBS@ -lcsage

# somethings terribly wrong with singular
matrix_mpolynomial_dense_la_CPPFLAGS += -DNDEBUG

noinst_HEADERS = matrix_modn_dense_template_header.pxi matrix_modn_dense_template.pxi
