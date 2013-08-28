# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
include $(top_srcdir)/common.mk
algebras_quatalgdir = $(pythondir)/sage/algebras/quatalg
algebras_quatalg_PYTHON = __init__.py all.py quaternion_algebra.py

AM_DEFAULT_SOURCE_EXT = .pyxx
algebras_quatalg_LTLIBRARIES = quaternion_algebra_cython.la quaternion_algebra_element.la

quaternion_algebra_cython_la_LIBADD = -lcsage -lflint @GMP_OR_MPIR_LIBS@ @GMP_OR_MPIR_XXLIBS@ -lm -lstdc++ -lntl
quaternion_algebra_cython_la_CPPFLAGS = @LIBFLINT_INCLUDES@ @PYTHON_INCLUDES@ $(AM_CPPFLAGS)

quaternion_algebra_element_la_LIBADD = -lcsage -lflint @GMP_OR_MPIR_LIBS@ @GMP_OR_MPIR_XXLIBS@ -lm -lstdc++ -lntl
