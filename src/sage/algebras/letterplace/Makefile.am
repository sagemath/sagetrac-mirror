# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
include $(top_srcdir)/common.mk
algebras_letterplacedir = $(pythondir)/sage/algebras/letterplace
algebras_letterplace_PYTHON = __init__.py
algebras_letterplace_LTLIBRARIES = \
    free_algebra_element_letterplace.la free_algebra_letterplace.la letterplace_ideal.la

AM_DEFAULT_SOURCE_EXT = .pyxx
free_algebra_element_letterplace_la_LIBADD = -lm -lreadline -lsingular -lgivaro -lntl @GMP_OR_MPIR_XXLIBS@ @GMP_OR_MPIR_LIBS@
free_algebra_element_letterplace_la_CPPFLAGS = @SINGULAR_INCLUDES@ $(AM_CPPFLAGS)
free_algebra_letterplace_la_LIBADD = -lm -lreadline -lsingular -lgivaro -lntl @GMP_OR_MPIR_XXLIBS@ @GMP_OR_MPIR_LIBS@
free_algebra_letterplace_la_CPPFLAGS = @SINGULAR_INCLUDES@ $(AM_CPPFLAGS)
letterplace_ideal_la_LIBADD = -lm -lreadline -lsingular -lgivaro -lntl @GMP_OR_MPIR_XXLIBS@ @GMP_OR_MPIR_LIBS@
letterplace_ideal_la_CPPFLAGS = @SINGULAR_INCLUDES@ $(AM_CPPFLAGS)
