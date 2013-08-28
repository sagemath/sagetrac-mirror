# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
include $(top_srcdir)/common.mk
libs_linboxdir = $(pythondir)/sage/libs/linbox
libs_linbox_PYTHON = __init__.py

AM_DEFAULT_SOURCE_EXT = .pyxx
libs_linbox_LTLIBRARIES = linbox.la
linbox_la_LIBADD = -llinboxsage -lntl -liml -llinbox -lstdc++ -lgivaro -lmpfr \
                   @GMP_OR_MPIR_LIBS@ @GMP_OR_MPIR_XXLIBS@ -l@BLAS@ -l@BLAS2@ -lcsage
linbox_la_CXXFLAGS = -D__STDC_LIMIT_MACROS $(AM_CXXFLAGS)
