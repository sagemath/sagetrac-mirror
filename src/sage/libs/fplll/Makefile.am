# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
include $(top_srcdir)/common.mk
libs_fpllldir = $(pythondir)/sage/libs/fplll
libs_fplll_PYTHON = __init__.py

AM_DEFAULT_SOURCE_EXT = .pyxx
libs_fplll_LTLIBRARIES = fplll.la

fplll_la_LIBADD = @GMP_OR_MPIR_LIBS@ -lmpfr -lfplll -lcsage
fplll_la_CPPFLAGS = @FPLLL_INCLUDES@ -DFPLLL_V3_COMPAT $(AM_CPPFLAGS)
