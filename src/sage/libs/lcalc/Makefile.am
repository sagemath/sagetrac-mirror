# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
include $(top_srcdir)/common.mk
libs_lcalcdir = $(pythondir)/sage/libs/lcalc
libs_lcalc_PYTHON = __init__.py

AM_DEFAULT_SOURCE_EXT = .pyxx
libs_lcalc_LTLIBRARIES = lcalc_Lfunction.la

lcalc_Lfunction_la_LIBADD = -lm -lntl -lmpfr @GMP_OR_MPIR_LIBS@ @GMP_OR_MPIR_XXLIBS@ -lLfunction -lstdc++
lcalc_Lfunction_la_CXXFLAGS = -O3 -ffast-math $(AM_CXXFLAGS)
lcalc_Lfunction_la_CPPFLAGS = @LIBLFUNCTION_INCLUDES@ $(AM_CPPFLAGS)
