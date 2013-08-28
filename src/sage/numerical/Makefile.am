# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
SUBDIRS = backends
include $(top_srcdir)/common.mk
numericaldir = $(pythondir)/sage/numerical
numerical_PYTHON = __init__.py all.py knapsack.py optimize.py test.py

AM_DEFAULT_SOURCE_EXT = .pyx
numerical_LTLIBRARIES = linear_functions.la mip.la
linear_functions_la_LIBADD = -lcsage -lstdc++
linear_functions_la_CPPFLAGS = -I@SAGE_INC@ $(AM_CPPFLAGS)
mip_la_LIBADD = -lcsage -lstdc++
mip_la_CPPFLAGS = -I@SAGE_INC@ $(AM_CPPFLAGS)
