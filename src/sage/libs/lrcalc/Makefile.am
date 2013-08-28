# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
include $(top_srcdir)/common.mk
libs_lrcalcdir = $(pythondir)/sage/libs/lrcalc
libs_lrcalc_PYTHON = __init__.py

AM_DEFAULT_SOURCE_EXT = .pyx
libs_lrcalc_LTLIBRARIES = lrcalc.la
lrcalc_la_LIBADD = -llrcalc
lrcalc_la_CPPFLAGS = @LIBLRCALC_INCLUDES@ $(AM_CPPFLAGS)
