# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
include $(top_srcdir)/common.mk
libs_symmetricadir = $(pythondir)/sage/libs/symmetrica
libs_symmetrica_PYTHON = __init__.py all.py

AM_DEFAULT_SOURCE_EXT = .pyx
libs_symmetrica_LTLIBRARIES = symmetrica.la

symmetrica_la_LIBADD = -lsymmetrica -lcsage
symmetrica_la_CPPFLAGS = $(AM_CPPFLAGS) -I/usr/include/malloc
