# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
include $(top_srcdir)/common.mk
libs_paridir = $(pythondir)/sage/libs/pari
libs_pari_PYTHON = __init__.py all.py gen_py.py

AM_DEFAULT_SOURCE_EXT = .pyx
libs_pari_LTLIBRARIES = gen.la
gen_la_LIBADD = -lpari @GMP_OR_MPIR_LIBS@ -lcsage
