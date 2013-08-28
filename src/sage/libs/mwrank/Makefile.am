# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
include $(top_srcdir)/common.mk
libs_mwrankdir = $(pythondir)/sage/libs/mwrank
libs_mwrank_PYTHON = __init__.py all.py interface.py

libs_mwrank_LTLIBRARIES = mwrank.la
mwrank_la_SOURCES = mwrank.pyxx wrap.cc
mwrank_la_LIBADD = -lec -lntl -lpari @GMP_OR_MPIR_LIBS@ @GMP_OR_MPIR_XXLIBS@ -lm -lcsage
mwrank_la_CPPFLAGS = -DNTL_ALL $(AM_CPPFLAGS)
