# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
include $(top_srcdir)/common.mk
libs_cremonadir = $(pythondir)/sage/libs/cremona
libs_cremona_PYTHON = __init__.py constructor.py all.py
libs_cremona_LTLIBRARIES = homspace.la mat.la newforms.la

AM_DEFAULT_SOURCE_EXT = .pyxx
homspace_la_LIBADD = -lec -lntl -lpari @GMP_OR_MPIR_XXLIBS@ @GMP_OR_MPIR_LIBS@ -lm -lcsage
homspace_la_CPPFLAGS = -DNTL_ALL $(AM_CPPFLAGS)
mat_la_LIBADD = -lec -lntl -lpari @GMP_OR_MPIR_XXLIBS@ @GMP_OR_MPIR_LIBS@ -lm
mat_la_CPPFLAGS = -DNTL_ALL $(AM_CPPFLAGS)
newforms_la_LIBADD = -lec -lntl -lpari @GMP_OR_MPIR_XXLIBS@ @GMP_OR_MPIR_LIBS@ -lm -lcsage
newforms_la_CPPFLAGS = -DNTL_ALL $(AM_CPPFLAGS)
