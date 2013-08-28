# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
include $(top_srcdir)/common.mk
libs_flintdir = $(pythondir)/sage/libs/flint
libs_flint_PYTHON = __init__.py

AM_DEFAULT_SOURCE_EXT = .pyx
libs_flint_LTLIBRARIES = arith.la flint.la fmpz_poly.la

arith_la_LIBADD = -lcsage -lflint @GMP_OR_MPIR_LIBS@ @GMP_OR_MPIR_XXLIBS@ -lm
arith_la_CFLAGS = -std=c99 -D_XPG6 $(AM_CFLAGS)

flint_la_LIBADD = -lcsage -lflint @GMP_OR_MPIR_LIBS@ @GMP_OR_MPIR_XXLIBS@ -lm
flint_la_CFLAGS = -std=c99 -D_XPG6 $(AM_CFLAGS)

fmpz_poly_la_LIBADD = -lcsage -lflint @GMP_OR_MPIR_LIBS@ @GMP_OR_MPIR_XXLIBS@ -lm
fmpz_poly_la_CFLAGS = -std=c99 -D_XPG6 $(AM_CFLAGS)
