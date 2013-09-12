# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
c_libdir = $(libdir)
c_lib_LTLIBRARIES = libcsage.la
libcsage_la_LIBADD = -lntl -lpari @GMP_OR_MPIR_LIBS@ -lpython@PYTHON_VERSION@
libcsage_la_CPPFLAGS = $(AM_CPPFLAGS) -I$(top_srcdir)/include @PYTHON_INCLUDES@
libcsage_la_SOURCES = \
    convert.c \
    gmp_globals.c \
    interrupt.c \
    memory.c \
    mpn_pylong.c \
    mpz_longlong.c \
    mpz_pylong.c \
    stdsage.c \
    ntl_wrap.cpp \
    ZZ_pylong.cpp
