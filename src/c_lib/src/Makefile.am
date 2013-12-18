lib_LTLIBRARIES = libcsage.la

libcsage_la_LIBADD = -lntl -lpari -lmpir $(PYTHON_LDFLAGS)
libcsage_la_CPPFLAGS = -I$(top_srcdir)/include $(PYTHON_CPPFLAGS)

libcsage_la_SOURCES = \
	convert.c \
	gmp_globals.c \
	interrupt.c \
	memory.c \
	mpn_pylong.c \
	mpz_longlong.c \
	mpz_pylong.c \
	stdsage.c \
	ZZ_pylong.cpp \
	ntl_wrap.cpp
