# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
if USE_LIBCOXETER3
COND_COXETER3 = coxeter3
endif
if USE_LIBLRCALC
COND_LRCALC = lrcalc
endif

if USE_LIBFES
FES_MODULE = fes.la
endif

SUBDIRS = \
    $(COND_COXETER3) cremona flint fplll gap lcalc \
    linbox $(COND_LRCALC) mpmath mwrank ntl pari singular symmetrica

include $(top_srcdir)/common.mk
libsdir = $(pythondir)/sage/libs
libs_PYTHON = __init__.py all.py

libs_polyboridir = $(libsdir)/polybori
libs_polybori_PYTHON = polybori/__init__.py

AM_DEFAULT_SOURCE_EXT = .pyx
libs_LTLIBRARIES = \
    ecl.la libsageecm.la ppl.la \
    ratpoints.la readline.la $(FES_MODULE)

fes_la_LIBADD = -lcsage -lfes

ecl_la_LIBADD = -lecl @GMP_OR_MPIR_LIBS@ -lcsage

# FIXME: use @LIBECL_INCLUDES@
ecl_la_CPPFLAGS = @LIBECL_INCLUDES@ $(AM_CPPFLAGS)
libsageecm_la_LIBADD = -lecm @GMP_OR_MPIR_LIBS@
ppl_la_SOURCES = ppl.pyxx ppl_shim.cc
ppl_la_LIBADD = -lppl @GMP_OR_MPIR_XXLIBS@ @GMP_OR_MPIR_LIBS@ -lm -lcsage
ratpoints_la_LIBADD = -lratpoints @GMP_OR_MPIR_LIBS@

libs_gmpdir = $(libsdir)/gmp
libs_gmp_PYTHON = gmp/__init__.py

noinst_HEADERS = \
    fplll/fplll.pxi coxeter3/decl.pxi pari/declinl.pxi pari/decl.pxi \
    pari/gen.pxi pari/pari_err.pxi ntl/ntl_ZZ_pEX_linkage.pxi ntl/decl.pxi \
    ntl/ntl_GF2X_linkage.pxi ntl/misc.pxi symmetrica/symmetrica.pxi \
    symmetrica/kostka.pxi symmetrica/sab.pxi symmetrica/sb.pxi \
    symmetrica/schur.pxi symmetrica/sc.pxi symmetrica/plet.pxi \
    symmetrica/part.pxi flint/nmod_poly_linkage.pxi flint/fmpz_poly.pxi \
    flint/fmpz.pxi cremona/defs.pxi singular/singular-cdefs.pxi \
    ppl_shim.hh polybori/decl.pxd
