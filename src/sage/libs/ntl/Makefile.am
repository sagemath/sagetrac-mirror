# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
include $(top_srcdir)/common.mk
libs_ntldir = $(pythondir)/sage/libs/ntl
libs_ntl_PYTHON = __init__.py all.py

AM_DEFAULT_SOURCE_EXT = .pyxx
libs_ntl_LTLIBRARIES = \
    ntl_GF2.la ntl_GF2E.la ntl_GF2EContext.la ntl_GF2EX.la ntl_GF2X.la ntl_ZZ.la \
    ntl_ZZX.la ntl_ZZ_p.la ntl_ZZ_pContext.la ntl_ZZ_pE.la ntl_ZZ_pEContext.la \
    ntl_ZZ_pEX.la ntl_ZZ_pX.la ntl_lzz_p.la ntl_lzz_pContext.la ntl_lzz_pX.la \
    ntl_mat_GF2.la ntl_mat_GF2E.la ntl_mat_ZZ.la

ntl_GF2_la_LIBADD = -lcsage -lntl @GMP_OR_MPIR_LIBS@ @GMP_OR_MPIR_XXLIBS@ -lstdc++
ntl_GF2E_la_LIBADD = -lcsage -lntl @GMP_OR_MPIR_LIBS@ @GMP_OR_MPIR_XXLIBS@ -lm -lstdc++
ntl_GF2EContext_la_LIBADD = -lcsage -lntl @GMP_OR_MPIR_LIBS@ @GMP_OR_MPIR_XXLIBS@ -lm -lstdc++
ntl_GF2EX_la_LIBADD = -lcsage -lntl @GMP_OR_MPIR_LIBS@ @GMP_OR_MPIR_XXLIBS@ -lm -lstdc++
ntl_GF2X_la_LIBADD = -lcsage -lntl @GMP_OR_MPIR_LIBS@ @GMP_OR_MPIR_XXLIBS@ -lm -lstdc++
ntl_ZZ_la_LIBADD = -lcsage -lntl @GMP_OR_MPIR_LIBS@ @GMP_OR_MPIR_XXLIBS@ -lm -lstdc++
ntl_ZZX_la_LIBADD = -lcsage -lntl @GMP_OR_MPIR_LIBS@ @GMP_OR_MPIR_XXLIBS@ -lm -lstdc++
ntl_ZZ_p_la_LIBADD = -lcsage -lntl @GMP_OR_MPIR_LIBS@ @GMP_OR_MPIR_XXLIBS@ -lm -lstdc++
ntl_ZZ_pContext_la_LIBADD = -lcsage -lntl @GMP_OR_MPIR_LIBS@ @GMP_OR_MPIR_XXLIBS@ -lm -lstdc++
ntl_ZZ_pE_la_LIBADD = -lcsage -lntl @GMP_OR_MPIR_LIBS@ @GMP_OR_MPIR_XXLIBS@ -lm -lstdc++
ntl_ZZ_pEContext_la_LIBADD = -lcsage -lntl @GMP_OR_MPIR_LIBS@ @GMP_OR_MPIR_XXLIBS@ -lm -lstdc++
ntl_ZZ_pEX_la_LIBADD = -lcsage -lntl @GMP_OR_MPIR_LIBS@ @GMP_OR_MPIR_XXLIBS@ -lm -lstdc++
ntl_ZZ_pX_la_LIBADD = -lcsage -lntl @GMP_OR_MPIR_LIBS@ @GMP_OR_MPIR_XXLIBS@ -lm -lstdc++
ntl_lzz_p_la_LIBADD = -lcsage -lntl @GMP_OR_MPIR_LIBS@ @GMP_OR_MPIR_XXLIBS@ -lm -lstdc++
ntl_lzz_pContext_la_LIBADD = -lcsage -lntl @GMP_OR_MPIR_LIBS@ @GMP_OR_MPIR_XXLIBS@ -lm -lstdc++
ntl_lzz_pX_la_LIBADD = -lcsage -lntl @GMP_OR_MPIR_LIBS@ @GMP_OR_MPIR_XXLIBS@ -lm -lstdc++
ntl_mat_GF2_la_LIBADD = -lcsage -lntl @GMP_OR_MPIR_LIBS@ @GMP_OR_MPIR_XXLIBS@ -lm -lstdc++
ntl_mat_GF2E_la_LIBADD = -lcsage -lntl @GMP_OR_MPIR_LIBS@ @GMP_OR_MPIR_XXLIBS@ -lm -lstdc++
ntl_mat_ZZ_la_LIBADD = -lcsage -lntl @GMP_OR_MPIR_LIBS@ @GMP_OR_MPIR_XXLIBS@ -lm -lstdc++
