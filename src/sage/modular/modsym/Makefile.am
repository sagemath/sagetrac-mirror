# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
include $(top_srcdir)/common.mk
modular_modsymdir = $(pythondir)/sage/modular/modsym
modular_modsym_PYTHON = \
    space.py ghlist.py __init__.py manin_symbols.py boundary.py tests.py modsym.py \
    p1list_nf.py ambient.py relation_matrix.py element.py all.py hecke_operator.py \
    subspace.py modular_symbols.py g1list.py

AM_DEFAULT_SOURCE_EXT = .pyx
modular_modsym_LTLIBRARIES = apply.la heilbronn.la p1list.la relation_matrix_pyx.la
apply_la_LIBADD = -lcsage -lflint @GMP_OR_MPIR_LIBS@ @GMP_OR_MPIR_XXLIBS@ -lm -lstdc++
apply_la_CFLAGS = -std=c99 -D_XPG6 $(AM_CFLAGS)
heilbronn_la_LIBADD = -lcsage -lflint @GMP_OR_MPIR_LIBS@ @GMP_OR_MPIR_XXLIBS@ -lm -lstdc++
heilbronn_la_CFLAGS = -std=c99 -D_XPG6 $(AM_CFLAGS)
p1list_la_LIBADD = @GMP_OR_MPIR_LIBS@ -lcsage
