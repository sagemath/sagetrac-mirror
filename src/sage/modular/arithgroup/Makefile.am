# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
include $(top_srcdir)/common.mk
modular_arithgroupdir = $(pythondir)/sage/modular/arithgroup
modular_arithgroup_PYTHON = \
    congroup_generic.py __init__.py congroup_sl2z.py arithgroup_generic.py tests.py \
    congroup_gammaH.py congroup_gamma0.py congroup_gamma1.py \
    all.py arithgroup_perm.py congroup_gamma.py

AM_DEFAULT_SOURCE_EXT = .pyx
modular_arithgroup_LTLIBRARIES = congroup_pyx.la farey_symbol.la arithgroup_element.la
congroup_pyx_la_LIBADD = -lcsage
farey_symbol_la_SOURCES = farey_symbol.pyxx farey.cpp sl2z.cpp
farey_symbol_la_LIBADD = @GMP_OR_MPIR_XXLIBS@ @GMP_OR_MPIR_LIBS@ -lcsage

noinst_HEADERS = sl2z.hpp farey.hpp
