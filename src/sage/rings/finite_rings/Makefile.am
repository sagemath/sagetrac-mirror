# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
include $(top_srcdir)/common.mk
rings_finite_ringsdir = $(pythondir)/sage/rings/finite_rings
rings_finite_rings_PYTHON = \
    all.py constructor.py element_ext_pari.py finite_field_ext_pari.py \
    finite_field_givaro.py finite_field_ntl_gf2e.py finite_field_prime_modn.py \
    homset.py __init__.py integer_mod_ring.py finite_field_pari_ffelt.py

AM_DEFAULT_SOURCE_EXT = .pyx
rings_finite_rings_LTLIBRARIES = \
    element_base.la element_givaro.la element_ntl_gf2e.la finite_field_base.la \
    integer_mod.la element_pari_ffelt.la

element_base_la_CPPFLAGS = -DINTEGER_MOD_INT32_LIMIT=46341 -DINTEGER_MOD_INT64_LIMIT=2147483647 $(AM_CPPFLAGS)
element_base_la_LIBADD = -lcsage

element_givaro_la_SOURCES = element_givaro.pyxx
element_givaro_la_LIBADD = -lgivaro -lntl @GMP_OR_MPIR_XXLIBS@ @GMP_OR_MPIR_LIBS@ -lm -lstdc++ -lcsage

element_givaro_la_CPPFLAGS = -D__STDC_LIMIT_MACROS \
                             -DINTEGER_MOD_INT32_LIMIT=46341 \
                             -DINTEGER_MOD_INT64_LIMIT=2147483647 \
                             $(AM_CPPFLAGS)

element_ntl_gf2e_la_SOURCES = element_ntl_gf2e.pyxx
element_ntl_gf2e_la_LIBADD = -lntl @GMP_OR_MPIR_LIBS@ -lcsage

element_ntl_gf2e_la_CPPFLAGS = $(AM_CPPFLAGS) \
                               -DINTEGER_MOD_INT32_LIMIT=46341 \
                               -DINTEGER_MOD_INT64_LIMIT=2147483647

element_pari_ffelt_la_LIBADD = -lcsage

finite_field_base_la_CPPFLAGS = $(AM_CPPFLAGS) -DINTEGER_MOD_INT32_LIMIT=46341 -DINTEGER_MOD_INT64_LIMIT=2147483647 $(AM_CPPFLAGS)
finite_field_base_la_LIBADD = -lcsage

integer_mod_la_LIBADD = @GMP_OR_MPIR_LIBS@ -lcsage
integer_mod_la_CPPFLAGS = -DINTEGER_MOD_INT32_LIMIT=46341 -DINTEGER_MOD_INT64_LIMIT=2147483647 $(AM_CPPFLAGS)

noinst_HEADERS = stdint.pxd element_pari_ffelt.pxd
