# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
include $(top_srcdir)/common.mk
rings_padicsdir = $(pythondir)/sage/rings/padics
rings_padics_PYTHON = \
    all.py eisenstein_extension_generic.py factory.py generic_nodes.py __init__.py \
    local_generic.py misc.py padic_base_generic.py padic_base_leaves.py \
    padic_extension_generic.py padic_extension_leaves.py padic_generic.py \
    precision_error.py tests.py tutorial.py unramified_extension_generic.py

AM_DEFAULT_SOURCE_EXT = .pyx
rings_padics_LTLIBRARIES = \
    local_generic_element.la padic_ZZ_pX_CA_element.la \
    padic_ZZ_pX_CR_element.la padic_ZZ_pX_FM_element.la padic_ZZ_pX_element.la \
    padic_base_coercion.la padic_base_generic_element.la \
    padic_capped_absolute_element.la padic_capped_relative_element.la \
    padic_ext_element.la padic_fixed_mod_element.la padic_generic_element.la \
    padic_printing.la pow_computer.la pow_computer_ext.la

padic_ZZ_pX_CA_element_la_SOURCES = padic_ZZ_pX_CA_element.pyxx
padic_ZZ_pX_CA_element_la_LIBADD = -lntl @GMP_OR_MPIR_LIBS@ -lcsage @GMP_OR_MPIR_XXLIBS@ -lm -lstdc++

padic_ZZ_pX_CR_element_la_SOURCES = padic_ZZ_pX_CR_element.pyxx
padic_ZZ_pX_CR_element_la_LIBADD = -lntl @GMP_OR_MPIR_LIBS@ -lcsage @GMP_OR_MPIR_XXLIBS@ -lm -lstdc++

padic_ZZ_pX_FM_element_la_SOURCES = padic_ZZ_pX_FM_element.pyxx
padic_ZZ_pX_FM_element_la_LIBADD = -lntl @GMP_OR_MPIR_LIBS@ -lcsage @GMP_OR_MPIR_XXLIBS@ -lm -lstdc++

padic_ZZ_pX_element_la_SOURCES = padic_ZZ_pX_element.pyxx
padic_ZZ_pX_element_la_LIBADD = -lntl @GMP_OR_MPIR_LIBS@ -lcsage @GMP_OR_MPIR_XXLIBS@ -lm -lstdc++

padic_base_coercion_la_LIBADD = @GMP_OR_MPIR_LIBS@ -lcsage

padic_base_generic_element_la_LIBADD = @GMP_OR_MPIR_LIBS@
padic_capped_absolute_element_la_LIBADD = @GMP_OR_MPIR_LIBS@ -lcsage
padic_capped_relative_element_la_LIBADD = @GMP_OR_MPIR_LIBS@ -lcsage

padic_ext_element_la_SOURCES = padic_ext_element.pyxx
padic_ext_element_la_LIBADD = -lntl @GMP_OR_MPIR_LIBS@ -lcsage @GMP_OR_MPIR_XXLIBS@ -lm -lstdc++

padic_fixed_mod_element_la_LIBADD = @GMP_OR_MPIR_LIBS@ -lcsage
padic_generic_element_la_LIBADD = @GMP_OR_MPIR_LIBS@ -lcsage

padic_printing_la_LIBADD = @GMP_OR_MPIR_LIBS@ -lntl -lcsage @GMP_OR_MPIR_XXLIBS@ -lm -lstdc++
padic_printing_la_SOURCES = padic_printing.pyxx

pow_computer_la_LIBADD = -lcsage -lntl @GMP_OR_MPIR_LIBS@ @GMP_OR_MPIR_XXLIBS@ -lm
pow_computer_la_SOURCES = pow_computer.pyxx

pow_computer_ext_la_LIBADD = -lcsage -lntl @GMP_OR_MPIR_LIBS@ @GMP_OR_MPIR_XXLIBS@ -lm
pow_computer_ext_la_SOURCES = pow_computer_ext.pyxx
