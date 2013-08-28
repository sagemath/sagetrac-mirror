# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
include $(top_srcdir)/common.mk
schemes_toricdir = $(pythondir)/sage/schemes/toric
schemes_toric_PYTHON = \
    __init__.py all.py chow_group.py divisor.py fano_variety.py homset.py \
	 ideal.py library.py morphism.py variety.py weierstrass.py \
    weierstrass_covering.py

AM_DEFAULT_SOURCE_EXT = .pyx
schemes_toric_LTLIBRARIES = divisor_class.la
divisor_class_la_LIBADD = @GMP_OR_MPIR_LIBS@
