# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
include $(top_srcdir)/common.mk
schemes_hyperelliptic_curvesdir = $(pythondir)/sage/schemes/hyperelliptic_curves
schemes_hyperelliptic_curves_PYTHON = \
    hyperelliptic_generic.py hyperelliptic_g2_finite_field.py __init__.py \
    invariants.py jacobian_morphism.py hyperelliptic_rational_field.py \
    hyperelliptic_g2_generic.py kummer_surface.py constructor.py all.py \
    hyperelliptic_g2_padic_field.py jacobian_generic.py \
    hyperelliptic_g2_rational_field.py hyperelliptic_padic_field.py \
    jacobian_homset.py hyperelliptic_finite_field.py jacobian_g2.py

schemes_hyperelliptic_curves_LTLIBRARIES = hypellfrob.la
hypellfrob_la_SOURCES = \
	 hypellfrob.pyxx hypellfrob/hypellfrob_cc.cpp hypellfrob/recurrences_ntl.cpp \
    hypellfrob/recurrences_zn_poly.cpp

hypellfrob_la_LIBADD = -lntl -lstdc++ @GMP_OR_MPIR_LIBS@ -lzn_poly -lcsage
hypellfrob_la_CPPFLAGS = -I@srcdir@/../../libs/ntl -I@srcdir@/hypellfrob $(AM_CPPFLAGS)
