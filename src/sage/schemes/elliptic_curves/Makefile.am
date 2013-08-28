# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
include $(top_srcdir)/common.mk
schemes_elliptic_curvesdir = $(pythondir)/sage/schemes/elliptic_curves
schemes_elliptic_curves_PYTHON = \
    gal_reps.py mod5family.py ell_point.py ell_rational_field.py padics.py \
    __init__.py ell_local_data.py ell_finite_field.py ell_padic_field.py \
    ell_generic.py ell_wp.py sha_tate.py heegner.py lseries_ell.py ec_database.py \
    ell_modular_symbols.py formal_group.py ell_tate_curve.py ell_field.py \
    ell_torsion.py weierstrass_morphism.py kodaira_symbol.py constructor.py all.py \
    ell_curve_isogeny.py cm.py padic_lseries.py gp_simon.py monsky_washnitzer.py \
    BSD.py ell_number_field.py modular_parametrization.py ell_egros.py \
    isogeny_class.py period_lattice.py \
    gal_reps_number_field.py weierstrass_transform.py jacobian.py

AM_DEFAULT_SOURCE_EXT = .pyx
schemes_elliptic_curves_LTLIBRARIES = descent_two_isogeny.la
descent_two_isogeny_la_LIBADD = -lflint @GMP_OR_MPIR_LIBS@ -lratpoints
descent_two_isogeny_la_CFLAGS = -std=c99 $(AM_CFLAGS)
