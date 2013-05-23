# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
SUBDIRS = proof
include $(top_srcdir)/common.mk

structuredir = $(pythondir)/sage/structure
structure_PYTHON = \
    unique_representation.py sequence.py dynamic_class.py __init__.py \
    element_verify.py test_factory.py global_options.py all.py \
    factorization_integer.py factorization.py nonexact.py \
    coerce_exceptions.py formal_sum.py list_clone_timings.py gens_py.py

AM_DEFAULT_SOURCE_EXT = .pyx
structure_LTLIBRARIES = \
    category_object.la coerce.la coerce_actions.la coerce_dict.la coerce_maps.la \
    element.la factory.la generators.la list_clone.la list_clone_timings_cy.la \
    misc.la mutability.la parent.la parent_base.la parent_gens.la parent_old.la \
    sage_object.la wrapper_parent.la debug_options.la element_wrapper.la \
    list_clone_demo.la

element_la_CFLAGS = -Os $(AM_CFLAGS)

coerce_actions_la_LIBADD = -lcsage
coerce_la_LIBADD = -lcsage
coerce_maps_la_LIBADD = -lcsage
element_la_LIBADD = -lcsage
generators_la_LIBADD = -lcsage
list_clone_la_LIBADD = -lcsage
list_clone_demo_la_LIBADD = -lcsage
misc_la_LIBADD = -lcsage
parent_base_la_LIBADD = -lcsage
parent_gens_la_LIBADD = -lcsage
parent_la_LIBADD = -lcsage
parent_old_la_LIBADD = -lcsage

noinst_HEADERS = coerce.pxi debug_options.pxd
