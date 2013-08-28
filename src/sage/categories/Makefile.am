# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
SUBDIRS = examples
include $(top_srcdir)/common.mk
categoriesdir = $(pythondir)/sage/categories
categories_PYTHON = \
    sets_with_grading.py \
    tutorial.py ring_ideals.py modules_with_basis.py finite_enumerated_sets.py \
    gcd_domains.py monoids.py affine_weyl_groups.py commutative_algebras.py \
    additive_magmas.py left_modules.py algebras.py commutative_rings.py \
    hopf_algebras.py subobjects.py coalgebras_with_basis.py posets.py polyhedra.py \
    rngs.py finite_dimensional_algebras_with_basis.py graded_modules_with_basis.py \
    primer.py category_types.py quotient_fields.py subquotients.py tensor.py \
    finite_lattice_posets.py __init__.py modular_abelian_varieties.py pushout.py \
    group_algebras.py monoid_algebras.py graded_algebras.py isomorphic_objects.py \
    euclidean_domains.py partially_ordered_monoids.py crystals.py \
    finite_weyl_groups.py finite_posets.py semigroups.py \
    graded_bialgebras_with_basis.py classical_crystals.py graded_coalgebras.py \
    poor_man_map.py covariant_functorial_construction.py quotients.py category.py \
    coxeter_groups.py commutative_additive_groups.py pointed_sets.py sets_cat.py \
    vector_spaces.py weyl_groups.py enumerated_sets.py finite_groups.py \
    finite_permutation_groups.py graded_hopf_algebras.py \
    finite_dimensional_hopf_algebras_with_basis.py algebra_ideals.py \
    graded_coalgebras_with_basis.py finite_fields.py schemes.py finite_monoids.py \
    algebra_functor.py hopf_algebras_with_basis.py commutative_additive_monoids.py \
    highest_weight_crystals.py fields.py realizations.py division_rings.py \
    matrix_algebras.py principal_ideal_domains.py graded_algebras_with_basis.py \
    with_realizations.py graded_bialgebras.py sets_with_partial_maps.py groups.py \
    rings.py right_modules.py algebras_with_basis.py groupoid.py \
    bialgebras_with_basis.py all.py g_sets.py function_fields.py modules.py \
    magmas.py commutative_algebra_ideals.py \
    finite_dimensional_bialgebras_with_basis.py basic.py algebra_modules.py \
    number_fields.py commutative_additive_semigroups.py graded_modules.py \
    facade_sets.py semirings.py dual.py graded_hopf_algebras_with_basis.py \
    commutative_ring_ideals.py infinite_enumerated_sets.py coalgebras.py \
    bimodules.py homset.py finite_semigroups.py integral_domains.py domains.py \
    cartesian_product.py unique_factorization_domains.py finite_coxeter_groups.py \
    regular_crystals.py bialgebras.py objects.py hecke_modules.py lattice_posets.py \
    finite_crystals.py finite_dimensional_modules_with_basis.py \
    finite_dimensional_coalgebras_with_basis.py

AM_DEFAULT_SOURCE_EXT = .pyx
categories_LTLIBRARIES = action.la category_singleton.la functor.la map.la morphism.la

map_la_LIBADD = -lcsage
morphism_la_LIBADD = -lcsage
action_la_LIBADD = -lcsage
