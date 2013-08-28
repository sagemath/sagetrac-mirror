# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
include $(top_srcdir)/common.mk
quadratic_formsdir = $(pythondir)/sage/quadratic_forms
quadratic_forms_PYTHON = \
   quadratic_form__local_density_congruence.py \
   quadratic_form__genus.py quadratic_form.py quadratic_form__automorphisms.py \
   __init__.py quadratic_form__mass.py quadratic_form__variable_substitutions.py \
   quadratic_form__theta.py quadratic_form__local_normal_form.py \
   quadratic_form__count_local_2.py quadratic_form__mass__Conway_Sloane_masses.py \
   quadratic_form__equivalence_testing.py \
   quadratic_form__local_representation_conditions.py random_quadraticform.py \
   quadratic_form__mass__Siegel_densities.py special_values.py \
   quadratic_form__siegel_product.py all.py \
   quadratic_form__local_density_interfaces.py \
   quadratic_form__local_field_invariants.py constructions.py \
   quadratic_form__reduction_theory.py quadratic_form__ternary_Tornaria.py \
   extras.py binary_qf.py quadratic_form__split_local_covering.py qfsolve.py \
   quadratic_form__neighbors.py ternary_qf.py

AM_DEFAULT_SOURCE_EXT = .pyx
quadratic_forms_LTLIBRARIES = count_local_2.la quadratic_form__evaluate.la ternary.la

count_local_2_la_LIBADD = @GMP_OR_MPIR_LIBS@

quadratic_forms_generadir = $(quadratic_formsdir)/genera
quadratic_forms_genera_PYTHON = \
    genera/all.py \
    genera/genus.py \
    genera/__init__.py
