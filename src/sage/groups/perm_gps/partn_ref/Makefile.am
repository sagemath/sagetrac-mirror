# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
include $(top_srcdir)/common.mk
groups_perm_gps_partn_refdir = $(pythondir)/sage/groups/perm_gps/partn_ref
groups_perm_gps_partn_ref_PYTHON = __init__.py

AM_DEFAULT_SOURCE_EXT = .pyx
groups_perm_gps_partn_ref_LTLIBRARIES = \
    automorphism_group_canonical_label.la canonical_augmentation.la double_coset.la \
    refinement_binary.la refinement_graphs.la refinement_lists.la \
    refinement_matrices.la refinement_python.la refinement_sets.la

# all of them (correct?)
LIBS += -lcsage

automorphism_group_canonical_label_la_LIBADD = @GMP_OR_MPIR_LIBS@ -lflint
automorphism_group_canonical_label_la_CFLAGS = -std=c99 $(AM_CFLAGS)
canonical_augmentation_la_LIBADD = @GMP_OR_MPIR_LIBS@ -lflint
canonical_augmentation_la_CFLAGS = -std=c99 $(AM_CFLAGS)
double_coset_la_LIBADD = @GMP_OR_MPIR_LIBS@ -lflint
double_coset_la_CFLAGS = -std=c99 $(AM_CFLAGS)
refinement_binary_la_LIBADD = @GMP_OR_MPIR_LIBS@ -lflint
refinement_binary_la_CFLAGS = -std=c99 $(AM_CFLAGS)
refinement_graphs_la_LIBADD = @GMP_OR_MPIR_LIBS@ -lflint
refinement_graphs_la_CFLAGS = -std=c99 $(AM_CFLAGS)
refinement_lists_la_LIBADD = @GMP_OR_MPIR_LIBS@ -lflint
refinement_lists_la_CFLAGS = -std=c99 $(AM_CFLAGS)
refinement_matrices_la_LIBADD = @GMP_OR_MPIR_LIBS@ -lflint
refinement_matrices_la_CFLAGS = -std=c99 $(AM_CFLAGS)
refinement_python_la_LIBADD = @GMP_OR_MPIR_LIBS@ -lflint
refinement_python_la_CFLAGS = -std=c99 $(AM_CFLAGS)
refinement_sets_la_LIBADD = @GMP_OR_MPIR_LIBS@ -lflint
refinement_sets_la_CFLAGS = -std=c99 $(AM_CFLAGS)

noinst_HEADERS = data_structures_pxd.pxi data_structures_pyx.pxi
