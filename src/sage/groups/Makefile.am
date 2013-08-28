# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
SUBDIRS = abelian_gps perm_gps matrix_gps
include $(top_srcdir)/common.mk
groupsdir = $(pythondir)/sage/groups
groups_PYTHON = \
    braid.py __init__.py class_function.py generic.py pari_group.py \
    free_group.py all.py group_homset.py groups_catalog.py conjugacy_classes.py \
    finitely_presented.py \
    libgap_group.py libgap_mixin.py

AM_DEFAULT_SOURCE_EXT = .pyx
groups_LTLIBRARIES = group.la libgap_wrapper.la old.la

groups_additive_abeliandir = $(pythondir)/sage/groups/additive_abelian
groups_additive_abelian_PYTHON = \
    additive_abelian/additive_abelian_group.py \
    additive_abelian/additive_abelian_wrapper.py \
    additive_abelian/all.py \
    additive_abelian/__init__.py

groups_misc_gpsdir = $(pythondir)/sage/groups/misc_gps
groups_misc_gps_PYTHON = \
    misc_gps/__init__.py \
    misc_gps/misc_groups_catalog.py \
    misc_gps/misc_groups.py

groups_affine_gpsdir = $(pythondir)/sage/groups/affine_gps
groups_affine_gps_PYTHON = \
    affine_gps/affine_group.py \
    affine_gps/euclidean_group.py \
    affine_gps/group_element.py \
    affine_gps/__init__.py
