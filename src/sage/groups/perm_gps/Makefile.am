# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
SUBDIRS = partn_ref
include $(top_srcdir)/common.mk
groups_perm_gps_PYTHON = \
    cubegroup.py __init__.py permgroup_morphism.py permgroup_named.py permgroup.py \
    all.py permutation_groups_catalog.py
groups_perm_gpsdir = $(pythondir)/sage/groups/perm_gps

AM_DEFAULT_SOURCE_EXT = .pyx
groups_perm_gps_LTLIBRARIES = permgroup_element.la

permgroup_element_la_LIBADD = -lcsage
