# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#

include $(top_srcdir)/common.mk
sagedir = $(pythondir)/sage
sage_matroidsdir = $(sagedir)/matroids
sage_matroids_PYTHON = \
    advanced.py \
    all.py \
    catalog.py \
    constructor.py \
    dual_matroid.py \
    __init__.py \
    matroids_catalog.py \
    minor_matroid.py \
    named_matroids.py \
    rank_matroid.py \
    utilities.py

noinst_HEADERS = \
    basis_exchange_matroid.pxd \
    basis_matroid.pxd \
    circuit_closures_matroid.pxd \
    extension.pxd \
    lean_matrix.pxd \
    linear_matroid.pxd \
    matroid.pxd \
    minorfix.h \
    set_system.pxd \
    unpickling.pxd

AM_DEFAULT_SOURCE_EXT = .pyx
sage_matroids_LTLIBRARIES = \
    matroid.la extension.la set_system.la lean_matrix.la \
    basis_exchange_matroid.la basis_matroid.la \
    linear_matroid.la circuit_closures_matroid.la \
    unpickling.la
