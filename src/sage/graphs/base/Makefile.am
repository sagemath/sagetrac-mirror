# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
include $(top_srcdir)/common.mk
inst_PYTHON = __init__.py all.py graph_backends.py

LIBS+= -lcsage

AM_DEFAULT_SOURCE_EXT = .pyx
inst_LTLIBRARIES = \
    c_graph.la \
    dense_graph.la \
    sparse_graph.la \
    static_sparse_graph.la \
    static_sparse_backend.la

noinst_HEADERS = static_sparse_backend.pxd
