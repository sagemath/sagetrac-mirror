# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
include $(top_srcdir)/common.mk
graphs_basedir = $(pythondir)/sage/graphs/base
graphs_base_PYTHON = __init__.py all.py graph_backends.py

AM_DEFAULT_SOURCE_EXT = .pyx
graphs_base_LTLIBRARIES = c_graph.la dense_graph.la sparse_graph.la static_sparse_graph.la

c_graph_la_LIBADD = -lcsage
dense_graph_la_LIBADD = -lcsage
sparse_graph_la_LIBADD = -lcsage
static_sparse_graph_la_LIBADD = -lcsage
