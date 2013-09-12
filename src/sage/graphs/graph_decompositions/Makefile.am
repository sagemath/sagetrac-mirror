# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
include $(top_srcdir)/common.mk
inst_PYTHON = __init__.py

AM_DEFAULT_SOURCE_EXT = .pyx
inst_LTLIBRARIES = \
    graph_products.la rankwidth.la \
    vertex_separation.la

rankwidth_la_SOURCES = rankwidth.pyx rankwidth_c/rw.c
rankwidth_la_LIBADD = -lcsage
vertex_separation_la_LIBADD = -lcsage

EXTRA_DIST = fast_digraph.pyx
