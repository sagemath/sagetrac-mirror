# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
include $(top_srcdir)/common.mk
graphs_modular_decompositiondir = $(pythondir)/sage/graphs/modular_decomposition
graphs_modular_decomposition_PYTHON = __init__.py

AM_DEFAULT_SOURCE_EXT = .pyx
graphs_modular_decomposition_LTLIBRARIES = modular_decomposition.la 
modular_decomposition_la_SOURCES = modular_decomposition.pyx src/dm.c
modular_decomposition_la_LIBADD = -lcsage
