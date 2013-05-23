# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
include $(top_srcdir)/common.mk
combinat_matricesdir = $(pythondir)/sage/combinat/matrices
combinat_matrices_PYTHON = __init__.py all.py hadamard_matrix.py dlxcpp.py latin.py

AM_DEFAULT_SOURCE_EXT = .pyxx
combinat_matrices_LTLIBRARIES = dancing_links.la
dancing_links_la_LIBADD = -lcsage
