# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
SUBDIRS = hmm
include $(top_srcdir)/common.mk
statsdir = $(pythondir)/sage/stats
stats_PYTHON = __init__.py all.py basic_stats.py r.py test.py

AM_DEFAULT_SOURCE_EXT = .pyx
stats_LTLIBRARIES = intlist.la

intlist_la_LIBADD = -lcsage
