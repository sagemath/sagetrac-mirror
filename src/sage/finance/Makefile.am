# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
include $(top_srcdir)/common.mk
financedir = $(pythondir)/sage/finance
finance_PYTHON = __init__.py all.py easter.py markov_multifractal.py stock.py

AM_DEFAULT_SOURCE_EXT = .pyx
finance_LTLIBRARIES = \
    fractal.la markov_multifractal_cython.la time_series.la \
    option.la

time_series_la_LIBADD = -lcsage
time_series_la_CPPFLAGS = $(AM_CPPFLAGS) @NUMPY_INCLUDES@
