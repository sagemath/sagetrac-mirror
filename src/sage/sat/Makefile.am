# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
SUBDIRS = solvers
include $(top_srcdir)/common.mk
satdir = $(pythondir)/sage/sat
sat_PYTHON = \
    boolean_polynomials.py \
    __init__.py

satconvertersdir = $(satdir)/converters
satconverters_PYTHON = \
    converters/anf2cnf.py \
    converters/__init__.py \
    converters/polybori.py
