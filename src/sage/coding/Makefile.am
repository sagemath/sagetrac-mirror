# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
include $(top_srcdir)/common.mk
codingdir = $(pythondir)/sage/coding
coding_PYTHON = \
    __init__.py ag_code.py all.py code_bounds.py code_constructions.py decoder.py \
    guava.py linear_code.py sd_codes.py delsarte_bounds.py

AM_DEFAULT_SOURCE_EXT = .pyx
coding_LTLIBRARIES = binary_code.la

source_codingdir = $(codingdir)/source_coding
    source_coding_PYTHON = \
    source_coding/all.py \
    source_coding/huffman.py \
    source_coding/__init__.py
