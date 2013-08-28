# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
SUBDIRS = notebook
include $(top_srcdir)/common.mk

serverdir = $(pythondir)/sage/server
server_PYTHON = \
    all.py __init__.py \
    introspect.py misc.py \
    nodoctest.py support.py

server_tracdir = $(serverdir)/trac
server_trac_PYTHON = \
    trac/all.py \
    trac/__init__.py \
    trac/trac.py

server_simpledir = $(serverdir)/simple
server_simple_PYTHON = \
    simple/__init__.py \
    simple/nodoctest.py \
    simple/twist.py


