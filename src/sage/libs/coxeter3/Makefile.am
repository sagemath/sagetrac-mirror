# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
include $(top_srcdir)/common.mk
libs_coxeter3dir = $(pythondir)/sage/libs/coxeter3
libs_coxeter3_PYTHON = __init__.py coxeter_group.py

AM_DEFAULT_SOURCE_EXT = .pyxx
libs_coxeter3_LTLIBRARIES = coxeter.la
coxeter_la_LIBADD = -lcsage -lcoxeter3 -lstdc++
coxeter_la_CPPFLAGS = -I@SAGE_LOCAL@/include/coxeter $(AM_CPPFLAGS)
