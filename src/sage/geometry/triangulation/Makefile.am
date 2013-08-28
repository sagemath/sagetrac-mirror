# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
include $(top_srcdir)/common.mk
geometry_triangulationdir = $(pythondir)/sage/geometry/triangulation
geometry_triangulation_PYTHON = __init__.py all.py element.py point_configuration.py
geometry_triangulation_LTLIBRARIES = base.la
base_la_SOURCES = base.pyxx functions.cc data.cc triangulations.cc
