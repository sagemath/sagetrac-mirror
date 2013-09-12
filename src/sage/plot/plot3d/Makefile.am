# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
include $(top_srcdir)/common.mk
plot_plot3ddir = $(pythondir)/sage/plot/plot3d
plot_plot3d_PYTHON = \
    __init__.py all.py examples.py implicit_plot3d.py list_plot3d.py \
    parametric_plot3d.py platonic.py plot3d.py plot_field3d.py revolution_plot3d.py \
    shapes2.py tachyon.py texture.py tri_plot.py

AM_DEFAULT_SOURCE_EXT = .pyx
plot_plot3d_LTLIBRARIES = \
    base.la implicit_surface.la index_face_set.la parametric_surface.la shapes.la \
    transform.la

base_la_CFLAGS = -std=c99 $(AM_CFLAGS)
implicit_surface_la_LIBADD = -lgsl -lcsage
implicit_surface_la_CPPFLAGS = $(AM_CPPFLAGS) @NUMPY_INCLUDES@
index_face_set_la_CFLAGS = -std=c99 $(AM_CFLAGS)
index_face_set_la_LIBADD = -lcsage
parametric_surface_la_LIBADD = -lcsage

noinst_HEADERS = point_c.pxi
