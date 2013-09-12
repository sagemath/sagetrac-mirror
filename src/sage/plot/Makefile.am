# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
SUBDIRS = plot3d
include $(top_srcdir)/common.mk
plotdir = $(pythondir)/sage/plot
plot_PYTHON = \
    __init__.py all.py animate.py arc.py arrow.py bar_chart.py bezier_path.py \
    circle.py colors.py contour_plot.py density_plot.py disk.py ellipse.py \
    graphics.py hyperbolic_arc.py hyperbolic_triangle.py java3d.py line.py \
    matrix_plot.py misc.py plot.py plot_field.py point.py polygon.py primitive.py \
    scatter_plot.py step.py text.py

AM_DEFAULT_SOURCE_EXT = .pyx
plot_LTLIBRARIES = complex_plot.la
complex_plot_la_CPPFLAGS = $(AM_CPPFLAGS) @NUMPY_INCLUDES@
complex_plot_la_LIBADD = -lcsage
