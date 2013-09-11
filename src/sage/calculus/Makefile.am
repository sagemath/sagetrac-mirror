# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
include $(top_srcdir)/common.mk
calculusdir = $(pythondir)/sage/calculus
calculus_PYTHON = \
    wester.py __init__.py tests.py functional.py all.py test_sympy.py functions.py \
    predefined.py desolvers.py calculus.py

AM_DEFAULT_SOURCE_EXT = .pyx
calculus_LTLIBRARIES = interpolators.la riemann.la var.la

interpolators_la_CPPFLAGS = @NUMPY_INCLUDES@ $(AM_CPPFLAGS)

riemann_la_CPPFLAGS = @NUMPY_INCLUDES@ $(AM_CPPFLAGS)
riemann_la_LIBADD = -lcsage
