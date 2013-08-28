# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
include $(top_srcdir)/common.mk
functionsdir = $(pythondir)/sage/functions
functions_PYTHON = \
    trig.py orthogonal_polys.py other.py piecewise.py __init__.py transcendental.py \
    wigner.py all.py spike_function.py exp_integral.py log.py hyperbolic.py \
    special.py min_max.py generalized.py bessel.py

AM_DEFAULT_SOURCE_EXT = .pyx
functions_LTLIBRARIES = prime_pi.la
prime_pi_la_LIBADD = -lcsage
