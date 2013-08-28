# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
libs_singulardir = $(pythondir)/sage/libs/singular
include $(top_srcdir)/common.mk
libs_singular_PYTHON = __init__.py function_factory.py

AM_DEFAULT_SOURCE_EXT = .pyxx
libs_singular_LTLIBRARIES = function.la groebner_strategy.la option.la polynomial.la ring.la singular.la 
function_la_LIBADD = -lm -lreadline -lsingular -lgivaro -lntl @GMP_OR_MPIR_XXLIBS@ @GMP_OR_MPIR_LIBS@ -lcsage
function_la_CPPFLAGS = $(AM_CPPFLAGS) @SINGULAR_INCLUDES@ \
                       -D__STDC_LIMIT_MACROS

groebner_strategy_la_LIBADD = -lm -lreadline -lsingular -lgivaro -lntl @GMP_OR_MPIR_XXLIBS@ @GMP_OR_MPIR_LIBS@
groebner_strategy_la_CXXFLAGS = -D__STDC_LIMIT_MACROS $(AM_CXXFLAGS)
groebner_strategy_la_CPPFLAGS = @SINGULAR_INCLUDES@ $(AM_CPPFLAGS)
option_la_LIBADD = -lm -lreadline -lsingular -lgivaro -lntl @GMP_OR_MPIR_XXLIBS@ @GMP_OR_MPIR_LIBS@
option_la_CXXFLAGS = -D__STDC_LIMIT_MACROS $(AM_CXXFLAGS)
option_la_CPPFLAGS = @SINGULAR_INCLUDES@ $(AM_CPPFLAGS)
polynomial_la_LIBADD = -lm -lreadline -lsingular -lgivaro -lntl @GMP_OR_MPIR_XXLIBS@ @GMP_OR_MPIR_LIBS@ -lcsage
polynomial_la_CPPFLAGS = $(AM_CPPFLAGS) @SINGULAR_INCLUDES@ \
                         -D__STDC_LIMIT_MACROS

ring_la_LIBADD = -lm -lreadline -lsingular -lgivaro -lntl @GMP_OR_MPIR_XXLIBS@ @GMP_OR_MPIR_LIBS@ -lcsage
ring_la_CPPFLAGS = $(AM_CPPFLAGS) @SINGULAR_INCLUDES@ \
                   -D__STDC_LIMIT_MACROS

singular_la_LIBADD = -lm -lreadline -lsingular -lgivaro -lntl @GMP_OR_MPIR_XXLIBS@ @GMP_OR_MPIR_LIBS@ -lcsage
singular_la_CXXFLAGS = -D__STDC_LIMIT_MACROS $(AM_CXXFLAGS)
singular_la_CPPFLAGS = @SINGULAR_INCLUDES@ $(AM_CPPFLAGS)

# FIXME: something extremely wrong
function_la_CPPFLAGS += -DNDEBUG
polynomial_la_CPPFLAGS += -DNDEBUG
ring_la_CPPFLAGS += -DNDEBUG
groebner_strategy_la_CPPFLAGS += -DNDEBUG
