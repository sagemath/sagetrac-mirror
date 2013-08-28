# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
SUBDIRS = interpreters
include $(top_srcdir)/common.mk
extdir = $(pythondir)/sage/ext
ext_PYTHON = gen_interpreters.py __init__.py

AM_DEFAULT_SOURCE_EXT = .pyx
ext_LTLIBRARIES = \
    c_lib.la fast_callable.la fast_eval.la \
    interactive_constructors_c.la multi_modular.la \
    pselect.la

c_lib_la_LIBADD = -lcsage

multi_modular_la_LIBADD = @GMP_OR_MPIR_LIBS@ -lcsage
multi_modular_la_CFLAGS = $(AM_CFLAGS) -std=c99

fast_callable_la_LIBADD = -lcsage
fast_eval_la_LIBADD = -lcsage

pselect_la_LIBADD = -lcsage

# interpreters-stamp:
#	$(PYTHON) @srcdir@/gen_interpreters.py
# this is needed to make cython find the newly generated .pxd files
# (cython -I doesnt work otherwise)
#.PHONY: interpreters-stamp

noinst_HEADERS = \
    cdefs.pxi gmp.pxi interrupt.pxi python_debug.pxi python.pxi \
    python_rich_object.pxi python_slice.pxi random.pxi signals.pxi \
    stdsage.pxi

CLEANFILES+=interpreters-stamp interpreters/wrapper*.p* $(INIT_PY_HERE)

py: $(INIT_PY_HERE)
all-local: $(INIT_PY_HERE)
