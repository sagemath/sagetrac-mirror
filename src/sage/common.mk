# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
# this file helps working around automake limitations
# eventually it can be removed, after
# - all default values are moved to configure.ac
# - automake supports cython
#   - dependencies
#   - vpath quirks
# - automake supports AM_RECURSIVE_TARGETS
#

# this is stupid, intended to approximate distutils.
# see http://www.gnu.org/software/automake/manual/html_node/Wildcards.html
#
#allpythonstuff = $(wildcard *.py *.pxd *.pxi)

# compile within the environment specified during configure
maybe_colon = $(my_LD_LIBRARY_PATH:%=:)
LIBTOOL = export LD_LIBRARY_PATH="$(my_LD_LIBRARY_PATH)$(LD_LIBRARY_PATH:%=$(maybe_colon)%)"; @LIBTOOL@

# -L@top_builddir@/../c_lib/src/.libs will be added automatically
AM_LDFLAGS = -module -avoid-version

AM_CPPFLAGS = @my_CPPFLAGS@

# how to do this right? LDADD?
AM_LDFLAGS+= -L@top_builddir@/../c_lib/src/.libs

# user specified MAKEFLAGS
# AM_CPPFLAGS += $(CPPFLAGS)

# just append globally, needed quite often.
AM_CPPFLAGS += @LIBCSAGE_INCLUDES@

# this is the directory that contains config.h
AM_CPPFLAGS += -I$(top_builddir)

# this is the directory that contains "sage"
# (might change later)
AM_CPPFLAGS += -I$(top_builddir)/.. -I$(top_srcdir)/..

# setup.py defines these
AM_CPPFLAGS += -I$(top_builddir)/ext -I$(top_srcdir)/ext

# BUG. that's what sage-upstream does
AM_CFLAGS = -fno-strict-aliasing -fwrapv
AM_CXXFLAGS = -fno-strict-aliasing -fwrapv

# FIXME: place where required (partly done)
LIBS += -lcsage

# dont spam terminal with error messages
CYTHONFLAGS = --fast-fail

# write cython_debug
# does not work. workdir must be @srcdir@
# CYTHON_GDBOPT = --gdb
CYTHON_GDBOPT = --gdb-outdir @abs_top_builddir@/..
if CYGDB
CYTHONFLAGS += $(CYTHON_GDBOPT)
endif

# this is required for inspection
CYTHONFLAGS += --embed-positions

# sage still relies on ...
CYTHONFLAGS += --old-style-globals

# more switches
CYTHONFLAGS += -X autotestdict=False
CYTHONFLAGS += -X cdivision=True
CYTHONFLAGS += -X fast_getattr=True

# problem: include paths...
CYTHONFLAGS += -I$(abs_top_srcdir)/.. -I$(abs_top_builddir)/..

# should look like this (but doesn't, because we are still in src/sage
# CYTHONFLAGS += -I$(abs_top_srcdir) -I$(abs_top_builddir)

CYTHON ?= cython

# simplify paths
here_rel = $(subst $(abs_top_srcdir),,$(abs_srcdir))
instdir = $(pythondir)/sage$(here_rel)

# does not work with subdir sources
install-data-hook: $(filter %.pyx,$(SOURCES))
	@echo installing more stuff
	@for i in $<; do \
		$(INSTALL) -t $(instdir) $$i; \
	done

# need to be adapted in during core modules build system merge
sagelib_abs_top_srcdir = $(abs_top_srcdir)/..
sagelib_abs_top_builddir = $(abs_top_builddir)/..

AM_V_CYT = $(am__v_CYT_$(V))
am__v_CYT_ = $(am__v_CYT_$(AM_DEFAULT_VERBOSITY))
am__v_CYT_0 = @echo "  CYTH    " $@;

AM_V_PYC = $(am__v_PYC_$(V))
am__v_PYC_ = $(am__v_PYC_$(AM_DEFAULT_VERBOSITY))
am__v_PYC_0 = @echo "  PYC     " $@;

AM_V_PYO = $(am__v_PYO_$(V))
am__v_PYO_ = $(am__v_PYO_$(AM_DEFAULT_VERBOSITY))
am__v_PYO_0 = @echo "  PYO     " $@;

# ouch, trailing colon.
# empty string will be (mis)interpreted as '.'.
PYTHONPATHENV = PYTHONPATH="@abs_top_builddir@/..$(PYTHONPATH:%=:%)"

# -w needed to force correct paths in docstrings
# wrong paths lead to mysterious sageinspect.py errors

define cython_call
	$(AM_V_CYT)$(PYTHONPATHENV) $(CYTHON) $(CYTHONFLAGS) \
	    -w $(if $(findstring @srcdir@,$<),$(sagelib_abs_top_srcdir),$(sagelib_abs_top_builddir)) \
	    $(abspath $<) \
	    -o $(abs_builddir)/$@ @AMDEP_TRUE@@PYDEP_TRUE@-MD -MP
	$(AM_V_at)@AMDEP_TRUE@@PYDEP_TRUE@mv $@.d $(DEPDIR)/$*.Pcython
endef

%.cc: CYTHONFLAGS+= --cplus

MANUAL_DEP_PYS = $(filter %.pyx %.pyxx,$(SOURCES))
MANUAL_DEP_PYX = $(MANUAL_DEP_PYS:%.pyxx=%.pyx)
MANUAL_DEP = $(MANUAL_DEP_PYX:%.pyx=$(DEPDIR)/%.Pcython)

CLEANFILES = $(MANUAL_DEP_PYX:%.pyx=%.c) \
             $(MANUAL_DEP_PYX:%.pyx=%.cc) \
             *.so *.pyc *.pyo

# BUG: this makes config.status create \*.Pcython ...
@AMDEP_TRUE@ifneq (,$(MANUAL_DEP))
@AMDEP_TRUE@@am__include@ $(DEPDIR)/*.Pcython
@AMDEP_TRUE@endif

.pyx.c:
.pyxx.cc:

.SECONDEXPANSION:
%.c: %.pyx $$(if $$(wildcard $$*.pxd),$$(wildcard $$*.pxd), \
	              $$(if $$(wildcard @srcdir@/$$*.pxd), \
                       $$(wildcard @srcdir@/$$*.pxd)))
	$(cython_call)

.SECONDEXPANSION:
%.cc: %.pyxx $$(if $$(wildcard $$*.pxd),$$(wildcard $$*.pxd), \
	                $$(if $$(wildcard @srcdir@/$$*.pxd), \
                         $$(wildcard @srcdir@/$$*.pxd)))
	$(cython_call)

PYS = $(filter %.py,$(DIST_COMMON))
PYCS = $(PYS:%.py=%.pyc)
PYOS = $(PYS:%.py=%.pyo)

py-local: $(LTLIBRARIES:%.la=%.so) $(PYCS) $(PYOS)
	
# FIXME: V
$(LTLIBRARIES:%.la=%.so): %.so: | %.la
	$(AM_V_at)[ -f $@ ] || $(LN_S) .libs/$@ .

# manually implementing AM_EXTRA_RECURSIVE_TARGETS([py pycheck])
# will be implemented in automake1.12
# (pycheck is not a good name...)
py: py-local py-recursive
py-recursive:
	for i in $(SUBDIRS); do $(MAKE) -C $$i py || exit 1; done
pycheck: pycheck-local pycheck-recursive
pycheck-recursive:
	for i in $(SUBDIRS); do $(MAKE) -C $$i pycheck PYLIST=../$(PYLIST); done

@am__leading_dot@PHONY: py-recursive py-local py \
                        pycheck-recursive pycheck-local pycheck

# create python module for checking
check-local: py-local

PYLIST = none

# compare registered .py's against existing
pycheck-local:
	echo $(PYS) $(noinst_HEADERS) $(filter %.pyx,$(DIST_SOURCES)) $(filter %.pyxx,$(DIST_SOURCES)) | tr ' ' '\n' | \
	    sed 's#^#@abs_builddir@/#' | \
	    sed 's#^@abs_top_builddir@#.#' >> $(PYLIST)

# this won't work as long as _PYTHON includes subdirectories
# @diff <(echo "$(PYS)" | tr ' ' '\n' | sort) <(cd $(VPATH); ls *.py | sort)

#don't delete .cc .c just because gcc fails.
@am__leading_dot@PRECIOUS: %.cc %.c

if VPATH_BUILD
# need some path trickery, because inspection does not work otherwise
%.pyc: SHELL=/usr/bin/env bash
%.pyo: SHELL=/usr/bin/env bash
%.pyc: %.py
	@$(MKDIR_P) $(dir $@)
	$(AM_V_PYC)echo -e "\n__file__='$<'" | cat "$<" - | $(PYTHON) -c 'import py_compile; py_compile.compile("/dev/stdin","$@","$<")'
%.pyo: %.py
	@$(MKDIR_P) $(dir $@)
	$(AM_V_PYO)echo -e "\n__file__='$<'" | cat "$<" - | $(PYTHON) -O -c 'import py_compile; py_compile.compile("/dev/stdin","$@","$<")'
else
# use plain python bytecompiler.
%.pyc: %.py
	$(AM_V_PYC)$(PYTHON) -c 'import py_compile; py_compile.compile("$<","$@","$<")'
%.pyo: %.py
	$(AM_V_PYO)$(PYTHON) -O -c 'import py_compile; py_compile.compile("$<","$@","$<")'
endif

# sometimes, __init__.py is required to make
# cython -I work.
if VPATH_BUILD
INIT_PY_HERE = @abs_builddir@/__init__.py
$(INIT_PY_HERE): @srcdir@/__init__.py
	[ -f __init__.py -o -L __init__.py ] ||\
	    $(LN_S) $< .
else
INIT_PY_HERE =
endif
