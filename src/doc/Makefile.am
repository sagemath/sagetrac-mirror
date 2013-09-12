# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#

SUBDIRS = common

install-data-local: docbuild

all: docbuild

EXTRA_DIST = en/reference/footer.txt \
             en/reference/index.rst \
             en/reference/conf_sub.py \
             en/reference/conf.py \
             $(shell (cd @srcdir@; ls en/reference/*/index.rst)) \
             $(shell (cd @srcdir@; ls en/reference/*/conf.py)) \
             $(shell (cd @srcdir@; ls en/reference/combinat/*.rst)) \
             $(shell (cd @srcdir@; ls en/reference/polynomial_rings/polynomial_rings_*.rst)) \
             en/reference/misc/sagetex.rst \
             en/reference/cmd/environ.rst \
             en/reference/cmd/options.rst \
             en/reference/cmd/startup.rst \
             en/reference/dynamics/interval_exchanges.rst \
             en/reference/dynamics/flat_surfaces.rst \
             en/reference/todolist.rst

LANGDIRS = de en fr ru tr
LANGSUBDIRS = $(shell (cd @srcdir@; for i in $(LANGDIRS); do ls $$i/* -d; done))
EXTRA_DIST+= $(filter-out %/reference, $(LANGSUBDIRS))

# hmm better build env?

docbuild: common/builder.py common/build_options.py
	for i in $(LANGDIRS); do \
	    $(MKDIR_P) $$i/reference; \
	    ( cd $$i/reference; [ -L footer.txt -o -f footer.txt ] ||\
	        $(LN_S) @abs_top_srcdir@/$$i/reference/footer.txt . ); \
	done
	export SAGE_DOC="@abs_top_builddir@"; \
	export SAGE_DOC_SRC="@abs_top_srcdir@"; \
	export SAGE_SRC="@abs_top_srcdir@/.."; \
	export DOT_SAGE="$$( mktemp -d )"; \
	export PATH="@abs_top_srcdir@/../bin:$(PATH)"; \
	export SAGE_MAXIMA_INIT="@abs_top_builddir@/../bin/sage-maxima.lisp"; \
	export LD_LIBRARY_PATH="@abs_top_builddir@/../c_lib/src/.libs$(LD_LIBRARY_PATH:%=:)$(LD_LIBRARY_PATH)"; \
	export PYTHONPATH=common:$$PYTHONPATH$(PYTHONPATH:%=:)@abs_top_builddir@/..:@abs_top_srcdir@; \
	$(PYTHON) $< --no-pdf-links all html; ret=$$?; \
	rm -r $$DOT_SAGE; \
	exit $$ret;

install-data-local:
	$(INSTALL) -d $(docdir)
	cp -r output/* $(docdir)

clean-local:
	@echo "Deleting generated docs..."
	for i in $(LANGDIRS); do \
	    ( cd $$i/reference; [ ! -L footer.txt ] || rm -f footer.txt ) \
	done
	rm -rf en/reference/*/sage
	rm -rf en/reference/*/sagenb
	rm -rf en/reference/sage
	rm -rf en/reference/sagenb
	rm -rf output

dist-hook:
	find $(LANGDIRS:%=$(distdir)/%) -name *~ -delete

.PHONY: docbuild
