# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#

if BUILD_DOC
DOC = doc
endif

SUBDIRS = c_lib sage $(DOC) bin ext # mac-app whatever

# FIXME: environment mustnt be complicated
# FIXME: SAGE_ETC points to bin?! maybe add etc in sagelib...
check-local: # py sage
	export SAGE_ROOT=/tmp; \
	export SAGE_SRC=@abs_top_srcdir@; \
	export SAGE_EXTCODE="@abs_top_builddir@/ext"; \
	export LIBCSAGE_INCLUDES="-I@abs_top_builddir@/c_lib/includes -I@abs_top_srcdir@/c_lib/includes"; \
	export PYTHONPATH="@abs_top_builddir@$(PYTHONPATH:%=:%)"; \
	export PATH="@abs_top_builddir@/bin:@abs_top_srcdir@/bin:$$PATH"; \
	export SAGE_MAXIMA_INIT=@abs_top_builddir@/bin/sage-maxima.lisp; \
	export SAGE_ETC=@abs_top_builddir@/bin; \
	./bin/sage -t --all

py:
	$(MAKE) -C sage py

EXTRA_DIST = autogen.sh install-sh
