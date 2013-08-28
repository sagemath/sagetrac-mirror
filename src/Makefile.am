# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#

if BUILD_DOC
DOC=doc
endif

SUBDIRS = c_lib sage $(DOC) bin ext # mac-app whatever

# FIXME: environment mustnt be complicated
check-local: # py sage
	export SAGE_ROOT=/tmp; \
	export SAGE_SRC=@abs_top_srcdir@; \
	export SAGE_EXTCODE="@abs_top_builddir@/ext"; \
	export PYTHONPATH="@abs_top_builddir@:$$PYTHONPATH"; \
	export PATH="@abs_top_builddir@/bin:@abs_top_srcdir@/bin:$$PATH"; \
	export SAGE_MAXIMA_INIT=/etc/sage/sage-maxima.lisp; \
	export SAGE_ETC=@abs_top_builddir@/bin; \
	./bin/sage -t --all

py:
	$(MAKE) -C sage py

EXTRA_DIST = autogen.sh install-sh
