# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
bin_SCRIPTS = sage_fortran sage-setup.py

EXTRA_DIST = sage-setup.py sage-dummy.sh sage-sh

all-local: sage-dist-make sage-dist-install

all-local: $(FOREIGN_ALL:%=%-dummy)

# sage already trusts gcc from $PATH.
# (is this a bug)?
# for now, just make sure there is no gcc dummy.
gcc-dummy:
	rm -f gcc

DUMMY_ALL=$(filter-out gcc,$(FOREIGN_ALL))
# move to configure is probably better
$(DUMMY_ALL:%=%-dummy): %-dummy:
	@if test -n "$(BINARY_$*)"; then \
	  if test x$(SAGE_$*) = xyes; then \
	    echo creating $(BINARY_$*) dummies for $*; \
	    for i in $(BINARY_$*); do \
			echo sage-dummy.sh $$i $* > $$i; \
			chmod +x $$i; \
	    done;\
	  else\
	    echo removing $(BINARY_$*) dummies;\
	    rm -f $(BINARY_$*);\
	  fi;\
	fi

.PHONY: $(FOREIGN_ALL:%=%-dummy)
