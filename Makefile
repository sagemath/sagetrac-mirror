# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#

# provide old, static build interface
# incomplete!

MODULES="src/sage src/c_lib src/bin src/doc src/ext"
CONFIGURE_ACS=$(MODULES:%=%/configure.ac)

all: prereq configuration
	$(MAKE) -f GNUmakefile

configuration: configure
	./configure # creates GNUmakefile

prereq:
	@$(MAKE) --version | grep -q GNU

check distclean: configure
	$(MAKE) $*

configure: $(CONFIGURE_ACS)
	./autogen.sh

.PHONY: check
