# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#

# provide old, static build interface

all: check
	./configure # creates GNUmakefile
	$(MAKE) -f GNUmakefile

check:
	@$(MAKE) --version | grep -q GNU

.PHONY: check
