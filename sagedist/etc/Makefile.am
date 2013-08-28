# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#

all-local: env.sh

# this is the particular env hack
all-local: dist.sh
	echo "# this file has been placed by toplevel Make" > 00dist.sh
	echo "# edit $$(pwd)/dist.sh instead" >> 00dist.sh
	echo "#" >> 00dist.sh
	cat $< >> 00dist.sh
	$(MKDIR_P) @top_builddir@/local/etc/env.d/
	mv 00dist.sh @top_builddir@/local/etc/env.d/
