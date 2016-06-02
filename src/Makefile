.PHONY: all sage clean

# Build a list of sources that should be generated for the
# sage.ext.interpreters package
INTERP_PATH:=sage/ext/interpreters
INTERP_SOURCES:=$(shell python -c 'from sage_setup.autogen.interpreters import *; list_sources("$(INTERP_PATH)")')

all: sage

sage: sage/libs/pari/auto_gen.pxi $(INTERP_SOURCES)
	python -u setup.py install

clean:
	@echo "Deleting Sage library build artifacts..."
	rm -rf c_lib build
	find . -name '*.pyc' | xargs rm -f
	rm -f sage/libs/pari/auto_*
	rm -rf sage/ext/interpreters

# Auto-generated files
sage/libs/pari/auto_gen.pxi: $(SAGE_LOCAL)/share/pari/pari.desc \
        sage/libs/pari/decl.pxi sage_setup/autogen/pari/*.py
	python -c "from sage_setup.autogen.pari import rebuild; rebuild()"


AUTOGEN_INTERP_SOURCES:=$(shell find sage_setup/autogen/interpreters -type f -name '*.py')
# One rule just for __init__.py so the generator is run once
# This will rebuild all the other sources though
$(filter %/__init__.py,$(INTERP_SOURCES)): $(AUTOGEN_INTERP_SOURCES)
	python -c "from sage_setup.autogen.interpreters import rebuild; rebuild('sage/ext/interpreters')"

$(filter-out %/__init__.py,$(INTERP_SOURCES)): $(INTERP_PATH)/__init__.py
