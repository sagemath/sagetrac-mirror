.PHONY: all sage clean

AUTOGEN_INTERPS=sage_setup/autogen/interpreters
SPEC_MODS=$(AUTOGEN_INTERPS)/specs/*.py
INTERP_PATH=sage/ext/interpreters
INTERPS=$(shell sed -n "s/\s\+name = '\(\w\+\)'/\1/p" $SPEC_MODS)
INTERP_SRCS=$(INTERP_PATH)/__init__.py \
            $(patsubst %,$(INTERP_PATH)/interp_%.c,$INTERPS) \
            $(patsubst %,$(INTERP_PATH)/wrapper_%.pxd,$INTERPS) \
            $(patsubst %,$(INTERP_PATH)/wrapper_%.pyx,$INTERPS)
INTERP_GENS=$(shell find $(AUTOGEN_INTERPS) -type f -name '*.py')

all: sage

sage: sage/libs/pari/auto_gen.pxi $(INTERP_SRCS)
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

$(filter %/__init__.py,$(INTERP_SRCS)): $(INTERP_GENS)
	# One rule just for __init__.py so the generator is run once
	# This will rebuild all the other sources though
	python -c "from sage_setup.autogen.interpreters import rebuild; rebuild('sage/ext/interpreters')"

$(filter-out %/__init__.py,$(INTERP_SRCS)): $(INTERP_PATH)/__init__.py
