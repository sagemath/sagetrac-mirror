.PHONY: all sage gp2c clean

all: sage

sage: sage/libs/pari/auto_gen.pxi sage/ext/interpreters/__init__.py
	python -u setup.py install

clean:
	@echo "Deleting Sage library build artifacts..."
	rm -rf c_lib build
	find . -name '*.pyc' | xargs rm -f
	rm -f sage/libs/pari/auto_*
	rm -rf sage/ext/interpreters


# Auto-generated files
sage/libs/pari/auto_gen.pxi: $(SAGE_LOCAL)/share/pari/pari.desc \
		sage_setup/autogen/pari/*.py $(SAGE_HAVE_GP2C) \
		sage/libs/pari/paridecl.pxd \
		sage/libs/pari/declinl.pxi \
		sage/libs/pari/scripts.pxi
	python -c "from sage_setup.autogen.pari import rebuild; rebuild()"

sage/ext/interpreters/__init__.py: sage_setup/autogen/interpreters.py
	python -c "from sage_setup.autogen.interpreters import rebuild; rebuild('sage/ext/interpreters')"


gp2c: build/gp2c/simon.c

build/gp2c/simon.c: build/gp2c/simon.gp
	gp2c -g $< -o $@

build/gp2c/simon.gp: \
		ext/pari/simon/qfsolve.gp \
		ext/pari/simon/ellQ.gp \
		ext/pari/simon/ell.gp \
		ext/pari/simon/resultant3.gp
	mkdir -p build/gp2c
	cat $^ >$@
