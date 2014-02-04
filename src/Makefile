clean:
	@echo "Deleting Sage library build artifacts..."
	-cd c_lib && $(MAKE) distclean
	rm -rf build
	find . -name '*.pyc' | xargs rm -f
	find . -name '*.pyx' | sed 's/pyx$$/c/' | xargs rm -f
	find . -name '*.pyx' | sed 's/pyx$$/cpp/' | xargs rm -f
	rm -f sage/modular/arithgroup/farey_symbol.h
	rm -f sage/rings/real_mpfi.h
	rm -f sage/symbolic/pynac.h
