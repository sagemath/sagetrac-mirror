# Main Makefile for Sage.

# The default target ("all") builds Sage and the whole (HTML) documentation.
#
# Target "build" just builds Sage.
#
# See below for targets to build the documentation in other formats,
# to run various types of test suites, and to remove parts of the build etc.

default: all

build: all-build

# The --stop flag below is just a random flag to induce graceful
# breakage with non-GNU versions of make.
# See https://trac.sagemath.org/ticket/24617

# Defer unknown targets to build/make/Makefile
%::
	@if [ -x relocate-once.py ]; then ./relocate-once.py; fi
	$(MAKE) build/make/Makefile --stop
	+build/bin/sage-logger \
		"cd build/make && ./install '$@'" logs/install.log

# If configure was run before, rerun it with the old arguments.
# Otherwise, run configure with argument $PREREQ_OPTIONS.
build/make/Makefile: configure build/make/deps build/pkgs/*/*
	rm -f config.log
	mkdir -p logs/pkgs
	ln -s logs/pkgs/config.log config.log
	@if [ -x config.status ]; then \
		./config.status --recheck && ./config.status; \
	else \
		./configure $$PREREQ_OPTIONS; \
	fi || ( \
		if [ "x$$SAGE_PORT" = x ]; then \
			echo "If you would like to try to build Sage anyway (to help porting),"; \
			echo "export the variable 'SAGE_PORT' to something non-empty."; \
			exit 1; \
		else \
			echo "Since 'SAGE_PORT' is set, we will try to build anyway."; \
		fi; )

# Preemptively download all standard upstream source tarballs.
download:
	export SAGE_ROOT=$$(pwd) && \
	export PATH=$$SAGE_ROOT/src/bin:$$PATH && \
	./src/bin/sage-download-upstream

dist: build/make/Makefile
	./sage --sdist

# ssl: build Sage, and also install pyOpenSSL. This is necessary for
# running the secure notebook. This make target requires internet
# access. Note that this requires that your system have OpenSSL
# libraries and headers installed. See README.txt for more
# information.
ssl: all
	./sage -i pyopenssl

misc-clean:
	@echo "Deleting miscellaneous artifacts generated by build system ..."
	rm -rf logs
	rm -rf dist
	rm -rf tmp
	rm -f aclocal.m4 config.log config.status confcache
	rm -rf autom4te.cache
	rm -f build/make/Makefile build/make/Makefile-auto
	rm -f .BUILDSTART

# TODO: What is a "bdist"? A binary distribution?
bdist-clean: clean
	$(MAKE) misc-clean

distclean: build-clean
	$(MAKE) misc-clean
	@echo "Deleting all remaining output from build system ..."
	rm -rf local
	rm -f src/bin/sage-env-config

# Delete all auto-generated files which are distributed as part of the
# source tarball
bootstrap-clean:
	rm -rf config configure build/make/Makefile-auto.in

# Remove absolutely everything which isn't part of the git repo
maintainer-clean: distclean bootstrap-clean
	rm -rf upstream

micro_release: bdist-clean sagelib-clean
	@echo "Stripping binaries ..."
	LC_ALL=C find local/lib local/bin -type f -exec strip '{}' ';' 2>&1 | grep -v "File format not recognized" |  grep -v "File truncated" || true
	@echo "Removing sphinx artifacts..."
	rm -rf local/share/doc/sage/doctrees local/share/doc/sage/inventory
	@echo "Removing documentation. Inspection in IPython still works."
	rm -rf local/share/doc local/share/*/doc local/share/*/examples local/share/singular/html
	@echo "Removing unnecessary files & directories - make will not be functional afterwards anymore"
	@# We need src/sage/ for introspection with "??"
	@# We need src/sage/bin/ for the scripts that invoke Sage
	@# We need sage, the script to start Sage
	@# We need local/, the dependencies and the built Sage library itself.
	@# We keep VERSION.txt.
	@# We keep COPYING.txt so we ship a license with this distribution.
	find . -name . -o -prune ! -name src ! -name sage ! -name local ! -name VERSION.txt ! -name COPYING.txt ! -name build -exec rm -rf \{\} \;
	cd src && find . -name . -o -prune ! -name sage ! -name bin -exec rm -rf \{\} \;
	if command -v rdfind > /dev/null; then \
		@echo "Hardlinking identical files."; \
		rdfind -makeresultsfile false -makehardlinks true .; \
	else \
		@echo "rdfind not installed. Not hardlinking identical files."; \
	fi

# Leaves everything that is needed to make the next "make" fast but removes
# all the cheap build artifacts that can be quickly regenerated.
fast-rebuild-clean: misc-clean bdist-clean
	rm -rf upstream/
	rm -rf src/build/temp.*
	# Without site-packages/sage sage does not start but copying/compiling
	# them from src/build is very fast.
	rm -rf local/lib/python*/site-packages/sage
	# The .py files in src/build are restored from src/sage without their
	# mtimes changed.
	find src/build -name '*.py' -exec rm \{\} \;

TESTALL = ./sage -t --all
PTESTALL = ./sage -t -p --all

test: all
	$(TESTALL) --logfile=logs/test.log

check: test

testall: all
	$(TESTALL) --optional=sage,optional,external --logfile=logs/testall.log

testlong: all
	$(TESTALL) --long --logfile=logs/testlong.log

testalllong: all
	$(TESTALL) --long --optional=sage,optional,external --logfile=logs/testalllong.log

ptest: all
	$(PTESTALL) --logfile=logs/ptest.log

ptestall: all
	$(PTESTALL) --optional=sage,optional,external --logfile=logs/ptestall.log

ptestlong: all
	$(PTESTALL) --long --logfile=logs/ptestlong.log

ptestalllong: all
	$(PTESTALL) --long --optional=sage,optional,external --logfile=logs/ptestalllong.log

testoptional: all
	$(TESTALL) --optional=sage,optional --logfile=logs/testoptional.log

testoptionallong: all
	$(TESTALL) --long --optional=sage,optional --logfile=logs/testoptionallong.log

ptestoptional: all
	$(PTESTALL) --optional=sage,optional --logfile=logs/ptestoptional.log

ptestoptionallong: all
	$(PTESTALL) --long --optional=sage,optional --logfile=logs/ptestoptionallong.log

configure: configure.ac src/bin/sage-version.sh m4/*.m4
	./bootstrap -d

install: all
	@echo "******************************************************************"
	@echo "The '$@' target is a no-op; 'make' already does 'make install'"
	@echo "You can change the install prefix from its default"
	@echo "(the subdirectory 'local') by using ./configure --prefix=PREFIX"
	@echo "You can also consider using the binary packaging scripts"
	@echo "from https://github.com/sagemath/binary-pkg"
	@echo "******************************************************************"

.PHONY: default build install micro_release \
	misc-clean bdist-clean distclean bootstrap-clean maintainer-clean \
	test check testoptional testall testlong testoptionallong testallong \
	ptest ptestoptional ptestall ptestlong ptestoptionallong ptestallong
