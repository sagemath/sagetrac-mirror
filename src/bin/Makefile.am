# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#

ACLOCAL_AMFLAGS = -I m4

dist_bin_SCRIPTS = \
    sage \
    sage-build \
    sage-CSI \
    sage-CSI-helper.py \
    sage-cleaner \
    sage-combinat \
    sage-coverageall \
    sage-cython \
    sage-download-file \
    sage-download-upstream \
    sage-eval \
    sage-fixdoctests \
    sage-gdb-commands \
    sage-grep \
    sage-grepdoc \
    sage-inline-fortran \
    sage-ipython \
    sage-location \
    sage-notebook \
    sage-num-threads.py \
    sage-native-execute \
    sage-omega \
    sage-open \
    sage-preparse \
    sage-python \
    sage-run \
    sage-run-cython \
    sage-runtests \
    sage-starts \
    sage-startuptime.py \
    sage-test-import

sageconfdir=$(sysconfdir)
sageconf_DATA = \
    env.sh \
    sage-banner \
    sage-started.txt \
    sage-current-location.txt

sage-started.txt:
	touch $@

# probably doesnt matter much for installed instances of sagelib
sage-current-location.txt:
	echo $(prefix) > $@

#FIXME: this does not exactly belong here
EXTRA_DIST = sage-banner

#FIXME: ext(code)?
sageconf_DATA+= sage-maxima.lisp
EXTRA_DIST+= sage-maxima.lisp

# fixme: create symlink?
# sage-bdist
# sage-cachegrind
# sage-callgrind
# sage-coverage
# sage-coverageall
# sage-crap
# sage-fix-pkg-checksums
# sage-list-experimental
# sage-list-optional
# sage-list-packages
# sage-list-standard
# sage-massif
# sage-notebook-insecure
# sage-pkg
# sage-pypkg-location
# sage-README-osx.txt
# sage-real-upgrade
# sage-rebaseall.bat
# sage-rebaseall.sh
# sage-rebase.bat
# sage-rebase.sh
# sage-rst2sws
# sage-rst2txt
# sage-rsyncdist
# sage-sdist
# sage-startuptime.py
# sage-sync-build.py
# sage-upgrade
# sage-valgrind
