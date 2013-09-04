# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
include $(top_srcdir)/common.mk
doctestdir = $(pythondir)/sage/doctest

doctest_PYTHON = \
    all.py control.py forker.py __init__.py parsing.py reporting.py sources.py \
    test.py util.py

doctest_testsdir = $(doctestdir)/tests
doctest_tests_PYTHON = tests/nodoctest.py

EXTRA_DIST = \
    tests/1second.rst \
    tests/99seconds.rst \
    tests/abort.rst \
    tests/atexit.rst \
    tests/fail_and_die.rst \
    tests/initial.rst \
    tests/interrupt_diehard.rst \
    tests/interrupt.rst \
    tests/keyboardinterrupt.rst \
    tests/longtime.rst \
    tests/show_skipped.rst \
    tests/sig_on.rst \
    tests/simple_failure.rst \
    tests/sleep_and_raise.rst
